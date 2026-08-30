"""Двовимірна відбивна ТВЗ ВВЕР-1000 UGD із 331 позицією решітки.

Картограма містить 300 твелів UO2, 12 твегів UO2-Gd2O3, 18 заповнених водою
напрямних труб і одну центральну трубу з водою. Відбивний аксіальний шар 1 см та
відбивна межа ТВЗ задають нескінченну решітку нескінченно високих ТВЗ. Усі
геометричні, матеріальні та розрахункові параметри читаються з JSON.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent


SCORE_UK = {
    "flux": "потік",
    "absorption": "поглинання",
    "fission": "поділ",
    "nu-fission": "нейтрони поділу",
    "kappa-fission": "енергія поділів",
    "(n,gamma)": "радіаційне захоплення",
}
TYPE_UK = {
    "fuel": "звичайний твел UO₂",
    "gd": "твег UO₂-Gd₂O₃",
    "guide": "напрямна труба",
    "central": "центральна труба",
}
COLUMN_UK = {
    "assembly_id": "номер_твз",
    "assembly_name": "назва_твз",
    "assembly_type": "тип_твз",
    "local_x_cm": "локальна_x_см",
    "local_y_cm": "локальна_y_см",
    "power": "потужність",
    "power_std_dev_approx": "похибка_потужності_наближено",
    "position_id": "номер_позиції",
    "fuel_cell_id": "ідентифікатор_паливної_комірки",
    "axial_x": "аксіальна_x",
    "axial_alpha": "аксіальна_альфа",
    "x_cm": "x_см",
    "y_cm": "y_см",
    "type": "тип",
    "score": "показник",
    "mean": "середнє",
    "std_dev": "стандартне_відхилення",
    "material": "матеріал",
    "nuclide": "нуклід",
    "energy_low_eV": "нижня_межа_енергії_еВ",
    "energy_high_eV": "верхня_межа_енергії_еВ",
    "flux_fraction": "частка_потоку",
    "flux_fraction_std_dev_approx": "похибка_частки_потоку_наближено",
    "kappa_fission": "енергія_поділів_еВ_на_джерельну_частинку",
    "kappa_fission_std_dev": "похибка_енергії_поділів_еВ",
    "fission": "поділи_на_джерельну_частинку",
    "fission_std_dev": "похибка_поділів",
    "flux": "потік_на_джерельну_частинку",
    "flux_std_dev": "похибка_потоку",
    "absorption": "поглинання_на_джерельну_частинку",
    "absorption_std_dev": "похибка_поглинання",
    "nu_fission": "нейтрони_поділу_на_джерельну_частинку",
    "nu_fission_std_dev": "похибка_нейтронів_поділу",
    "power_fraction": "частка_потужності",
    "relative_power": "відносна_потужність",
    "relative_power_std_dev_approx": "похибка_відносної_потужності_наближено",
}


def ukrainian_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Повернути копію таблиці з українськими заголовками та категоріями."""
    result = dataframe.copy()
    if "score" in result:
        result["score"] = result["score"].replace(SCORE_UK)
    if "type" in result:
        result["type"] = result["type"].replace(TYPE_UK)
    return result.rename(columns=COLUMN_UK)


def ukrainian_summary(summary: dict) -> dict:
    """Підготувати україномовний підсумок без зміни внутрішніх ключів моделі."""
    maximum = summary["maximum_power_pin"]
    gd_absorption = {
        isotope: {
            "середнє_на_джерельну_частинку": values["mean_per_source"],
            "похибка_на_джерельну_частинку_наближено": values["std_dev_per_source_approx"],
        }
        for isotope, values in summary["gd_absorption"].items()
    }
    energy_flux = [
        {
            "нижня_межа_енергії_еВ": group["energy_low_eV"],
            "верхня_межа_енергії_еВ": group["energy_high_eV"],
            "середній_потік": group["mean"],
            "похибка_потоку": group["std_dev"],
            "частка_потоку": group["flux_fraction"],
            "похибка_частки_потоку_наближено": group["flux_fraction_std_dev_approx"],
        }
        for group in summary["energy_group_flux"]
    ]
    return {
        "k_inf": summary["keff"],
        "похибка_k_inf": summary["keff_std_dev"],
        "очікувана_частка_витоку": summary["leakage_expected"],
        "кількість_позицій": summary["number_of_positions"],
        "кількість_паливних_твелів": summary["number_of_fuel_pins"],
        "кількість_звичайних_твелів_uo2": summary["number_of_uo2_pins"],
        "кількість_твегів": summary["number_of_gd_pins"],
        "найбільш_навантажений_твел": {
            "номер_позиції": maximum["position_id"],
            "тип": TYPE_UK.get(maximum["type"], maximum["type"]),
            "аксіальна_x": maximum["axial_x"],
            "аксіальна_альфа": maximum["axial_alpha"],
            "x_см": maximum["x_cm"],
            "y_см": maximum["y_cm"],
            "відносна_потужність": maximum["relative_power"],
            "похибка_відносної_потужності_наближено": maximum["relative_power_std_dev_approx"],
        },
        "K_pin": summary["K_pin"],
        "похибка_K_pin_наближено": summary["K_pin_std_dev_approx"],
        "середня_відносна_потужність_звичайних_твелів": summary["ordinary_pin_mean_relative_power"],
        "похибка_середньої_потужності_звичайних_твелів_наближено": summary["ordinary_pin_mean_relative_power_std_dev_approx"],
        "середня_відносна_потужність_твегів": summary["gd_pin_mean_relative_power"],
        "похибка_середньої_потужності_твегів_наближено": summary["gd_pin_mean_relative_power_std_dev_approx"],
        "коефіцієнт_варіації_потужності_твелів": summary["pin_power_coefficient_of_variation"],
        "відношення_захоплень_u238_до_поділів_u235": summary["conversion_ratio_U238_capture_over_U235_fission"],
        "похибка_відношення_захоплень_наближено": summary["conversion_ratio_std_dev_approx"],
        "відношення_нейтронів_поділу_до_поглинань": summary["nu_fission_over_absorption"],
        "похибка_відношення_нейтронів_поділу_наближено": summary["nu_fission_over_absorption_std_dev_approx"],
        "поглинання_гадолінієм": gd_absorption,
        "потік_за_енергетичними_групами": energy_flux,
        "примітка_щодо_похибок": "Похибки відношень і середніх не враховують коваріацію між інтервалами лічильника.",
    }


def validate_parameters(p: dict) -> dict:
    """Перевірити структуру, числові межі та геометричну сумісність ТВЗ."""
    if not isinstance(p, dict):
        raise ValueError("Кореневий елемент JSON має бути об'єктом")

    required = {
        "fuel": [
            "uo2_u235_enrichment_wt_percent", "ugd_u235_enrichment_wt_percent",
            "gd2o3_wt_percent", "temperature_K", "uo2_density_g_cm3",
            "gd2o3_density_g_cm3",
        ],
        "moderator": ["temperature_K", "density_g_cm3", "boron_ppm"],
        "cladding": ["temperature_K"],
        "helium": ["temperature_K", "density_g_cm3"],
        "geometry": [
            "pin_pitch_cm", "assembly_flat_to_flat_cm", "fuel_hole_radius_cm",
            "fuel_radius_cm", "clad_inner_radius_cm", "clad_outer_radius_cm",
            "guide_inner_radius_cm", "guide_outer_radius_cm",
            "central_inner_radius_cm", "central_outer_radius_cm", "z_min_cm",
            "z_max_cm", "radial_boundary_type", "axial_boundary_type",
        ],
        "layout": [
            "hex_radius", "central_coordinate", "gd_orbit_seeds", "guide_orbit_seeds",
        ],
        "simulation": ["particles", "batches", "inactive", "threads", "seed"],
    }
    missing_sections = [section for section in required if section not in p]
    if missing_sections:
        raise ValueError(f"Відсутні секції параметрів: {', '.join(missing_sections)}")
    for section, keys in required.items():
        if not isinstance(p[section], dict):
            raise ValueError(f"Секція {section} має бути JSON-об'єктом")
        missing = [key for key in keys if key not in p[section]]
        if missing:
            raise ValueError(f"Відсутні параметри {section}: {', '.join(missing)}")

    def finite(section: str, key: str) -> float:
        try:
            value = float(p[section][key])
        except (TypeError, ValueError) as error:
            raise ValueError(f"{section}.{key} має бути числом") from error
        if not math.isfinite(value):
            raise ValueError(f"{section}.{key} має бути скінченним числом")
        return value

    def positive_density(section: str, key: str, upper: float) -> float:
        value = finite(section, key)
        if not 0.0 < value <= upper:
            raise ValueError(f"{section}.{key} має бути в інтервалі (0, {upper}] г/см³")
        return value

    def positive_integer(section: str, key: str) -> int:
        value = p[section][key]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"{section}.{key} має бути додатним цілим числом")
        return value

    for key in ("uo2_u235_enrichment_wt_percent", "ugd_u235_enrichment_wt_percent"):
        enrichment = finite("fuel", key)
        if not 0.0 < enrichment <= 20.0:
            raise ValueError(f"fuel.{key} має бути в інтервалі (0, 20] мас.%")
    gd_fraction = finite("fuel", "gd2o3_wt_percent")
    if not 0.0 < gd_fraction < 100.0:
        raise ValueError("fuel.gd2o3_wt_percent має бути в інтервалі (0, 100) мас.%")
    positive_density("fuel", "uo2_density_g_cm3", 11.5)
    positive_density("fuel", "gd2o3_density_g_cm3", 10.0)
    positive_density("moderator", "density_g_cm3", 1.5)
    positive_density("helium", "density_g_cm3", 0.1)
    boron_ppm = finite("moderator", "boron_ppm")
    if not 0.0 <= boron_ppm <= 10000.0:
        raise ValueError("moderator.boron_ppm має бути в інтервалі [0, 10000] ppm")

    for section in ("fuel", "moderator", "cladding", "helium"):
        temperature = finite(section, "temperature_K")
        if not 294.0 <= temperature <= 1200.0:
            raise ValueError(f"{section}.temperature_K має бути в інтервалі [294, 1200] K")

    g = p["geometry"]
    pitch = finite("geometry", "pin_pitch_cm")
    flat_to_flat = finite("geometry", "assembly_flat_to_flat_cm")
    hole = finite("geometry", "fuel_hole_radius_cm")
    fuel_radius = finite("geometry", "fuel_radius_cm")
    clad_inner = finite("geometry", "clad_inner_radius_cm")
    clad_outer = finite("geometry", "clad_outer_radius_cm")
    guide_inner = finite("geometry", "guide_inner_radius_cm")
    guide_outer = finite("geometry", "guide_outer_radius_cm")
    central_inner = finite("geometry", "central_inner_radius_cm")
    central_outer = finite("geometry", "central_outer_radius_cm")
    z_min = finite("geometry", "z_min_cm")
    z_max = finite("geometry", "z_max_cm")
    if pitch <= 0.0 or flat_to_flat <= 0.0:
        raise ValueError("Крок решітки та розмір ТВЗ мають бути додатними")
    if not 0.0 <= hole < fuel_radius < clad_inner < clad_outer < 0.5 * pitch:
        raise ValueError("Радіуси твела не вміщуються в одну позицію решітки")
    if not 0.0 < guide_inner < guide_outer < 0.5 * pitch:
        raise ValueError("Радіуси напрямної труби не вміщуються в одну позицію решітки")
    if not 0.0 < central_inner < central_outer < 0.5 * pitch:
        raise ValueError("Радіуси центральної труби не вміщуються в одну позицію решітки")
    if z_min >= z_max:
        raise ValueError("geometry.z_min_cm має бути меншим за geometry.z_max_cm")
    allowed_boundaries = {"vacuum", "reflective", "periodic", "white"}
    for key in ("radial_boundary_type", "axial_boundary_type"):
        if g[key] not in allowed_boundaries:
            raise ValueError(
                f"geometry.{key} має бути одним із: {', '.join(sorted(allowed_boundaries))}"
            )

    layout = p["layout"]
    radius = layout["hex_radius"]
    if isinstance(radius, bool) or not isinstance(radius, int) or radius != 10:
        raise ValueError("Повний гексагон на 331 позицію вимагає layout.hex_radius = 10")
    central = layout["central_coordinate"]
    if not (isinstance(central, (list, tuple)) and len(central) == 2
            and all(isinstance(value, int) and not isinstance(value, bool) for value in central)):
        raise ValueError("layout.central_coordinate має містити два цілі числа")
    if tuple(central) != (0, 0):
        raise ValueError("Центральна вимірювальна труба має бути в координаті [0, 0]")
    for key in ("gd_orbit_seeds", "guide_orbit_seeds"):
        seeds = layout[key]
        if not isinstance(seeds, list) or not seeds:
            raise ValueError(f"layout.{key} має бути непорожнім списком координат")
        for seed in seeds:
            if not (isinstance(seed, (list, tuple)) and len(seed) == 2
                    and all(isinstance(value, int) and not isinstance(value, bool) for value in seed)):
                raise ValueError(f"Кожен seed у layout.{key} має містити два цілі числа")

    largest_outer_radius = max(clad_outer, guide_outer, central_outer)
    required_flat_to_flat = math.sqrt(3.0) * radius * pitch + 2.0 * largest_outer_radius
    if flat_to_flat + 1.0e-9 < required_flat_to_flat:
        raise ValueError(
            "geometry.assembly_flat_to_flat_cm замалий: для заданої решітки та труб "
            f"потрібно щонайменше {required_flat_to_flat:.6f} см"
        )

    sim = p["simulation"]
    particles = positive_integer("simulation", "particles")
    batches = positive_integer("simulation", "batches")
    threads = positive_integer("simulation", "threads")
    positive_integer("simulation", "seed")
    inactive = sim["inactive"]
    if isinstance(inactive, bool) or not isinstance(inactive, int) or not 0 <= inactive < batches:
        raise ValueError("simulation.inactive має бути цілим числом у межах [0, batches)")
    if batches - inactive < 30:
        raise ValueError("Потрібно щонайменше 30 активних пакетів для оцінки статистики")
    if threads > particles:
        raise ValueError("simulation.threads не може перевищувати simulation.particles")

    edges = p.get("energy_groups_eV")
    if not isinstance(edges, list) or len(edges) < 2:
        raise ValueError("energy_groups_eV має містити щонайменше дві межі")
    try:
        numeric_edges = [float(value) for value in edges]
    except (TypeError, ValueError) as error:
        raise ValueError("Усі межі energy_groups_eV мають бути числами") from error
    if any(not math.isfinite(value) or value < 0.0 for value in numeric_edges):
        raise ValueError("Межі energy_groups_eV мають бути скінченними та невід'ємними")
    if any(high <= low for low, high in zip(numeric_edges, numeric_edges[1:])):
        raise ValueError("Межі energy_groups_eV мають строго зростати")
    if "run_name" in p and (not isinstance(p["run_name"], str) or not p["run_name"].strip()):
        raise ValueError("run_name має бути непорожнім рядком")

    # Ця перевірка також ловить перекриття орбіт, координати за межами ТВЗ
    # та неправильні кількості 300/12/18/1.
    layout_records(p)
    return p


def load_parameters(path: Path) -> dict:
    """Прочитати JSON та виконати повну перевірку параметрів ТВЗ."""
    with path.open(encoding="utf-8") as stream:
        return validate_parameters(json.load(stream))


def rotate60(coord: tuple[int, int]) -> tuple[int, int]:
    x, a = coord
    return -a, x + a


def orbit(seed: list[int] | tuple[int, int]) -> set[tuple[int, int]]:
    current = (int(seed[0]), int(seed[1]))
    result = set()
    for _ in range(6):
        result.add(current)
        current = rotate60(current)
    if len(result) != 6:
        raise ValueError(f"Orbit seed {seed} does not generate six distinct positions")
    return result


def hex_positions(radius: int) -> list[tuple[int, int]]:
    return [
        (x, a)
        for x in range(-radius, radius + 1)
        for a in range(-radius, radius + 1)
        if max(abs(x), abs(a), abs(-x - a)) <= radius
    ]


def lattice_xy(
    coord: tuple[int, int], pitch: float, orientation: str = "y"
) -> tuple[float, float]:
    """Декартові координати центра HexLattice для заданої орієнтації OpenMC."""
    x, a = coord
    if orientation == "y":
        return math.sqrt(0.75) * pitch * x, pitch * (0.5 * x + a)
    if orientation == "x":
        # Поворот orientation='y' на -90°: перший елемент кільця розташований
        # праворуч, а наступні йдуть за годинниковою стрілкою, як очікує OpenMC.
        return pitch * (0.5 * x + a), -math.sqrt(0.75) * pitch * x
    raise ValueError("orientation має бути 'x' або 'y'")


def ring_index(coord: tuple[int, int], number_of_rings: int) -> tuple[int, int]:
    """Перетворити логічні координати OpenMC на індекси кілець від краю до центра."""
    x, a = coord
    z = -a - x
    radius = max(abs(x), abs(a), abs(z))
    ring = number_of_rings - 1 - radius
    if radius == 0:
        return ring, 0
    if x >= 0:
        within = x if a >= 0 else 2 * radius + z
    else:
        within = 3 * radius - x if a <= 0 else 5 * radius - z
    return ring, within


def layout_records(p: dict) -> list[dict]:
    layout = p["layout"]
    pitch = float(p["geometry"]["pin_pitch_cm"])
    radius = int(layout["hex_radius"])
    central = tuple(layout["central_coordinate"])
    gd = set().union(*(orbit(seed) for seed in layout["gd_orbit_seeds"]))
    guides = set().union(*(orbit(seed) for seed in layout["guide_orbit_seeds"]))
    if gd & guides or central in gd or central in guides:
        raise ValueError("Central, Gd, and guide coordinates must not overlap")

    all_coords = set(hex_positions(radius))
    if not gd <= all_coords or not guides <= all_coords or central not in all_coords:
        raise ValueError("At least one special coordinate lies outside the 331-position hexagon")
    if len(gd) != 12 or len(guides) != 18 or len(all_coords) != 331:
        raise ValueError(
            f"Expected 331 positions/12 Gd/18 guides, got {len(all_coords)}/{len(gd)}/{len(guides)}"
        )

    # Зручна нумерація: рядки згори вниз, у кожному рядку зліва направо.
    ordered = sorted(all_coords, key=lambda c: (-lattice_xy(c, pitch)[1], lattice_xy(c, pitch)[0]))
    records = []
    for position_id, coord in enumerate(ordered, start=1):
        if coord == central:
            kind = "central"
        elif coord in guides:
            kind = "guide"
        elif coord in gd:
            kind = "gd"
        else:
            kind = "fuel"
        x_cm, y_cm = lattice_xy(coord, pitch)
        records.append({
            "position_id": position_id,
            "axial_x": coord[0],
            "axial_alpha": coord[1],
            "x_cm": x_cm,
            "y_cm": y_cm,
            "type": kind,
        })
    counts = pd.Series([r["type"] for r in records]).value_counts().to_dict()
    expected = {"fuel": 300, "gd": 12, "guide": 18, "central": 1}
    if counts != expected:
        raise ValueError(f"Incorrect layout counts: {counts}; expected {expected}")
    return records


def make_materials(p: dict) -> dict[str, openmc.Material]:
    f = p["fuel"]
    uo2 = openmc.Material(name="UO2 fuel")
    uo2.add_element("U", 1.0, enrichment=f["uo2_u235_enrichment_wt_percent"])
    uo2.add_element("O", 2.0)
    uo2.set_density("g/cm3", f["uo2_density_g_cm3"])
    uo2.temperature = f["temperature_K"]

    ugd_uo2 = openmc.Material(name="UO2 component for U-Gd fuel")
    ugd_uo2.add_element("U", 1.0, enrichment=f["ugd_u235_enrichment_wt_percent"])
    ugd_uo2.add_element("O", 2.0)
    ugd_uo2.set_density("g/cm3", f["uo2_density_g_cm3"])
    ugd_uo2.temperature = f["temperature_K"]
    gd2o3 = openmc.Material(name="Natural Gd2O3 component")
    gd2o3.add_element("Gd", 2.0)
    gd2o3.add_element("O", 3.0)
    gd2o3.set_density("g/cm3", f["gd2o3_density_g_cm3"])
    gd2o3.temperature = f["temperature_K"]
    gd_fraction = f["gd2o3_wt_percent"] / 100.0
    ugd = openmc.Material.mix_materials(
        [ugd_uo2, gd2o3], [1.0 - gd_fraction, gd_fraction], percent_type="wo"
    )
    ugd.name = "UO2-Gd2O3 fuel"
    ugd.temperature = f["temperature_K"]

    helium = openmc.Material(name="Helium")
    helium.add_nuclide("He4", 1.0)
    helium.set_density("g/cm3", p["helium"]["density_g_cm3"])
    helium.temperature = p["helium"]["temperature_K"]

    e110 = openmc.Material(name="E110 Zr-Nb alloy")
    e110.add_element("Zr", 4.259e-2)
    e110.add_nuclide("Nb93", 4.225e-4)
    e110.add_element("Hf", 6.597e-6)
    e110.set_density("sum")
    e110.temperature = p["cladding"]["temperature_K"]

    mod = p["moderator"]
    water = openmc.model.borated_water(
        boron_ppm=mod["boron_ppm"],
        temperature=mod["temperature_K"],
        density=mod["density_g_cm3"],
        name="Water moderator",
    )
    water.temperature = mod["temperature_K"]
    return {"uo2": uo2, "ugd": ugd, "helium": helium, "e110": e110, "water": water}


def fuel_universe(record: dict, fuel: openmc.Material, mat: dict, p: dict) -> tuple[openmc.Universe, openmc.Cell]:
    g = p["geometry"]
    hole = openmc.ZCylinder(r=g["fuel_hole_radius_cm"])
    pellet = openmc.ZCylinder(r=g["fuel_radius_cm"])
    clad_i = openmc.ZCylinder(r=g["clad_inner_radius_cm"])
    clad_o = openmc.ZCylinder(r=g["clad_outer_radius_cm"])
    tag = f"position {record['position_id']}"
    fuel_cell = openmc.Cell(name=f"fuel {tag}", fill=fuel, region=+hole & -pellet)
    cells = [
        openmc.Cell(name=f"central helium {tag}", fill=mat["helium"], region=-hole),
        fuel_cell,
        openmc.Cell(name=f"helium gap {tag}", fill=mat["helium"], region=+pellet & -clad_i),
        openmc.Cell(name=f"cladding {tag}", fill=mat["e110"], region=+clad_i & -clad_o),
        openmc.Cell(name=f"outer water {tag}", fill=mat["water"], region=+clad_o),
    ]
    return openmc.Universe(name=f"{record['type']} universe {tag}", cells=cells), fuel_cell


def tube_universe(record: dict, tube_type: str, mat: dict, p: dict) -> openmc.Universe:
    g = p["geometry"]
    inner_key = "guide_inner_radius_cm" if tube_type == "guide" else "central_inner_radius_cm"
    outer_key = "guide_outer_radius_cm" if tube_type == "guide" else "central_outer_radius_cm"
    inner = openmc.ZCylinder(r=g[inner_key])
    outer = openmc.ZCylinder(r=g[outer_key])
    tag = f"{tube_type} position {record['position_id']}"
    return openmc.Universe(name=tag, cells=[
        openmc.Cell(name=f"inner water {tag}", fill=mat["water"], region=-inner),
        openmc.Cell(name=f"tube wall {tag}", fill=mat["e110"], region=+inner & -outer),
        openmc.Cell(name=f"outer water {tag}", fill=mat["water"], region=+outer),
    ])


def build_model(p: dict) -> tuple[openmc.Model, dict]:
    mat = make_materials(p)
    records = layout_records(p)
    radius = p["layout"]["hex_radius"]
    num_rings = radius + 1
    rings: list[list[openmc.Universe | None]] = [
        [None] * (6 * (radius - ring_index_) if ring_index_ < radius else 1)
        for ring_index_ in range(num_rings)
    ]
    fuel_cells = []
    fuel_records = []
    for record in records:
        if record["type"] in ("fuel", "gd"):
            universe, fuel_cell = fuel_universe(
                record, mat["uo2"] if record["type"] == "fuel" else mat["ugd"], mat, p
            )
            record["fuel_cell_id"] = fuel_cell.id
            fuel_cells.append(fuel_cell)
            fuel_records.append(record)
        else:
            universe = tube_universe(record, record["type"], mat, p)
        ring, within = ring_index((record["axial_x"], record["axial_alpha"]), num_rings)
        rings[ring][within] = universe
    if any(universe is None for ring in rings for universe in ring):
        raise RuntimeError("Not all HexLattice positions were filled")

    outer_water = openmc.Universe(
        name="outer lattice water",
        cells=[openmc.Cell(fill=mat["water"], name="outer lattice water cell")],
    )
    lattice = openmc.HexLattice(name="VVER-1000 331-position lattice")
    lattice.center = (0.0, 0.0)
    lattice.pitch = (p["geometry"]["pin_pitch_cm"],)
    lattice.orientation = "y"
    lattice.universes = rings
    lattice.outer = outer_water

    g = p["geometry"]
    boundary = openmc.model.HexagonalPrism(
        edge_length=g["assembly_flat_to_flat_cm"] / math.sqrt(3.0),
        orientation="y",
        boundary_type=g["radial_boundary_type"],
    )
    z_low = openmc.ZPlane(z0=g["z_min_cm"], boundary_type=g["axial_boundary_type"])
    z_high = openmc.ZPlane(z0=g["z_max_cm"], boundary_type=g["axial_boundary_type"])
    root_cell = openmc.Cell(name="assembly", fill=lattice, region=-boundary & +z_low & -z_high)
    geometry = openmc.Geometry(openmc.Universe(cells=[root_cell]))

    sim = p["simulation"]
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = sim["particles"]
    settings.batches = sim["batches"]
    settings.inactive = sim["inactive"]
    settings.seed = sim["seed"]
    settings.temperature = {"method": "interpolation", "range": (294.0, 1200.0)}
    half = 0.5 * g["assembly_flat_to_flat_cm"]
    source_box = openmc.stats.Box(
        [-half, -half, g["z_min_cm"]], [half, half, g["z_max_cm"]]
    )
    settings.source = openmc.IndependentSource(
        space=source_box, constraints={"fissionable": True}
    )

    pin_tally = openmc.Tally(name="pin_results")
    pin_tally.filters = [openmc.CellFilter(fuel_cells)]
    pin_tally.scores = ["kappa-fission", "fission", "flux", "absorption", "nu-fission"]

    assembly_tally = openmc.Tally(name="assembly_rates")
    assembly_tally.scores = [
        "flux", "absorption", "fission", "nu-fission", "kappa-fission", "(n,gamma)"
    ]

    spectrum = openmc.Tally(name="energy_group_flux")
    spectrum.filters = [openmc.EnergyFilter(p["energy_groups_eV"])]
    spectrum.scores = ["flux"]

    isotope = openmc.Tally(name="fuel_isotope_rates")
    isotope.filters = [openmc.MaterialFilter([mat["uo2"], mat["ugd"]])]
    isotope.nuclides = ["U235", "U238", "Gd155", "Gd157"]
    isotope.scores = ["absorption", "fission", "(n,gamma)"]

    plot = openmc.Plot(name="assembly_geometry")
    plot.filename = "heometriia_tvz_331"
    plot.basis = "xy"
    plot.origin = (0.0, 0.0, 0.0)
    plot.width = (25.0, 28.0)
    plot.pixels = (1000, 1120)
    plot.color_by = "material"
    plot.colors = {
        mat["uo2"]: "gold",
        mat["ugd"]: "forestgreen",
        mat["helium"]: "white",
        mat["e110"]: "slategray",
        mat["water"]: "deepskyblue",
    }

    model = openmc.Model(
        geometry=geometry,
        materials=openmc.Materials(list(mat.values())),
        settings=settings,
        tallies=openmc.Tallies([pin_tally, assembly_tally, spectrum, isotope]),
        plots=openmc.Plots([plot]),
    )
    return model, {
        "records": records,
        "fuel_records": fuel_records,
        "fuel_cells": fuel_cells,
        "materials": mat,
    }


def scalar_rows(tally: openmc.Tally) -> list[dict]:
    return [
        {
            "score": score,
            "mean": float(tally.mean[0, 0, index]),
            "std_dev": float(tally.std_dev[0, 0, index]),
        }
        for index, score in enumerate(tally.scores)
    ]


def ratio_with_error(a: float, sa: float, b: float, sb: float) -> tuple[float, float]:
    if b == 0.0:
        raise ZeroDivisionError("Неможливо обчислити відношення: знаменник tally дорівнює нулю")
    if a == 0.0:
        # Лінійне поширення похибки для r=a/b у точці a=0.
        return 0.0, abs(sa / b)
    ratio = a / b
    # OpenMC не зберігає коваріацію між цими bins tally.
    error = abs(ratio) * math.sqrt((sa / a) ** 2 + (sb / b) ** 2)
    return ratio, error


def mean_with_independent_error(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    return float(values.mean()), float(np.sqrt(np.square(errors).sum()) / len(values))


def make_power_map(df: pd.DataFrame, all_records: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 13.5))
    fuel = df[df["type"].isin(["fuel", "gd"])]
    vmax = max(1.15, float(fuel["relative_power"].max()))
    points = ax.scatter(
        fuel["x_cm"], fuel["y_cm"], c=fuel["relative_power"], cmap="turbo",
        vmin=0.0, vmax=vmax, marker="h", s=310, edgecolors="0.35", linewidths=0.45,
    )
    gd = fuel[fuel["type"] == "gd"]
    ax.scatter(gd["x_cm"], gd["y_cm"], facecolors="none", edgecolors="lime",
               marker="h", s=390, linewidths=2.0, label="Твег UO₂-Gd₂O₃")
    nonfuel = pd.DataFrame([r for r in all_records if r["type"] in ("guide", "central")])
    guides = nonfuel[nonfuel["type"] == "guide"]
    central = nonfuel[nonfuel["type"] == "central"]
    ax.scatter(guides["x_cm"], guides["y_cm"], c="#7c3aed", marker="h", s=310,
               edgecolors="#3b0764", linewidths=1.0, label="Напрямна труба")
    ax.scatter(central["x_cm"], central["y_cm"], c="#dc2626", marker="h", s=330,
               edgecolors="#7f1d1d", linewidths=1.2, label="Центральна труба")
    for row in fuel.itertuples():
        ax.text(row.x_cm, row.y_cm, str(row.position_id), ha="center", va="center",
                fontsize=4.0, color="black")
    max_pin = fuel.loc[fuel["relative_power"].idxmax()]
    ax.scatter([max_pin["x_cm"]], [max_pin["y_cm"]], facecolors="none", edgecolors="black",
               marker="h", s=470, linewidths=2.5)
    offset_x = -18 if max_pin["x_cm"] > 0 else 18
    horizontal_alignment = "right" if max_pin["x_cm"] > 0 else "left"
    ax.annotate(
        f"макс. №{int(max_pin['position_id'])}: {max_pin['relative_power']:.3f}",
        (max_pin["x_cm"], max_pin["y_cm"]), xytext=(offset_x, 18), textcoords="offset points",
        ha=horizontal_alignment, fontsize=9, weight="bold",
        arrowprops={"arrowstyle": "->", "color": "black"},
    )
    cbar = fig.colorbar(points, ax=ax, shrink=0.78, pad=0.02)
    cbar.set_label(r"Відносна потужність твела $P_i/\overline{P}_{312}$")
    ax.set_title("ТВЗ ВВЕР-1000 UGD: відносна потужність твелів")
    ax.set_xlabel("x, см")
    ax.set_ylabel("y, см")
    ax.set_aspect("equal")
    ax.grid(alpha=0.12)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def process_results(statepoint_path: Path, metadata: dict, run_dir: Path, p: dict) -> dict:
    with openmc.StatePoint(statepoint_path) as sp:
        keff = sp.keff
        pin_tally = sp.get_tally(name="pin_results")
        cell_to_record = {r["fuel_cell_id"]: r for r in metadata["fuel_records"]}
        pin_rows = []
        for cell_index, cell_bin in enumerate(pin_tally.filters[0].bins):
            cell_id = int(cell_bin[0] if isinstance(cell_bin, tuple) else cell_bin)
            record = cell_to_record[cell_id]
            row = dict(record)
            for score_index, score in enumerate(pin_tally.scores):
                key = score.replace("-", "_")
                row[key] = float(pin_tally.mean[cell_index, 0, score_index])
                row[f"{key}_std_dev"] = float(pin_tally.std_dev[cell_index, 0, score_index])
            pin_rows.append(row)
        pin_df = pd.DataFrame(pin_rows).sort_values("position_id")

        invalid_pins = pin_df[
            (pin_df["flux"] <= 0.0) | (pin_df["kappa_fission"] <= 0.0)
        ]
        if not invalid_pins.empty:
            examples = ", ".join(
                str(int(position_id)) for position_id in invalid_pins["position_id"].head(12)
            )
            raise RuntimeError(
                f"Виявлено {len(invalid_pins)} паливних комірок із нульовим або "
                f"від'ємним потоком/енерговиділенням (позиції {examples}). "
                "Перевірте геометрію, відповідність tally і статистику розрахунку; "
                "такі результати не можна нормувати."
            )

        power = pin_df["kappa_fission"].to_numpy()
        power_err = pin_df["kappa_fission_std_dev"].to_numpy()
        mean_power, mean_power_err = mean_with_independent_error(power, power_err)
        pin_df["power_fraction"] = power / power.sum()
        pin_df["relative_power"] = power / mean_power
        pin_df["relative_power_std_dev_approx"] = pin_df["relative_power"] * np.sqrt(
            np.square(power_err / power) + (mean_power_err / mean_power) ** 2
        )
        ukrainian_table(pin_df).to_csv(
            run_dir / "potuzhnist_tveliv_tvz_331.csv", index=False
        )

        assembly_rows = scalar_rows(sp.get_tally(name="assembly_rates"))
        ukrainian_table(pd.DataFrame(assembly_rows)).to_csv(
            run_dir / "intehralni_reaktsii_tvz_331.csv", index=False
        )
        assembly = {row["score"]: row for row in assembly_rows}

        spectrum_tally = sp.get_tally(name="energy_group_flux")
        edges = p["energy_groups_eV"]
        spectrum_rows = []
        for index in range(len(edges) - 1):
            spectrum_rows.append({
                "energy_low_eV": edges[index],
                "energy_high_eV": edges[index + 1],
                "mean": float(spectrum_tally.mean[index, 0, 0]),
                "std_dev": float(spectrum_tally.std_dev[index, 0, 0]),
            })
        spectrum_df = pd.DataFrame(spectrum_rows)
        total_flux = float(spectrum_df["mean"].sum())
        if not math.isfinite(total_flux) or total_flux <= 0.0:
            raise RuntimeError("Сумарний потік у спектральному tally має бути додатним")
        total_flux_error = float(math.sqrt(np.square(spectrum_df["std_dev"]).sum()))
        spectrum_df["flux_fraction"] = spectrum_df["mean"] / total_flux
        # Форма через похідні не ділить на mean окремої групи та коректно
        # обробляє групу з нульовою оцінкою потоку.
        spectrum_df["flux_fraction_std_dev_approx"] = np.sqrt(
            np.square(spectrum_df["std_dev"] / total_flux)
            + np.square(spectrum_df["mean"] * total_flux_error / total_flux**2)
        )
        spectrum_rows = spectrum_df.to_dict(orient="records")
        ukrainian_table(spectrum_df).to_csv(
            run_dir / "spektr_neitronnoho_potoku_tvz_331.csv", index=False
        )

        isotope_tally = sp.get_tally(name="fuel_isotope_rates")
        isotope_rows = []
        for material_index, material_bin in enumerate(isotope_tally.filters[0].bins):
            material_id = int(material_bin[0] if isinstance(material_bin, tuple) else material_bin)
            material_name = metadata["materials"]["uo2"].name if material_id == metadata["materials"]["uo2"].id else metadata["materials"]["ugd"].name
            for nuclide_index, nuclide in enumerate(isotope_tally.nuclides):
                for score_index, score in enumerate(isotope_tally.scores):
                    isotope_rows.append({
                        "material": material_name,
                        "nuclide": nuclide,
                        "score": score,
                        "mean": float(isotope_tally.mean[material_index, nuclide_index, score_index]),
                        "std_dev": float(isotope_tally.std_dev[material_index, nuclide_index, score_index]),
                    })
        isotope_df = pd.DataFrame(isotope_rows)
        ukrainian_table(isotope_df).to_csv(
            run_dir / "reaktsii_izotopiv_palyva_tvz_331.csv", index=False
        )

    ordinary = pin_df[pin_df["type"] == "fuel"]
    gd = pin_df[pin_df["type"] == "gd"]
    ordinary_mean, ordinary_err = mean_with_independent_error(
        ordinary["kappa_fission"].to_numpy(), ordinary["kappa_fission_std_dev"].to_numpy()
    )
    gd_mean, gd_err = mean_with_independent_error(
        gd["kappa_fission"].to_numpy(), gd["kappa_fission_std_dev"].to_numpy()
    )
    ordinary_relative, ordinary_relative_err = ratio_with_error(
        ordinary_mean, ordinary_err, mean_power, mean_power_err
    )
    gd_relative, gd_relative_err = ratio_with_error(gd_mean, gd_err, mean_power, mean_power_err)
    max_row = pin_df.loc[pin_df["relative_power"].idxmax()]
    kpin = float(max_row["relative_power"])
    kpin_err = float(max_row["relative_power_std_dev_approx"])

    u238_capture = isotope_df[(isotope_df["nuclide"] == "U238") & (isotope_df["score"] == "(n,gamma)")]
    u235_fission = isotope_df[(isotope_df["nuclide"] == "U235") & (isotope_df["score"] == "fission")]
    cr, cr_err = ratio_with_error(
        u238_capture["mean"].sum(), math.sqrt(np.square(u238_capture["std_dev"]).sum()),
        u235_fission["mean"].sum(), math.sqrt(np.square(u235_fission["std_dev"]).sum()),
    )
    eta_like, eta_like_err = ratio_with_error(
        assembly["nu-fission"]["mean"], assembly["nu-fission"]["std_dev"],
        assembly["absorption"]["mean"], assembly["absorption"]["std_dev"],
    )
    gd_abs = {}
    for nuclide in ("Gd155", "Gd157"):
        rows = isotope_df[(isotope_df["nuclide"] == nuclide) & (isotope_df["score"] == "absorption")]
        gd_abs[nuclide] = {
            "mean_per_source": float(rows["mean"].sum()),
            "std_dev_per_source_approx": float(math.sqrt(np.square(rows["std_dev"]).sum())),
        }

    summary = {
        "keff": float(keff.nominal_value),
        "keff_std_dev": float(keff.std_dev),
        "leakage_expected": 0.0 if p["geometry"]["radial_boundary_type"] == "reflective" and p["geometry"]["axial_boundary_type"] == "reflective" else None,
        "number_of_positions": 331,
        "number_of_fuel_pins": 312,
        "number_of_uo2_pins": 300,
        "number_of_gd_pins": 12,
        "maximum_power_pin": {
            "position_id": int(max_row["position_id"]),
            "type": max_row["type"],
            "axial_x": int(max_row["axial_x"]),
            "axial_alpha": int(max_row["axial_alpha"]),
            "x_cm": float(max_row["x_cm"]),
            "y_cm": float(max_row["y_cm"]),
            "relative_power": kpin,
            "relative_power_std_dev_approx": kpin_err,
        },
        "K_pin": kpin,
        "K_pin_std_dev_approx": kpin_err,
        "maximum_pin_kappa_fission_per_source_eV": float(max_row["kappa_fission"]),
        "maximum_pin_kappa_fission_std_dev_eV": float(max_row["kappa_fission_std_dev"]),
        "all_fuel_mean_kappa_fission_per_source_eV": mean_power,
        "all_fuel_mean_std_dev_approx": mean_power_err,
        "ordinary_pin_mean_kappa_fission_per_source_eV": ordinary_mean,
        "ordinary_pin_mean_kappa_fission_std_dev_approx_eV": ordinary_err,
        "ordinary_pin_mean_relative_power": ordinary_relative,
        "ordinary_pin_mean_relative_power_std_dev_approx": ordinary_relative_err,
        "gd_pin_mean_kappa_fission_per_source_eV": gd_mean,
        "gd_pin_mean_kappa_fission_std_dev_approx_eV": gd_err,
        "gd_pin_mean_relative_power": gd_relative,
        "gd_pin_mean_relative_power_std_dev_approx": gd_relative_err,
        "pin_power_coefficient_of_variation": float(power.std(ddof=1) / mean_power),
        "conversion_ratio_U238_capture_over_U235_fission": cr,
        "conversion_ratio_std_dev_approx": cr_err,
        "nu_fission_over_absorption": eta_like,
        "nu_fission_over_absorption_std_dev_approx": eta_like_err,
        "gd_absorption": gd_abs,
        "energy_group_flux": spectrum_rows,
        "uncertainty_note": "Ratio/mean uncertainties neglect covariance between tally bins.",
    }
    (run_dir / "pidsumok_rozrakhunku.json").write_text(
        json.dumps(ukrainian_summary(summary), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_power_map(
        pin_df, metadata["records"], run_dir / "karta_vidnosnoi_potuzhnosti_tvz_331.png"
    )
    return summary


def safe_name(value: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in "-_" else "_" for char in value)
    return cleaned.strip("_-") or "assembly_331"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "assembly_331_parameters.json",
    )
    parser.add_argument("--run-name", help="перевизначити run_name з JSON")
    parser.add_argument(
        "--output-root", type=Path,
        default=PROJECT_ROOT / "results" / "02_assembly_331_2d" / "runs",
    )
    parser.add_argument(
        "--build-only", action="store_true",
        help="експортувати XML і рисунки без запуску перенесення нейтронів"
    )
    parser.add_argument(
        "--process-existing", type=Path,
        help="повторно обробити наявну папку без нового перенесення нейтронів",
    )
    parser.add_argument(
        "--replot-existing", type=Path,
        help="перегенерувати геометричний рисунок у наявній папці",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.replot_existing:
        run_dir = args.replot_existing
        parameter_file = run_dir / "vykorystani_parametry.json"
        if not parameter_file.exists():
            parameter_file = run_dir / "parameters_used.json"
        p = load_parameters(parameter_file)
        model, metadata = build_model(p)
        model.plot_geometry(cwd=run_dir)
        position_table = pd.DataFrame(metadata["records"]).drop(
            columns=["fuel_cell_id"], errors="ignore"
        )
        ukrainian_table(position_table).to_csv(
            run_dir / "karta_pozytsii_tvz_331.csv", index=False
        )
        (run_dir / "vykorystani_parametry.json").write_text(
            json.dumps(p, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print(f"Геометричний рисунок перегенеровано: {run_dir.resolve()}")
        return
    if args.process_existing:
        run_dir = args.process_existing
        parameter_file = run_dir / "vykorystani_parametry.json"
        if not parameter_file.exists():
            parameter_file = run_dir / "parameters_used.json"
        p = load_parameters(parameter_file)
        model, metadata = build_model(p)
        model.plot_geometry(cwd=run_dir)
        (run_dir / "vykorystani_parametry.json").write_text(
            json.dumps(p, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        statepoint = run_dir / f"statepoint.{p['simulation']['batches']}.h5"
        summary = process_results(statepoint, metadata, run_dir, p)
        print(f"Результати повторно оброблено: {run_dir.resolve()}")
        print(f"k_inf = {summary['keff']:.6f} +/- {summary['keff_std_dev']:.6f}")
        return
    p = load_parameters(args.config)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    label = safe_name(args.run_name or p.get("run_name", "assembly_331"))
    run_dir = args.output_root / f"{timestamp}_{label}"
    suffix = 2
    while run_dir.exists():
        run_dir = args.output_root / f"{timestamp}_{label}_{suffix}"
        suffix += 1
    run_dir.mkdir(parents=True)
    p["run_name"] = label
    (run_dir / "vykorystani_parametry.json").write_text(
        json.dumps(p, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    model, metadata = build_model(p)
    position_table = pd.DataFrame(metadata["records"]).drop(
        columns=["fuel_cell_id"], errors="ignore"
    )
    ukrainian_table(position_table).to_csv(
        run_dir / "karta_pozytsii_tvz_331.csv", index=False
    )
    model.export_to_model_xml(path=run_dir / "model.xml")
    model.plot_geometry(cwd=run_dir)
    if args.build_only:
        print(f"Геометрію побудовано: {run_dir.resolve()}")
        return
    statepoint = Path(model.run(cwd=run_dir, threads=p["simulation"]["threads"]))
    summary = process_results(statepoint, metadata, run_dir, p)
    print("\nРезультати")
    print(f"  k_inf = {summary['keff']:.6f} +/- {summary['keff_std_dev']:.6f}")
    print(
        f"  max pin #{summary['maximum_power_pin']['position_id']}: "
        f"K_pin = {summary['K_pin']:.4f} +/- {summary['K_pin_std_dev_approx']:.4f}"
    )
    print(
        f"  ordinary mean = {summary['ordinary_pin_mean_relative_power']:.4f}; "
        f"Gd mean = {summary['gd_pin_mean_relative_power']:.4f}"
    )
    print(f"  папка результатів: {run_dir.resolve()}")


if __name__ == "__main__":
    main()
