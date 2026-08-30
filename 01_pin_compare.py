"""Порівняння нескінченних комірок ВВЕР із паливом UO2 та UO2-Gd2O3.

Модель використовує відбивний шар заввишки 1 см, еквівалентний нескінченній
періодичній моделі в аксіальному напрямку. OpenMC нормує всі інтегральні
результати tally на одну джерельну частинку.
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


SCORE_UK = {
    "flux": "потік",
    "absorption": "поглинання",
    "fission": "поділ",
    "(n,gamma)": "радіаційне захоплення",
}
REGION_UK = {
    "hole": "центральний гелієвий отвір",
    "fuel": "паливна таблетка",
    "gap": "гелієвий зазор",
    "clad": "оболонка E110",
    "water": "вода",
}
COLUMN_UK = {
    "cell": "номер_комірки",
    "region": "область",
    "score": "показник",
    "mean": "середнє",
    "std_dev": "стандартне_відхилення",
    "volume_cm3": "обʼєм_см3",
    "mean_per_cm3_source": "середнє_на_см3_та_джерельну_частинку",
    "std_dev_per_cm3_source": "похибка_на_см3_та_джерельну_частинку",
    "nuclide": "нуклід",
}


def ukrainian_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Повернути копію таблиці з українськими заголовками та категоріями."""
    result = dataframe.copy()
    if "score" in result:
        result["score"] = result["score"].replace(SCORE_UK)
    if "region" in result:
        result["region"] = result["region"].replace(REGION_UK)
    return result.rename(columns=COLUMN_UK)


def ukrainian_case_summary(summary: dict) -> dict:
    """Підготувати україномовний підсумок одного варіанта комірки."""
    reaction_names = {
        "flux": "потік",
        "absorption": "поглинання",
        "fission": "поділи",
        "(n,gamma)": "радіаційне_захоплення_n_gamma",
    }

    def rate_block(data: dict) -> dict:
        return {
            reaction_names.get(reaction, reaction): {
                "середнє_на_джерельну_частинку": values["mean_per_source"],
                "похибка_на_джерельну_частинку": values["std_dev_per_source"],
            }
            for reaction, values in data.items()
        }

    parameters = summary["parameters"]
    return {
        "варіант": "UO₂" if summary["case"] == "uo2" else "UO₂–Gd₂O₃",
        "k_inf": summary["keff"],
        "похибка_k_inf": summary["keff_std_dev"],
        "інтегральні_реакції": rate_block(summary["global_rates"]),
        "поглинання_гадолінієм": rate_block(summary["gd_absorption"]),
        "параметри": {
            "кількість_частинок": parameters["particles"],
            "кількість_пакетів": parameters["batches"],
            "неактивні_пакети": parameters["inactive"],
            "активні_покоління": parameters["active_generations"],
            "крок_решітки_см": parameters["pitch_cm"],
            "радіус_отвору_см": parameters["hole_radius_cm"],
            "радіус_палива_см": parameters["fuel_radius_cm"],
            "внутрішній_радіус_оболонки_см": parameters["clad_inner_radius_cm"],
            "зовнішній_радіус_оболонки_см": parameters["clad_outer_radius_cm"],
            "температура_палива_K": parameters["fuel_temperature_K"],
            "температура_сповільнювача_K": parameters["moderator_temperature_K"],
        },
    }


def load_parameters(path: Path) -> dict:
    """Прочитати, нормалізувати та фізично перевірити параметри JSON."""
    with path.open(encoding="utf-8") as stream:
        p = json.load(stream)

    if not isinstance(p, dict):
        raise ValueError("Кореневий елемент JSON має бути об'єктом")

    # Сумісність зі старими файлами 01: надалі використовуємо ті самі назви,
    # що й у моделях 02/03.
    if "coolant" in p:
        if "moderator" in p and p["moderator"] != p["coolant"]:
            raise ValueError("Одночасно задано різні секції moderator і coolant")
        p["moderator"] = p.pop("coolant")
    geometry = p.get("geometry")
    if isinstance(geometry, dict) and "hole_radius_cm" in geometry:
        if ("fuel_hole_radius_cm" in geometry
                and geometry["fuel_hole_radius_cm"] != geometry["hole_radius_cm"]):
            raise ValueError(
                "Одночасно задано різні geometry.fuel_hole_radius_cm і geometry.hole_radius_cm"
            )
        geometry["fuel_hole_radius_cm"] = geometry.pop("hole_radius_cm")

    required = {
        "fuel": [
            "uo2_u235_enrichment_wt_percent", "ugd_u235_enrichment_wt_percent",
            "gd2o3_wt_percent", "temperature_K", "uo2_density_g_cm3",
            "gd2o3_density_g_cm3",
        ],
        "moderator": ["temperature_K", "density_g_cm3"],
        "cladding": ["temperature_K"],
        "helium": ["temperature_K", "density_g_cm3"],
        "geometry": [
            "pin_pitch_cm", "fuel_hole_radius_cm", "fuel_radius_cm",
            "clad_inner_radius_cm", "clad_outer_radius_cm", "z_min_cm", "z_max_cm",
        ],
        "simulation": ["particles", "batches", "inactive", "threads"],
    }
    for section, keys in required.items():
        if section not in p:
            raise ValueError(f"Відсутня секція параметрів: {section}")
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

    fuel = p["fuel"]
    for key in ("uo2_u235_enrichment_wt_percent", "ugd_u235_enrichment_wt_percent"):
        enrichment = finite("fuel", key)
        if not 0.0 < enrichment <= 20.0:
            raise ValueError(f"fuel.{key} має бути в інтервалі (0, 20] мас.%")
    gd_fraction = finite("fuel", "gd2o3_wt_percent")
    if not 0.0 <= gd_fraction < 100.0:
        raise ValueError("fuel.gd2o3_wt_percent має бути в інтервалі [0, 100) мас.%")
    positive_density("fuel", "uo2_density_g_cm3", 11.5)
    positive_density("fuel", "gd2o3_density_g_cm3", 10.0)
    positive_density("moderator", "density_g_cm3", 1.5)
    positive_density("helium", "density_g_cm3", 0.1)

    for section in ("fuel", "moderator", "cladding", "helium"):
        temperature = finite(section, "temperature_K")
        if not 294.0 <= temperature <= 1200.0:
            raise ValueError(f"{section}.temperature_K має бути в інтервалі [294, 1200] K")

    g = p["geometry"]
    pitch = finite("geometry", "pin_pitch_cm")
    hole = finite("geometry", "fuel_hole_radius_cm")
    fuel_radius = finite("geometry", "fuel_radius_cm")
    clad_inner = finite("geometry", "clad_inner_radius_cm")
    clad_outer = finite("geometry", "clad_outer_radius_cm")
    z_min = finite("geometry", "z_min_cm")
    z_max = finite("geometry", "z_max_cm")
    if not (0.0 <= hole < fuel_radius < clad_inner < clad_outer < 0.5 * pitch):
        raise ValueError(
            "Радіуси мають задовольняти: 0 <= отвір < паливо < внутрішня оболонка "
            "< зовнішня оболонка < половина кроку решітки"
        )
    if z_min >= z_max:
        raise ValueError("geometry.z_min_cm має бути меншим за geometry.z_max_cm")

    sim = p["simulation"]
    particles = positive_integer("simulation", "particles")
    batches = positive_integer("simulation", "batches")
    threads = positive_integer("simulation", "threads")
    inactive = sim["inactive"]
    if isinstance(inactive, bool) or not isinstance(inactive, int) or not 0 <= inactive < batches:
        raise ValueError("simulation.inactive має бути цілим числом у межах [0, batches)")
    if batches - inactive < 30:
        raise ValueError("Потрібно щонайменше 30 активних пакетів для оцінки статистики")
    if threads > particles:
        raise ValueError("simulation.threads не може перевищувати simulation.particles")
    if "run_name" in p and (not isinstance(p["run_name"], str) or not p["run_name"].strip()):
        raise ValueError("run_name має бути непорожнім рядком")
    return p


def summary_parameters(p: dict) -> dict:
    """Повернути плоский набір параметрів, очікуваний форматом підсумку 01."""
    g = p["geometry"]
    sim = p["simulation"]
    return {
        "particles": sim["particles"],
        "batches": sim["batches"],
        "inactive": sim["inactive"],
        "active_generations": sim["batches"] - sim["inactive"],
        "pitch_cm": g["pin_pitch_cm"],
        "hole_radius_cm": g["fuel_hole_radius_cm"],
        "fuel_radius_cm": g["fuel_radius_cm"],
        "clad_inner_radius_cm": g["clad_inner_radius_cm"],
        "clad_outer_radius_cm": g["clad_outer_radius_cm"],
        "fuel_temperature_K": p["fuel"]["temperature_K"],
        "moderator_temperature_K": p["moderator"]["temperature_K"],
    }


def make_uo2(enrichment: float, name: str, p: dict) -> openmc.Material:
    """Створити діоксид урану; збагачення задано в масових відсотках U-235."""
    fuel = openmc.Material(name=name)
    fuel.add_element("U", 1.0, enrichment=enrichment)
    fuel.add_element("O", 2.0)
    fuel.set_density("g/cm3", p["fuel"]["uo2_density_g_cm3"])
    fuel.temperature = p["fuel"]["temperature_K"]
    return fuel


def make_materials(case: str, p: dict) -> dict[str, openmc.Material]:
    """Створити матеріали для одного варіанта елементарної комірки."""
    helium = openmc.Material(name="Helium")
    helium.add_nuclide("He4", 1.0)
    helium.set_density("g/cm3", p["helium"]["density_g_cm3"])
    helium.temperature = p["helium"]["temperature_K"]

    # Наближення цирконієвого сплаву E110 за атомними густинами бенчмарка.
    e110 = openmc.Material(name="E110 cladding")
    e110.add_element("Zr", 4.259e-2)
    e110.add_nuclide("Nb93", 4.225e-4)
    e110.add_element("Hf", 6.597e-6)
    e110.set_density("sum")
    e110.temperature = p["cladding"]["temperature_K"]

    water = openmc.Material(name="Water")
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", p["moderator"]["density_g_cm3"])
    water.add_s_alpha_beta("c_H_in_H2O")
    water.temperature = p["moderator"]["temperature_K"]

    if case == "uo2":
        enrichment = p["fuel"]["uo2_u235_enrichment_wt_percent"]
        fuel = make_uo2(enrichment, f"UO2, {enrichment:g} wt% U-235", p)
    elif case == "ugd":
        enrichment = p["fuel"]["ugd_u235_enrichment_wt_percent"]
        gd_percent = p["fuel"]["gd2o3_wt_percent"]
        uo2 = make_uo2(enrichment, f"UO2 component, {enrichment:g} wt% U-235", p)
        gd2o3 = openmc.Material(name="Natural Gd2O3 component")
        gd2o3.add_element("Gd", 2.0)
        gd2o3.add_element("O", 3.0)
        gd2o3.set_density("g/cm3", p["fuel"]["gd2o3_density_g_cm3"])
        gd2o3.temperature = p["fuel"]["temperature_K"]
        fuel = openmc.Material.mix_materials(
            [uo2, gd2o3], [1.0 - gd_percent / 100.0, gd_percent / 100.0],
            percent_type="wo"
        )
        fuel.name = f"UO2-Gd2O3: {enrichment:g} wt% U-235, {gd_percent:g} wt% Gd2O3"
        fuel.temperature = p["fuel"]["temperature_K"]
    else:
        raise ValueError(f"Unknown case: {case}")

    return {"fuel": fuel, "helium": helium, "e110": e110, "water": water}


def analytic_volumes(p: dict) -> dict[str, float]:
    """Обчислити аналітичні об'єми областей шару 1 см, см3."""
    g = p["geometry"]
    height = g["z_max_cm"] - g["z_min_cm"]
    side = g["pin_pitch_cm"] / math.sqrt(3.0)
    hex_area = 3.0 * math.sqrt(3.0) * side**2 / 2.0
    return {
        "hole": math.pi * g["fuel_hole_radius_cm"]**2 * height,
        "fuel": math.pi * (g["fuel_radius_cm"]**2 - g["fuel_hole_radius_cm"]**2) * height,
        "gap": math.pi * (g["clad_inner_radius_cm"]**2 - g["fuel_radius_cm"]**2) * height,
        "clad": math.pi * (g["clad_outer_radius_cm"]**2 - g["clad_inner_radius_cm"]**2) * height,
        "water": (hex_area - math.pi * g["clad_outer_radius_cm"]**2) * height,
    }


def build_model(case: str, p: dict) -> tuple[openmc.Model, dict]:
    """Побудувати відбивну шестикутну модель елементарної комірки."""
    mat = make_materials(case, p)
    g = p["geometry"]
    sim = p["simulation"]

    hole_cyl = openmc.ZCylinder(r=g["fuel_hole_radius_cm"])
    fuel_cyl = openmc.ZCylinder(r=g["fuel_radius_cm"])
    clad_inner = openmc.ZCylinder(r=g["clad_inner_radius_cm"])
    clad_outer = openmc.ZCylinder(r=g["clad_outer_radius_cm"])
    z_low = openmc.ZPlane(z0=g["z_min_cm"], boundary_type="reflective")
    z_high = openmc.ZPlane(z0=g["z_max_cm"], boundary_type="reflective")
    axial = +z_low & -z_high

    side = g["pin_pitch_cm"] / math.sqrt(3.0)
    boundary = openmc.model.HexagonalPrism(
        edge_length=side, orientation="x", boundary_type="reflective"
    )
    inside_hex = -boundary

    cells = {
        "hole": openmc.Cell(name="Central helium hole", fill=mat["helium"], region=-hole_cyl & axial),
        "fuel": openmc.Cell(name="Fuel pellet", fill=mat["fuel"], region=+hole_cyl & -fuel_cyl & axial),
        "gap": openmc.Cell(name="Helium gap", fill=mat["helium"], region=+fuel_cyl & -clad_inner & axial),
        "clad": openmc.Cell(name="E110 cladding", fill=mat["e110"], region=+clad_inner & -clad_outer & axial),
        "water": openmc.Cell(name="Water moderator", fill=mat["water"], region=+clad_outer & inside_hex & axial),
    }
    geometry = openmc.Geometry(openmc.Universe(cells=list(cells.values())))

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = sim["batches"]
    settings.inactive = sim["inactive"]
    settings.particles = sim["particles"]
    settings.seed = 1701 if case == "uo2" else 1702
    settings.temperature = {
        "method": "interpolation",
        "range": (294.0, 1200.0),
    }
    source_box = openmc.stats.Box(
        [-g["fuel_radius_cm"], -g["fuel_radius_cm"], g["z_min_cm"]],
        [g["fuel_radius_cm"], g["fuel_radius_cm"], g["z_max_cm"]],
    )
    settings.source = openmc.IndependentSource(
        space=source_box, constraints={"fissionable": True}
    )

    global_tally = openmc.Tally(name="global_rates")
    global_tally.scores = ["flux", "absorption", "fission", "(n,gamma)"]

    region_tally = openmc.Tally(name="rates_by_region")
    region_tally.filters = [openmc.CellFilter(list(cells.values()))]
    region_tally.scores = ["flux", "absorption", "fission", "(n,gamma)"]

    tallies = openmc.Tallies([global_tally, region_tally])
    if case == "ugd":
        gd_tally = openmc.Tally(name="gd_absorption")
        gd_tally.filters = [openmc.CellFilter([cells["fuel"]])]
        gd_tally.nuclides = ["Gd155", "Gd157"]
        gd_tally.scores = ["absorption"]
        tallies.append(gd_tally)

    plot = openmc.Plot(name=f"{case}_geometry")
    plot.filename = (
        "heometriia_komirky_uo2" if case == "uo2"
        else "heometriia_komirky_uo2_gd2o3"
    )
    plot.basis = "xy"
    plot.origin = (0.0, 0.0, 0.0)
    plot.width = (1.45, 1.45)
    plot.pixels = (900, 900)
    plot.color_by = "material"
    plot.colors = {
        mat["fuel"]: "gold" if case == "uo2" else "forestgreen",
        mat["helium"]: "white",
        mat["e110"]: "silver",
        mat["water"]: "deepskyblue",
    }

    model = openmc.Model(
        geometry=geometry,
        materials=openmc.Materials(list(mat.values())),
        settings=settings,
        tallies=tallies,
        plots=openmc.Plots([plot]),
    )
    metadata = {"cells": cells, "materials": mat, "volumes": analytic_volumes(p)}
    return model, metadata


def run_case(case: str, output_dir: Path, p: dict) -> dict:
    """Розрахувати один варіант, зберегти рисунки й таблиці та повернути підсумок."""
    case_dir = output_dir / case
    case_dir.mkdir(parents=True, exist_ok=True)
    sim = p["simulation"]
    model, metadata = build_model(case, p)
    model.export_to_model_xml(path=case_dir / "model.xml")
    model.plot_geometry(cwd=case_dir)
    statepoint_path = model.run(cwd=case_dir, threads=sim["threads"])

    with openmc.StatePoint(statepoint_path) as sp:
        keff = sp.keff
        global_tally = sp.get_tally(name="global_rates")
        global_rows = []
        for score_index, score in enumerate(global_tally.scores):
            global_rows.append({
                "score": score,
                "mean": float(global_tally.mean[0, 0, score_index]),
                "std_dev": float(global_tally.std_dev[0, 0, score_index]),
            })
        global_df = pd.DataFrame(global_rows)
        ukrainian_table(global_df).to_csv(
            case_dir / "intehralni_reaktsii_komirky.csv", index=False
        )

        region_tally = sp.get_tally(name="rates_by_region")
        cell_names = {cell.id: key for key, cell in metadata["cells"].items()}
        cell_bins = region_tally.filters[0].bins
        region_rows = []
        for cell_index, cell_bin in enumerate(cell_bins):
            cell_id = int(cell_bin[0] if isinstance(cell_bin, tuple) else cell_bin)
            region = cell_names[cell_id]
            volume = metadata["volumes"][region]
            for score_index, score in enumerate(region_tally.scores):
                mean = float(region_tally.mean[cell_index, 0, score_index])
                std_dev = float(region_tally.std_dev[cell_index, 0, score_index])
                region_rows.append({
                    "cell": cell_id,
                    "region": region,
                    "score": score,
                    "mean": mean,
                    "std_dev": std_dev,
                    "volume_cm3": volume,
                    "mean_per_cm3_source": mean / volume,
                    "std_dev_per_cm3_source": std_dev / volume,
                })
        region_df = pd.DataFrame(region_rows)
        ukrainian_table(region_df).to_csv(
            case_dir / "reaktsii_za_oblastiamy_komirky.csv", index=False
        )

        gd_rates: dict[str, dict[str, float]] = {}
        if case == "ugd":
            gd_tally = sp.get_tally(name="gd_absorption")
            gd_rows = []
            for nuclide_index, nuclide in enumerate(gd_tally.nuclides):
                mean = float(gd_tally.mean[0, nuclide_index, 0])
                std_dev = float(gd_tally.std_dev[0, nuclide_index, 0])
                gd_rows.append({
                    "nuclide": nuclide,
                    "score": "absorption",
                    "mean": mean,
                    "std_dev": std_dev,
                })
                gd_rates[nuclide] = {
                    "mean_per_source": mean,
                    "std_dev_per_source": std_dev,
                }
            gd_df = pd.DataFrame(gd_rows)
            ukrainian_table(gd_df).to_csv(
                case_dir / "pohlynannia_izotopamy_hadoliniiu.csv", index=False
            )

    global_rates = {
        str(row["score"]): {
            "mean_per_source": float(row["mean"]),
            "std_dev_per_source": float(row["std_dev"]),
        }
        for _, row in global_df.iterrows()
    }
    summary = {
        "case": case,
        "keff": float(keff.nominal_value),
        "keff_std_dev": float(keff.std_dev),
        "global_rates": global_rates,
        "gd_absorption": gd_rates,
        "parameters": summary_parameters(p),
    }
    (case_dir / "pidsumok_rozrakhunku.json").write_text(
        json.dumps(ukrainian_case_summary(summary), indent=2, ensure_ascii=False),
        encoding="utf-8"
    )
    return summary


def save_comparison(uo2: dict, ugd: dict, output_dir: Path, p: dict) -> dict:
    """Обчислити ефективність поглинача та створити матеріали порівняння."""
    ku, su = uo2["keff"], uo2["keff_std_dev"]
    kg, sg = ugd["keff"], ugd["keff_std_dev"]
    delta_k = ku - kg
    absorber_worth = 1.0 / kg - 1.0 / ku
    worth_pcm = absorber_worth * 1.0e5
    worth_std_pcm = math.sqrt((su / ku**2) ** 2 + (sg / kg**2) ** 2) * 1.0e5

    comparison = {
        "uo2_keff": ku,
        "uo2_keff_std_dev": su,
        "ugd_keff": kg,
        "ugd_keff_std_dev": sg,
        "delta_k_uo2_minus_ugd": delta_k,
        "gd_absorber_worth_delta_rho": absorber_worth,
        "gd_absorber_worth_pcm": worth_pcm,
        "gd_absorber_worth_std_dev_pcm": worth_std_pcm,
    }
    comparison_uk = pd.DataFrame([comparison]).rename(columns={
        "uo2_keff": "k_inf_uo2",
        "uo2_keff_std_dev": "похибка_k_inf_uo2",
        "ugd_keff": "k_inf_uo2_gd2o3",
        "ugd_keff_std_dev": "похибка_k_inf_uo2_gd2o3",
        "delta_k_uo2_minus_ugd": "різниця_k_uo2_мінус_uo2_gd2o3",
        "gd_absorber_worth_delta_rho": "ефективність_гадолінію_delta_rho",
        "gd_absorber_worth_pcm": "ефективність_гадолінію_pcm",
        "gd_absorber_worth_std_dev_pcm": "похибка_ефективності_гадолінію_pcm",
    })
    comparison_uk.to_csv(output_dir / "porivniannia_k_inf.csv", index=False)
    comparison_json_uk = comparison_uk.iloc[0].to_dict()
    (output_dir / "porivniannia_k_inf.json").write_text(
        json.dumps(comparison_json_uk, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    fuel = p["fuel"]
    labels = [
        f"UO₂, {fuel['uo2_u235_enrichment_wt_percent']:g}% U-235",
        "UO₂-Gd₂O₃\n"
        f"{fuel['ugd_u235_enrichment_wt_percent']:g}% U-235; "
        f"{fuel['gd2o3_wt_percent']:g}% Gd₂O₃",
    ]
    values = np.array([ku, kg])
    errors = np.array([su, sg])
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    bars = ax.bar(labels, values, yerr=errors, capsize=5, color=["#e9b949", "#2e8b57"])
    ax.set_ylabel(r"$k_\infty$")
    ax.set_title("Порівняння нескінченних елементарних комірок")
    ax.grid(axis="y", alpha=0.25)
    lower = max(0.0, min(values) - 0.12)
    ax.set_ylim(lower, max(values) + 0.08)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.012, f"{value:.5f}", ha="center")
    fig.tight_layout()
    fig.savefig(output_dir / "porivniannia_k_inf_komirok_uo2_ta_uo2_gd2o3.png", dpi=180)
    plt.close(fig)
    return comparison


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("pin_parameters.json"),
        help="JSON-файл з усіма параметрами моделі та розрахунку",
    )
    parser.add_argument(
        "--run-name", help="перевизначити run_name з файла параметрів"
    )
    parser.add_argument(
        "--output-root", type=Path, default=Path("results/01_pin_compare/runs"),
        help="коренева папка; всередині створюється окрема папка запуску",
    )
    parser.add_argument(
        "--replot-existing-root", type=Path,
        help="перегенерувати графіки у наявній папці без нового розрахунку",
    )
    return parser.parse_args()


def safe_run_name(value: str) -> str:
    """Повернути безпечну для файлової системи назву запуску."""
    cleaned = "".join(char if char.isalnum() or char in "-_" else "_" for char in value)
    return cleaned.strip("_-") or "pin_compare"


def main() -> None:
    args = parse_args()
    p = load_parameters(args.config)
    if args.replot_existing_root:
        root = args.replot_existing_root
        summaries = {}
        for case in ("uo2", "ugd"):
            model, _ = build_model(case, p)
            model.plot_geometry(cwd=root / case)
            summary_file = root / case / "pidsumok_rozrakhunku.json"
            if not summary_file.exists():
                summary_file = root / case / "summary.json"
            summaries[case] = json.loads(summary_file.read_text(encoding="utf-8"))
            if "keff" not in summaries[case]:
                summaries[case] = {
                    "case": case,
                    "keff": summaries[case]["k_inf"],
                    "keff_std_dev": summaries[case]["похибка_k_inf"],
                }
        save_comparison(summaries["uo2"], summaries["ugd"], root, p)
        print(f"Графіки перегенеровано: {root.resolve()}")
        return
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    label = safe_run_name(args.run_name or p.get("run_name", "pin_compare"))
    run_dir = args.output_root / f"{timestamp}_{label}"
    # Збіг у межах однієї секунди малоймовірний, але результати не можна змішувати.
    suffix = 2
    while run_dir.exists():
        run_dir = args.output_root / f"{timestamp}_{label}_{suffix}"
        suffix += 1
    run_dir.mkdir(parents=True)
    p["run_name"] = label
    (run_dir / "vykorystani_parametry.json").write_text(
        json.dumps(p, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("Розрахунок комірки UO2...")
    uo2 = run_case("uo2", run_dir, p)
    print("Розрахунок комірки UO2-Gd2O3...")
    ugd = run_case("ugd", run_dir, p)
    comparison = save_comparison(uo2, ugd, run_dir, p)
    print("\nРезультати")
    print(f"  UO2 k_inf       = {uo2['keff']:.6f} +/- {uo2['keff_std_dev']:.6f}")
    print(f"  UO2-Gd2O3 k_inf = {ugd['keff']:.6f} +/- {ugd['keff_std_dev']:.6f}")
    print(
        "  Ефективність Gd = "
        f"{comparison['gd_absorber_worth_pcm']:.1f} +/- "
        f"{comparison['gd_absorber_worth_std_dev_pcm']:.1f} pcm"
    )
    print(f"  Папка результатів: {run_dir.resolve()}")


if __name__ == "__main__":
    main()
