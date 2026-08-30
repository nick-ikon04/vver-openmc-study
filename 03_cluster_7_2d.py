"""Параметризована 2D-модель кластера з центральної та шести сусідніх ТВЗ."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import openmc
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent


def load_assembly_module():
    """Завантажити перевірені будівельні функції моделі однієї ТВЗ."""
    path = Path(__file__).with_name("02_assembly_331_2d.py")
    spec = importlib.util.spec_from_file_location("assembly_331", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Не вдалося завантажити модуль геометрії ТВЗ")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ASSEMBLY = load_assembly_module()


def load_parameters(path: Path) -> dict:
    """Прочитати параметри та перевірити фізику ТВЗ і геометрію кластера."""
    with path.open(encoding="utf-8") as stream:
        p = json.load(stream)

    if not isinstance(p, dict):
        raise ValueError("Кореневий елемент JSON має бути об'єктом")
    required = {
        "fuel": [],
        "moderator": [],
        "cladding": [],
        "helium": [],
        "geometry": [
            "pin_pitch_cm", "assembly_pitch_cm", "cluster_flat_to_flat_cm",
            "fuel_hole_radius_cm", "fuel_radius_cm", "clad_inner_radius_cm",
            "clad_outer_radius_cm", "guide_inner_radius_cm", "guide_outer_radius_cm",
            "central_inner_radius_cm", "central_outer_radius_cm", "z_min_cm", "z_max_cm",
            "radial_boundary_type", "axial_boundary_type",
        ],
        "assembly_layout": ["hex_radius", "central_coordinate", "guide_orbit_seeds"],
        "assembly_types": [],
        "cluster": ["hex_radius", "assemblies"],
        "simulation": [],
    }
    missing_sections = [section for section in required if section not in p]
    if missing_sections:
        raise ValueError(f"Відсутні розділи параметрів: {', '.join(missing_sections)}")
    if "energy_groups_eV" not in p:
        raise ValueError("Відсутній параметр energy_groups_eV")
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

    cluster_radius = p["cluster"]["hex_radius"]
    if isinstance(cluster_radius, bool) or not isinstance(cluster_radius, int) or cluster_radius != 1:
        raise ValueError("Семикасетний кластер вимагає cluster.hex_radius = 1")

    assembly_types = p["assembly_types"]
    if not assembly_types:
        raise ValueError("assembly_types має містити щонайменше один тип ТВЗ")
    for name, type_parameters in assembly_types.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(type_parameters, dict):
            raise ValueError("Кожен тип ТВЗ має мати непорожню назву та об'єкт параметрів")
        if "gd_orbit_seeds" not in type_parameters:
            raise ValueError(f"Відсутній assembly_types.{name}.gd_orbit_seeds")
        # Повна перевірка матеріалів, радіусів, орбіт, меж, статистики та спектра
        # виконується тим самим валідатором, що й для окремої ТВЗ.
        ASSEMBLY.validate_parameters(assembly_parameters(p, name))

    assemblies = p["cluster"]["assemblies"]
    if not isinstance(assemblies, list):
        raise ValueError("cluster.assemblies має бути списком")
    for item in assemblies:
        if not isinstance(item, dict):
            raise ValueError("Кожен елемент cluster.assemblies має бути JSON-об'єктом")
        missing = [key for key in ("assembly_id", "name", "coordinate", "type") if key not in item]
        if missing:
            raise ValueError(f"У записі ТВЗ відсутні поля: {', '.join(missing)}")
        if (isinstance(item["assembly_id"], bool) or not isinstance(item["assembly_id"], int)
                or item["assembly_id"] <= 0):
            raise ValueError("assembly_id має бути додатним цілим числом")
        if not isinstance(item["name"], str) or not item["name"].strip():
            raise ValueError("Назва ТВЗ має бути непорожнім рядком")
        coordinate = item["coordinate"]
        if not (isinstance(coordinate, (list, tuple)) and len(coordinate) == 2
                and all(isinstance(value, int) and not isinstance(value, bool) for value in coordinate)):
            raise ValueError("Координата кожної ТВЗ має містити два цілі числа")

    coordinates = [tuple(item["coordinate"]) for item in assemblies]
    if len(assemblies) != 7 or len(set(coordinates)) != 7 or (0, 0) not in coordinates:
        raise ValueError("Кластер має містити сім унікальних ТВЗ, включно з координатою [0, 0]")
    if set(coordinates) != set(ASSEMBLY.hex_positions(1)):
        raise ValueError("Координати ТВЗ мають утворювати повний шестикутник радіуса 1")
    if len({item["assembly_id"] for item in assemblies}) != 7:
        raise ValueError("Номери assembly_id мають бути унікальними")
    if len({item["name"] for item in assemblies}) != 7:
        raise ValueError("Назви ТВЗ мають бути унікальними")
    for item in assemblies:
        if item["type"] not in p["assembly_types"]:
            raise ValueError(f"Невідомий тип ТВЗ: {item['type']}")

    assembly_pitch = finite("geometry", "assembly_pitch_cm")
    cluster_flat_to_flat = finite("geometry", "cluster_flat_to_flat_cm")
    if assembly_pitch <= 0.0 or cluster_flat_to_flat <= 0.0:
        raise ValueError("Крок ТВЗ і розмір кластера мають бути додатними")
    expected_cluster_flat_to_flat = (2 * cluster_radius + 1) * assembly_pitch
    if not math.isclose(cluster_flat_to_flat, expected_cluster_flat_to_flat, rel_tol=0.0, abs_tol=1.0e-8):
        raise ValueError(
            "geometry.cluster_flat_to_flat_cm має дорівнювати "
            f"(2*cluster.hex_radius+1)*assembly_pitch_cm = {expected_cluster_flat_to_flat:.6f} см"
        )
    return p


def assembly_parameters(p: dict, assembly_type: str) -> dict:
    """Підготувати параметри внутрішньої решітки заданого типу ТВЗ."""
    g = p["geometry"]
    layout = p["assembly_layout"]
    return {
        "fuel": p["fuel"], "moderator": p["moderator"],
        "cladding": p["cladding"], "helium": p["helium"],
        "geometry": {
            "pin_pitch_cm": g["pin_pitch_cm"],
            "assembly_flat_to_flat_cm": g["assembly_pitch_cm"],
            "fuel_hole_radius_cm": g["fuel_hole_radius_cm"],
            "fuel_radius_cm": g["fuel_radius_cm"],
            "clad_inner_radius_cm": g["clad_inner_radius_cm"],
            "clad_outer_radius_cm": g["clad_outer_radius_cm"],
            "guide_inner_radius_cm": g["guide_inner_radius_cm"],
            "guide_outer_radius_cm": g["guide_outer_radius_cm"],
            "central_inner_radius_cm": g["central_inner_radius_cm"],
            "central_outer_radius_cm": g["central_outer_radius_cm"],
            "z_min_cm": g["z_min_cm"], "z_max_cm": g["z_max_cm"],
            "radial_boundary_type": g["radial_boundary_type"],
            "axial_boundary_type": g["axial_boundary_type"],
        },
        "layout": {
            "hex_radius": layout["hex_radius"],
            "central_coordinate": layout["central_coordinate"],
            "gd_orbit_seeds": p["assembly_types"][assembly_type]["gd_orbit_seeds"],
            "guide_orbit_seeds": layout["guide_orbit_seeds"],
        },
        "simulation": p["simulation"], "energy_groups_eV": p["energy_groups_eV"],
    }


def build_assembly_universe(item: dict, p: dict, materials: dict) -> tuple[openmc.Universe, list[dict], list[openmc.Cell]]:
    """Побудувати одну унікальну ТВЗ для макрорешітки кластера."""
    ap = assembly_parameters(p, item["type"])
    local_records = ASSEMBLY.layout_records(ap)
    radius = ap["layout"]["hex_radius"]
    rings = [[None] * (6 * (radius - r) if r < radius else 1) for r in range(radius + 1)]
    assembly_x, assembly_y = ASSEMBLY.lattice_xy(
        tuple(item["coordinate"]), p["geometry"]["assembly_pitch_cm"], orientation="x"
    )
    records: list[dict] = []
    fuel_cells: list[openmc.Cell] = []
    for local in local_records:
        record = dict(local)
        record.update({
            "assembly_id": item["assembly_id"], "assembly_name": item["name"],
            "assembly_type": item["type"], "local_x_cm": local["x_cm"],
            "local_y_cm": local["y_cm"], "x_cm": local["x_cm"] + assembly_x,
            "y_cm": local["y_cm"] + assembly_y,
        })
        if record["type"] in ("fuel", "gd"):
            universe, fuel_cell = ASSEMBLY.fuel_universe(
                record, materials["uo2"] if record["type"] == "fuel" else materials["ugd"],
                materials, ap,
            )
            record["fuel_cell_id"] = fuel_cell.id
            fuel_cells.append(fuel_cell)
        else:
            universe = ASSEMBLY.tube_universe(record, record["type"], materials, ap)
        ring, within = ASSEMBLY.ring_index(
            (record["axial_x"], record["axial_alpha"]), radius + 1
        )
        rings[ring][within] = universe
        records.append(record)
    if any(universe is None for ring in rings for universe in ring):
        raise RuntimeError(f"Не всі 331 позиції ТВЗ {item['assembly_id']} заповнено")
    lattice = openmc.HexLattice(name=f"Решітка ТВЗ {item['assembly_id']}")
    lattice.center = (0.0, 0.0)
    lattice.pitch = (p["geometry"]["pin_pitch_cm"],)
    lattice.orientation = "y"
    lattice.universes = rings
    lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=materials["water"])])
    return openmc.Universe(
        name=f"ТВЗ {item['assembly_id']} {item['name']}",
        cells=[openmc.Cell(fill=lattice)],
    ), records, fuel_cells


def build_model(p: dict) -> tuple[openmc.Model, dict]:
    """Побудувати дворівневу решітку: твели всередині семи ТВЗ."""
    materials = ASSEMBLY.make_materials(assembly_parameters(p, next(iter(p["assembly_types"]))))
    for key, name in {
        "uo2": "Паливо UO₂", "ugd": "Паливо UO₂–Gd₂O₃", "helium": "Гелій",
        "e110": "Оболонка E110", "water": "Борована вода",
    }.items():
        materials[key].name = name
    macro_rings: list[list[openmc.Universe | None]] = [[None] * 6, [None]]
    records: list[dict] = []
    fuel_cells: list[openmc.Cell] = []
    for item in p["cluster"]["assemblies"]:
        universe, assembly_records, assembly_fuel_cells = build_assembly_universe(item, p, materials)
        ring, within = ASSEMBLY.ring_index(tuple(item["coordinate"]), 2)
        macro_rings[ring][within] = universe
        records.extend(assembly_records)
        fuel_cells.extend(assembly_fuel_cells)
    if any(universe is None for ring in macro_rings for universe in ring):
        raise RuntimeError("Не всі позиції макрорешітки заповнено")

    cluster_lattice = openmc.HexLattice(name="Кластер із семи ТВЗ")
    cluster_lattice.center = (0.0, 0.0)
    cluster_lattice.pitch = (p["geometry"]["assembly_pitch_cm"],)
    # Внутрішня картограма кожної ТВЗ має гострі вершини зверху/знизу.
    # Для HexLattice така форма макрокомірки відповідає orientation='x'.
    # orientation='y' обрізала шість вершин кожної ТВЗ, бо в HexLattice та
    # HexagonalPrism однакові літери позначають різні геометричні домовленості.
    cluster_lattice.orientation = "x"
    cluster_lattice.universes = macro_rings
    cluster_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=materials["water"])])

    g = p["geometry"]
    boundary = openmc.model.HexagonalPrism(
        edge_length=g["cluster_flat_to_flat_cm"] / math.sqrt(3.0), orientation="y",
        boundary_type=g["radial_boundary_type"],
    )
    z_low = openmc.ZPlane(z0=g["z_min_cm"], boundary_type=g["axial_boundary_type"])
    z_high = openmc.ZPlane(z0=g["z_max_cm"], boundary_type=g["axial_boundary_type"])
    root = openmc.Cell(fill=cluster_lattice, region=-boundary & +z_low & -z_high)

    settings = openmc.Settings()
    sim = p["simulation"]
    settings.run_mode = "eigenvalue"
    settings.particles = sim["particles"]
    settings.batches = sim["batches"]
    settings.inactive = sim["inactive"]
    settings.seed = sim["seed"]
    settings.temperature = {"method": "interpolation", "range": (294.0, 1200.0)}
    half = 0.5 * g["cluster_flat_to_flat_cm"]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box([-half, -half, g["z_min_cm"]], [half, half, g["z_max_cm"]]),
        constraints={"fissionable": True},
    )

    pin_tally = openmc.Tally(name="potuzhnist_tveliv")
    pin_tally.filters = [openmc.CellFilter(fuel_cells)]
    pin_tally.scores = ["kappa-fission", "fission", "flux", "absorption"]
    total_tally = openmc.Tally(name="intehralni_reaktsii_klastera")
    total_tally.scores = ["flux", "absorption", "fission", "nu-fission", "kappa-fission"]
    spectrum = openmc.Tally(name="spektr_neitronnoho_potoku")
    spectrum.filters = [openmc.EnergyFilter(p["energy_groups_eV"])]
    spectrum.scores = ["flux"]

    plot = openmc.Plot(name="Геометрія кластера")
    plot.filename = "heometriia_klastera_iz_semy_tvz"
    plot.basis = "xy"
    plot.origin = (0.0, 0.0, 0.0)
    plot.width = (74.0, 82.0)
    plot.pixels = (1200, 1320)
    plot.color_by = "material"
    plot.colors = {
        materials["uo2"]: "gold", materials["ugd"]: "forestgreen",
        materials["helium"]: "white", materials["e110"]: "slategray",
        materials["water"]: "deepskyblue",
    }
    model = openmc.Model(
        geometry=openmc.Geometry(openmc.Universe(cells=[root])),
        materials=openmc.Materials(list(materials.values())), settings=settings,
        tallies=openmc.Tallies([pin_tally, total_tally, spectrum]), plots=openmc.Plots([plot]),
    )
    return model, {"records": records, "fuel_cells": fuel_cells, "materials": materials}


def make_power_map(pin_df: pd.DataFrame, p: dict, output: Path) -> None:
    """Зберегти кольорову карту потужності 2184 паливних твелів."""
    fig, ax = plt.subplots(figsize=(15, 16))
    points = ax.scatter(
        pin_df["x_cm"], pin_df["y_cm"], c=pin_df["relative_power"], cmap="turbo",
        marker="h", s=31, linewidths=0.08, edgecolors="0.3",
    )
    gd = pin_df[pin_df["type"] == "gd"]
    ax.scatter(gd["x_cm"], gd["y_cm"], marker="h", s=48, facecolors="none",
               edgecolors="lime", linewidths=0.7, label="Твеги UO₂–Gd₂O₃")
    side = p["geometry"]["assembly_pitch_cm"] / math.sqrt(3.0)
    for item in p["cluster"]["assemblies"]:
        x, y = ASSEMBLY.lattice_xy(
            tuple(item["coordinate"]), p["geometry"]["assembly_pitch_cm"], orientation="x"
        )
        ax.add_patch(RegularPolygon((x, y), 6, radius=side, orientation=0.0,
                                    fill=False, edgecolor="black", linewidth=1.0))
        ax.text(x, y + 0.38 * p["geometry"]["assembly_pitch_cm"],
                f"ТВЗ-{item['assembly_id']}", ha="center", fontsize=7, weight="bold")
    maximum = pin_df.loc[pin_df["relative_power"].idxmax()]
    ax.scatter([maximum["x_cm"]], [maximum["y_cm"]], marker="h", s=90,
               facecolors="none", edgecolors="black", linewidths=1.7)
    cbar = fig.colorbar(points, ax=ax, shrink=0.82)
    cbar.set_label("Відносна потужність твела")
    ax.set_title("Кластер із семи ТВЗ: карта відносної потужності твелів")
    ax.set_xlabel("x, см")
    ax.set_ylabel("y, см")
    ax.set_aspect("equal")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def process_results(statepoint: Path, metadata: dict, p: dict, run_dir: Path) -> dict:
    """Обробити коефіцієнт розмноження і розподіл потужності кластера."""
    with openmc.StatePoint(statepoint) as sp:
        keff = sp.keff
        tally = sp.get_tally(name="potuzhnist_tveliv")
        by_cell = {record["fuel_cell_id"]: record for record in metadata["records"] if "fuel_cell_id" in record}
        rows = []
        for index, cell_bin in enumerate(tally.filters[0].bins):
            cell_id = int(cell_bin[0] if isinstance(cell_bin, tuple) else cell_bin)
            row = dict(by_cell[cell_id])
            for score_index, score in enumerate(tally.scores):
                key = score.replace("-", "_")
                row[key] = float(tally.mean[index, 0, score_index])
                row[f"{key}_std_dev"] = float(tally.std_dev[index, 0, score_index])
            rows.append(row)
        pin_df = pd.DataFrame(rows).sort_values(["assembly_id", "position_id"])
        invalid_pins = pin_df[
            (pin_df["flux"] <= 0.0) | (pin_df["kappa_fission"] <= 0.0)
        ]
        if not invalid_pins.empty:
            examples = ", ".join(
                f"ТВЗ {int(row.assembly_id)} / позиція {int(row.position_id)}"
                for row in invalid_pins.head(12).itertuples()
            )
            raise RuntimeError(
                f"Виявлено {len(invalid_pins)} паливних комірок із нульовим або "
                f"від'ємним потоком/енерговиділенням "
                f"({examples}). Геометрію кластера обрізано або tally не охоплює всі твели; "
                "результати не можна нормувати чи публікувати."
            )
        mean_power = float(pin_df["kappa_fission"].mean())
        if not math.isfinite(mean_power) or mean_power <= 0.0:
            raise RuntimeError("Середня енергія поділів має бути додатною та скінченною")
        mean_power_error = float(
            np.sqrt(np.square(pin_df["kappa_fission_std_dev"]).sum()) / len(pin_df)
        )
        pin_df["relative_power"] = pin_df["kappa_fission"] / mean_power
        pin_df["relative_power_std_dev_approx"] = pin_df["relative_power"] * np.sqrt(
            np.square(pin_df["kappa_fission_std_dev"] / pin_df["kappa_fission"])
            + (mean_power_error / mean_power) ** 2
        )
        assembly_df = pin_df.groupby(["assembly_id", "assembly_name"], as_index=False).agg(
            power=("kappa_fission", "sum"),
            power_std_dev_approx=("kappa_fission_std_dev", lambda x: float(np.sqrt(np.square(x).sum()))),
        )
        assembly_df["relative_power"] = assembly_df["power"] / assembly_df["power"].mean()
        integral = sp.get_tally(name="intehralni_reaktsii_klastera")
        integral_rows = ASSEMBLY.scalar_rows(integral)
        spectrum = sp.get_tally(name="spektr_neitronnoho_potoku")
        edges = p["energy_groups_eV"]
        spectrum_rows = [
            {"energy_low_eV": edges[i], "energy_high_eV": edges[i + 1],
             "mean": float(spectrum.mean[i, 0, 0]), "std_dev": float(spectrum.std_dev[i, 0, 0])}
            for i in range(len(edges) - 1)
        ]

    ASSEMBLY.ukrainian_table(pin_df).to_csv(run_dir / "potuzhnist_tveliv_klastera_7_tvz.csv", index=False)
    ASSEMBLY.ukrainian_table(assembly_df).to_csv(run_dir / "potuzhnist_okremykh_tvz_klastera.csv", index=False)
    ASSEMBLY.ukrainian_table(pd.DataFrame(integral_rows)).to_csv(run_dir / "intehralni_reaktsii_klastera.csv", index=False)
    ASSEMBLY.ukrainian_table(pd.DataFrame(spectrum_rows)).to_csv(run_dir / "spektr_neitronnoho_potoku_klastera.csv", index=False)
    maximum = pin_df.loc[pin_df["relative_power"].idxmax()]
    summary = {
        "k_inf": float(keff.nominal_value), "похибка_k_inf": float(keff.std_dev),
        "кількість_твз": 7, "кількість_позицій": 2317,
        "кількість_паливних_твелів": int(len(pin_df)),
        "кількість_звичайних_твелів": int((pin_df["type"] == "fuel").sum()),
        "кількість_твегів": int((pin_df["type"] == "gd").sum()),
        "найбільш_навантажений_твел": {
            "номер_твз": int(maximum["assembly_id"]),
            "номер_позиції_у_твз": int(maximum["position_id"]),
            "відносна_потужність": float(maximum["relative_power"]),
            "похибка_відносної_потужності_наближено": float(maximum["relative_power_std_dev_approx"]),
        },
        "K_pin": float(maximum["relative_power"]),
        "похибка_K_pin_наближено": float(maximum["relative_power_std_dev_approx"]),
        "потужність_окремих_твз": ASSEMBLY.ukrainian_table(assembly_df).to_dict(orient="records"),
        "примітка": (
            "Зовнішні радіальні та аксіальні межі відбивні. Це відбивний "
            "семикасетний суперелемент, у якому центральна та крайні ТВЗ не є "
            "нейтронно еквівалентними; результат не слід ототожнювати з моделлю "
            "однієї ТВЗ з відбивними межами на кожній її грані."
        ),
    }
    (run_dir / "pidsumok_rozrakhunku.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_power_map(pin_df, p, run_dir / "karta_vidnosnoi_potuzhnosti_klastera_7_tvz.png")
    return summary


def safe_name(value: str) -> str:
    """Зробити безпечну частину назви папки."""
    return "".join(char if char.isalnum() or char in "-_" else "_" for char in value).strip("_-")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=PROJECT_ROOT / "cluster_7_parameters.json",
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=PROJECT_ROOT / "results" / "03_cluster_7_2d" / "runs",
    )
    parser.add_argument("--build-only", action="store_true", help="лише побудувати геометрію")
    parser.add_argument(
        "--process-existing", type=Path,
        help="повторно обробити готовий statepoint без нового розрахунку",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.process_existing:
        run_dir = args.process_existing
        p = load_parameters(run_dir / "vykorystani_parametry.json")
        _, metadata = build_model(p)
        summary = process_results(run_dir / f"statepoint.{p['simulation']['batches']}.h5", metadata, p, run_dir)
        print(f"Результати повторно оброблено: {run_dir.resolve()}")
        print(f"k_inf = {summary['k_inf']:.6f} +/- {summary['похибка_k_inf']:.6f}")
        return
    p = load_parameters(args.config)
    label = safe_name(p.get("run_name", "klaster_7_tvz"))
    run_dir = args.output_root / f"{datetime.now():%Y%m%d_%H%M%S}_{label}"
    run_dir.mkdir(parents=True)
    (run_dir / "vykorystani_parametry.json").write_text(
        json.dumps(p, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    model, metadata = build_model(p)
    position_table = pd.DataFrame(metadata["records"]).drop(columns=["fuel_cell_id"], errors="ignore")
    ASSEMBLY.ukrainian_table(position_table).to_csv(run_dir / "karta_pozytsii_klastera_7_tvz.csv", index=False)
    model.export_to_model_xml(path=run_dir / "model.xml")
    model.plot_geometry(cwd=run_dir)
    if args.build_only:
        print(f"Геометрію кластера побудовано: {run_dir.resolve()}")
        return
    statepoint = Path(model.run(cwd=run_dir, threads=p["simulation"]["threads"]))
    summary = process_results(statepoint, metadata, p, run_dir)
    print(f"k_inf = {summary['k_inf']:.6f} +/- {summary['похибка_k_inf']:.6f}")
    print(f"Папка результатів: {run_dir.resolve()}")


if __name__ == "__main__":
    main()
