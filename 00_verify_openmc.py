"""Перевірка встановлення OpenMC і бібліотеки ENDF/B-VIII.0.

Запуск з Ubuntu/WSL:
    conda activate openmc-env
    cd /mnt/c/Work/universitet/task_openmc
    python 00_verify_openmc.py
"""

from pathlib import Path
import os

import openmc


def check_nuclear_data() -> Path:
    """Перевірити змінну середовища та наявність потрібних даних."""
    value = os.environ.get("OPENMC_CROSS_SECTIONS")
    if not value:
        raise RuntimeError("Не задано змінну OPENMC_CROSS_SECTIONS")

    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(path)

    library = openmc.data.DataLibrary.from_xml(path)
    available = set()
    for item in library:
        available.update(item["materials"])

    required = {
        "H1",
        "O16",
        "U235",
        "U238",
        "Zr90",
        "Nb93",
        "Gd155",
        "Gd157",
        "c_H_in_H2O",
    }
    missing = sorted(required - available)
    if missing:
        raise RuntimeError(f"У бібліотеці немає даних: {missing}")

    print(f"OpenMC: {openmc.__version__}")
    print(f"Бібліотека: {path}")
    print(f"Записів у cross_sections.xml: {len(library)}")
    print("Потрібні дані U, Gd, Zr, Nb, H, O та S(a,b) для води знайдено.")
    return path


def build_test_model() -> openmc.Model:
    """Створити дуже малу тестову модель елементарної комірки."""
    fuel = openmc.Material(name="UO2 3.7%")
    fuel.add_element("U", 1.0, enrichment=3.7)
    fuel.add_element("O", 2.0)
    fuel.set_density("g/cm3", 10.4)

    water = openmc.Material(name="Light water")
    water.add_element("H", 2.0)
    water.add_element("O", 1.0)
    water.set_density("g/cm3", 0.997)
    water.add_s_alpha_beta("c_H_in_H2O")

    fuel_surface = openmc.ZCylinder(r=0.386)
    half_pitch = 1.275 / 2.0
    xmin = openmc.XPlane(x0=-half_pitch, boundary_type="reflective")
    xmax = openmc.XPlane(x0=half_pitch, boundary_type="reflective")
    ymin = openmc.YPlane(y0=-half_pitch, boundary_type="reflective")
    ymax = openmc.YPlane(y0=half_pitch, boundary_type="reflective")
    zmin = openmc.ZPlane(z0=-0.5, boundary_type="reflective")
    zmax = openmc.ZPlane(z0=0.5, boundary_type="reflective")
    box = +xmin & -xmax & +ymin & -ymax & +zmin & -zmax

    fuel_cell = openmc.Cell(
        name="Fuel",
        fill=fuel,
        region=box & -fuel_surface,
    )
    moderator_cell = openmc.Cell(
        name="Moderator",
        fill=water,
        region=box & +fuel_surface,
    )

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = 15
    settings.inactive = 5
    settings.particles = 1_000
    settings.seed = 20260812
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            [-0.38, -0.38, -0.49],
            [0.38, 0.38, 0.49],
        ),
        constraints={"fissionable": True},
    )

    geometry = openmc.Geometry([fuel_cell, moderator_cell])
    return openmc.Model(
        geometry=geometry,
        materials=openmc.Materials([fuel, water]),
        settings=settings,
    )


def main() -> None:
    check_nuclear_data()
    output = Path("results") / "installation_check"
    output.mkdir(parents=True, exist_ok=True)

    model = build_test_model()
    statepoint = model.run(cwd=output, threads=2)
    with openmc.StatePoint(statepoint) as sp:
        print(f"Тестовий розрахунок завершено: keff = {sp.keff}")
    print(f"Результати: {output.resolve()}")


if __name__ == "__main__":
    main()
