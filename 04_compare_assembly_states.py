"""Зібрати україномовну таблицю та графік п'яти станів однієї ТВЗ."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT_ROOT = (
    PROJECT_ROOT / "results" / "02_assembly_331_2d" / "temperature_states"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "results" / "02_assembly_331_2d" / "porivniannia_staniv"
)

STATE_NAMES = {
    "tvz331_kholodnyi_T300_B1300": "Холодний стан",
    "tvz331_hariachyi_bez_potuzhnosti_T575_B1300": "Гарячий стан без потужності",
    "tvz331_robochyi_Tfuel1027_Twater575_B1300": "Робочий стан на потужності",
    "tvz331_robochyi_Tfuel1027_Twater575_B0": "Робочий стан без бору",
    "tvz331_pidvyshchena_Tfuel1200_Twater575_B1300": "Підвищена температура палива",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root", type=Path, default=DEFAULT_INPUT_ROOT,
        help="папка із завершеними розрахунками п'яти станів",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="папка для підсумкової таблиці та графіка",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.input_root.is_dir():
        raise FileNotFoundError(
            f"Папку з результатами станів не знайдено: {args.input_root}"
        )
    rows = []
    for run_dir in args.input_root.iterdir():
        summary_path = run_dir / "pidsumok_rozrakhunku.json"
        parameters_path = run_dir / "vykorystani_parametry.json"
        if not summary_path.exists() or not parameters_path.exists():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
        run_name = parameters["run_name"]
        if run_name not in STATE_NAMES:
            continue
        rows.append({
            "стан": STATE_NAMES[run_name],
            "температура_палива_K": parameters["fuel"]["temperature_K"],
            "температура_води_K": parameters["moderator"]["temperature_K"],
            "густина_води_г_см3": parameters["moderator"]["density_g_cm3"],
            "температура_оболонки_K": parameters["cladding"]["temperature_K"],
            "бор_ppm": parameters["moderator"]["boron_ppm"],
            "k_inf": summary["k_inf"],
            "похибка_k_inf": summary["похибка_k_inf"],
            "K_pin": summary["K_pin"],
            "похибка_K_pin_наближено": summary["похибка_K_pin_наближено"],
            "номер_найбільш_навантаженого_твела": summary["найбільш_навантажений_твел"]["номер_позиції"],
            "середня_потужність_звичайних_твелів": summary["середня_відносна_потужність_звичайних_твелів"],
            "середня_потужність_твегів": summary["середня_відносна_потужність_твегів"],
        })
    table = pd.DataFrame(rows)
    if len(table) != 5:
        raise RuntimeError(f"Очікувалося п'ять завершених станів, знайдено {len(table)}")
    order = list(STATE_NAMES.values())
    table["порядок"] = table["стан"].map({name: index for index, name in enumerate(order)})
    table = table.sort_values("порядок").drop(columns="порядок")
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "porivniannia_piaty_staniv_tvz.csv", index=False)

    fig, ax = plt.subplots(figsize=(12, 6.8))
    labels = [
        "Холодний\n300 K, B=1300",
        "Гарячий без\nпотужності, 575 K",
        "Робочий\n1027 K, B=1300",
        "Робочий\n1027 K, B=0",
        "Паливо 1200 K\nB=1300",
    ]
    ax.errorbar(
        range(5), table["k_inf"], yerr=table["похибка_k_inf"],
        fmt="o-", color="#2563eb", capsize=5, linewidth=2, markersize=7,
    )
    for index, value in enumerate(table["k_inf"]):
        ax.annotate(f"{value:.5f}", (index, value), xytext=(0, 10),
                    textcoords="offset points", ha="center", fontsize=9)
    ax.axhline(1.0, color="#dc2626", linestyle="--", linewidth=1.2, label="Критичний рівень k=1")
    ax.set_xticks(range(5), labels)
    ax.set_ylabel(r"Коефіцієнт розмноження $k_\infty$")
    ax.set_title("Порівняння п'яти температурно-борних станів ТВЗ")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output / "porivniannia_k_inf_piaty_staniv_tvz.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
