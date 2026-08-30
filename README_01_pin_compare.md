# Скрипт 01 - порівняння UO2 та UO2-Gd2O3

Усі змінювані параметри містяться у файлі `pin_parameters.json`. Температуру,
збагачення і вміст гадолінію можна змінювати без редагування Python-коду.

## Значення концентрацій

- `uo2_u235_enrichment_wt_percent` - масова частка U-235 в урані звичайного твела;
- `ugd_u235_enrichment_wt_percent` - масова частка U-235 в урані гадолінієвого палива;
- `gd2o3_wt_percent` - масова частка Gd2O3 у всій суміші UO2-Gd2O3.

Наприклад, `gd2o3_wt_percent: 4.0` означає суміш із 96 мас.% UO2 та
4 мас.% Gd2O3. OpenMC автоматично використовує природний ізотопний склад
гадолінію. Температури задаються в K, густини - у г/см3, розміри - у см.

## Запуск

```bash
python 01_pin_compare.py --config pin_parameters.json
```

Назву серії можна змінити через `run_name` у JSON або командним параметром:

```bash
python 01_pin_compare.py --config pin_parameters.json --run-name T900_Gd6
```

## Результати

Кожен запуск створює окрему папку в `results/01_pin_compare/runs/` і не
перезаписує попередні результати. Файл `vykorystani_parametry.json` зберігає
точний знімок параметрів запуску.

Основні вихідні файли мають змістовні назви:

- `heometriia_komirky_uo2.png`;
- `heometriia_komirky_uo2_gd2o3.png`;
- `porivniannia_k_inf_komirok_uo2_ta_uo2_gd2o3.png`;
- `intehralni_reaktsii_komirky.csv`;
- `reaktsii_za_oblastiamy_komirky.csv`;
- `pohlynannia_izotopamy_hadoliniiu.csv`;
- `pidsumok_rozrakhunku.json`.

Допустимий діапазон температур у поточній бібліотеці даних - 294-1200 K.
Скрипт перевіряє порядок радіусів, концентрації та кількість активних і
неактивних поколінь до запуску тривалого розрахунку.
