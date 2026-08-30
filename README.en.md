# VVER OpenMC Study

[Українська](README.md) | [English](README.en.md)

Parameterized OpenMC models for neutronics analysis of a fuel pin, a VVER-1000
331-position fuel assembly, and a seven-assembly cluster.

## Features

- OpenMC installation and nuclear-data verification;
- comparison of UO₂ and UO₂–Gd₂O₃ fuel-pin cells;
- a 2D VVER-1000 fuel-assembly model with 331 lattice positions;
- a seven-assembly cluster model;
- comparison of five temperature and boron states;
- CSV tables, JSON summaries, and PNG plots.

## Project files

| File | Purpose |
| --- | --- |
| `00_verify_openmc.py` | Checks OpenMC and nuclear data, then runs a small test model. |
| `01_pin_compare.py` | Compares UO₂ and UO₂–Gd₂O₃ fuel-pin cells. |
| `02_assembly_331_2d.py` | Models one 331-position fuel assembly. |
| `03_cluster_7_2d.py` | Models a seven-assembly cluster. |
| `04_compare_assembly_states.py` | Compares five completed assembly states. |
| `parameters/assembly_states/` | Stores versioned state configurations. |
| `environment.yml` | Recreates the `openmc-env` Conda environment. |

More detailed model documentation is available for the
[fuel-pin comparison](README_01_pin_compare.md) and the
[331-position assembly](README_02_assembly_331_2d.md).

## Report and example results

- The [five-page Ukrainian technical report](docs/zvit_openmc_tvel_tvz_7tvz_20260823_uk.docx)
  summarizes the fuel-pin, 331-position assembly, seven-assembly cluster, and
  five temperature and boron state calculations.
- The [verified example results](examples/README.md) provide compact CSV, JSON,
  and PNG files without large HDF5 data, generated XML, or OpenMC logs.

## Requirements

- Windows 10 or 11 with WSL 2;
- Ubuntu on WSL;
- Conda or Mamba;
- an OpenMC-format nuclear-data library.

Create the environment:

```bash
conda env create -f environment.yml
conda activate openmc-env
```

Set `OPENMC_CROSS_SECTIONS` to the `cross_sections.xml` file of the installed
nuclear-data library:

```bash
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
```

The nuclear-data library is not included because of its size.

## Launching on Windows

- `start-openmc.ps1` or `start-openmc.cmd` opens the project environment;
- `start-jupyter.cmd` starts JupyterLab.

The launchers derive the repository path automatically and support paths that
contain spaces. They do not depend on an old computer-specific path.

Verify the environment and run a small OpenMC test:

```bash
python 00_verify_openmc.py
```

Run the fuel-pin comparison, one assembly, or the seven-assembly cluster:

```bash
python 01_pin_compare.py
python 02_assembly_331_2d.py
python 03_cluster_7_2d.py
```

Build geometry without particle transport:

```bash
python 02_assembly_331_2d.py --build-only
python 03_cluster_7_2d.py --build-only
```

Use `--config` to select another configuration and `--output-root` to select
an output directory. Default paths are derived from the script location, so
the commands can be run from any working directory.

## Comparing assembly states

Run all five state configurations into one directory, then build the summary:

```bash
for config in parameters/assembly_states/*.json; do
  python 02_assembly_331_2d.py \
    --config "$config" \
    --output-root results/02_assembly_331_2d/temperature_states
done
python 04_compare_assembly_states.py
```

## Files excluded from Git

`.gitignore` excludes complete OpenMC run directories, HDF5 and generated XML
files, most reports, archives, caches, virtual environments, and local settings.
Only selected compact examples and one published technical report are tracked.
Store large output data separately, for example in Google Drive.

Do not publish a personal `cross_sections.xml`, `.env` files, or local paths.
