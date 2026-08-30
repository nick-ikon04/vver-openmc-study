# VVER OpenMC Study

OpenMC models for neutronics analysis of a VVER fuel pin, a VVER-1000
331-pin fuel assembly, and a seven-assembly cluster.

## Contents

- `00_verify_openmc.py` checks the OpenMC installation and nuclear data.
- `01_pin_compare.py` compares UO2 and UGD fuel-pin models.
- `02_assembly_331_2d.py` models a 331-pin VVER fuel assembly.
- `03_cluster_7_2d.py` models a seven-assembly cluster.
- `04_compare_assembly_states.py` compares operating states.
- `parameters/` contains versioned JSON input configurations.

Generated HDF5 files, calculation results, reports, archives, and local nuclear
data are intentionally excluded from Git. They are stored separately in the
project's Google Drive archive.

## Environment

The launch scripts expect WSL Ubuntu and a Conda environment named
`openmc-env`. The environment must provide OpenMC, NumPy, pandas, and
Matplotlib, together with a configured OpenMC nuclear-data library.

Run `start-openmc.ps1` or `start-openmc.cmd` to open the project environment.
Use `start-jupyter.cmd` to start JupyterLab.
