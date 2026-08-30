$project = (wsl -d Ubuntu -- wslpath -a $PSScriptRoot).Trim()
wsl -d Ubuntu -- bash -lic "conda activate openmc-env && cd '$project' && exec bash"
