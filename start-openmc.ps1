$project = wsl -d Ubuntu -- wslpath -a $PSScriptRoot
if ($LASTEXITCODE -ne 0) {
    throw "Не вдалося перетворити шлях проєкту для WSL."
}
$project = $project.Trim()
wsl -d Ubuntu -- bash -lic "conda activate openmc-env && cd '$project' && exec bash"
exit $LASTEXITCODE
