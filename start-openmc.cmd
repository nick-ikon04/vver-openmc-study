@echo off
setlocal
chcp 65001 >nul
title OpenMC environment
for /f "usebackq delims=" %%I in (`wsl -d Ubuntu -- wslpath -a "%~dp0"`) do set "PROJECT=%%I"
if not defined PROJECT (
  echo Не вдалося перетворити шлях проєкту для WSL.
  exit /b 1
)
wsl -d Ubuntu -- bash -lic "conda activate openmc-env && cd '%PROJECT%' && exec bash"
exit /b %ERRORLEVEL%
