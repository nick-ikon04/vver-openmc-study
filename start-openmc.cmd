@echo off
title OpenMC environment
for /f "usebackq delims=" %%I in (`wsl -d Ubuntu -- wslpath -a "%~dp0"`) do set "PROJECT=%%I"
wsl -d Ubuntu -- bash -lic "conda activate openmc-env && cd '%PROJECT%' && exec bash"
