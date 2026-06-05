# mapseq2-py installer for Windows.
#
# This script:
#   1. Verifies that WSL2 + Ubuntu is installed (offers to install if missing).
#   2. Delegates the actual install to install.sh inside Ubuntu.
#
# Run from PowerShell (does NOT need Administrator unless WSL needs installing):
#   .\install.ps1
#
# To install into a custom-named conda env:
#   $env:ENV_NAME = "foo"; .\install.ps1

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition

function Write-Info  ($msg) { Write-Host "[install.ps1] $msg" -ForegroundColor Cyan }
function Write-WarnX ($msg) { Write-Host "[install.ps1 WARN] $msg" -ForegroundColor Yellow }
function Write-Err   ($msg) { Write-Host "[install.ps1 ERROR] $msg" -ForegroundColor Red }

# ---------- 1. WSL ----------------------------------------------------------
$wslPath = (Get-Command wsl.exe -ErrorAction SilentlyContinue)
if (-not $wslPath) {
    Write-Err "WSL is not installed."
    Write-Host "  Open an Administrator PowerShell and run:  wsl --install -d Ubuntu"
    Write-Host "  Reboot, finish the Ubuntu first-time setup, then re-run this script."
    exit 1
}

$distros = (wsl.exe --list --quiet) -replace "`0", ""
if (-not ($distros -match "Ubuntu")) {
    Write-Err "No Ubuntu distro found in WSL."
    Write-Host "  Run (Administrator):  wsl --install -d Ubuntu"
    exit 1
}

# Pick the default Ubuntu distro
$distro = ($distros -split "`r?`n" | Where-Object { $_ -match "Ubuntu" } | Select-Object -First 1).Trim()
Write-Info "Using WSL distro: $distro"

# ---------- 2. Convert this script's path to a WSL path --------------------
# C:\Users\rocke\MAPseq\mapseq2-py -> /mnt/c/Users/rocke/MAPseq/mapseq2-py
$driveLetter = $ScriptDir.Substring(0, 1).ToLower()
$wslProjPath = "/mnt/$driveLetter" + ($ScriptDir.Substring(2) -replace '\\', '/')
Write-Info "Project dir (WSL view): $wslProjPath"

# ---------- 3. Run install.sh inside WSL -----------------------------------
$envName = if ($env:ENV_NAME) { $env:ENV_NAME } else { "mapseq" }
Write-Info "Delegating to install.sh inside $distro (ENV_NAME=$envName)..."

wsl.exe -d $distro -- bash -lc "ENV_NAME='$envName' bash '$wslProjPath/install.sh'"
if ($LASTEXITCODE -ne 0) {
    Write-Err "install.sh failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

Write-Info "DONE."
Write-Host @"

To use mapseq2 from now on, open a WSL shell and run:
    conda activate $envName
    mapseq2 --help

See MANUAL.md for the full workflow.
"@
