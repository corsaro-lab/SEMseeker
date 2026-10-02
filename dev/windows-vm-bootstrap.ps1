# ---------------------------------------------------------------------------
# Windows VM - bootstrap (run FIRST, inside Windows, in an elevated PowerShell)
#
#   Right-click Start -> "Terminal (Admin)", then:
#     Set-ExecutionPolicy -Scope Process Bypass -Force
#     \\Mac\Home\Documents\Progetti\RICERCA\semseeker\SEMseeker\dev\windows-vm-bootstrap.ps1
#
# Installs the toolchain, then hands over to dev/windows-vm-setup.R for the
# package dependencies.
#
# ASCII ONLY, deliberately. Windows PowerShell 5.1 reads a file without a BOM as
# cp1252, so a UTF-8 em dash arrives as three bytes whose last one is a curly
# quote - which closes a string early and produces parse errors pointing at
# lines that are perfectly fine. Keep this file to plain ASCII.
#
# The R version is pinned to what CI uses ("release"). Matching it matters: the
# point of this VM is to reproduce the Windows job, and a different R would
# reproduce something else.
#
# On Apple Silicon, Windows 11 is ARM and R is x64, so R runs under the OS
# emulation layer. That is expected and fine for what Windows actually breaks:
# path separators, case-insensitive file systems, file locking, encoding,
# `parallel` without fork. It is NOT equivalent for an architecture-specific
# defect; only the CI can settle those.
# ---------------------------------------------------------------------------

$ErrorActionPreference = "Stop"

$RVersion = "4.6.1"
$RBin     = "C:\Program Files\R\R-$RVersion\bin"

Write-Host "== R $RVersion ==" -ForegroundColor Cyan
winget install --id RProject.R --version $RVersion --silent --accept-package-agreements --accept-source-agreements

Write-Host "== Rtools ==" -ForegroundColor Cyan
winget install --id RProject.Rtools --silent --accept-package-agreements --accept-source-agreements

Write-Host "== Git (needed by devtools) ==" -ForegroundColor Cyan
winget install --id Git.Git --silent --accept-package-agreements --accept-source-agreements

# Do not assume where R landed. On Windows ARM winget installs native aarch64
# builds (Rtools45 arrives as rtools45-aarch64-*.exe), and the install root is
# not necessarily "C:\Program Files\R\R-<version>". Find Rscript.exe instead of
# predicting it.
$roots = @("C:\Program Files\R", "C:\Program Files (x86)\R",
           "$env:LOCALAPPDATA\Programs\R") | Where-Object { Test-Path $_ }

$rscript = $null
if ($roots) {
    $rscript = Get-ChildItem -Path $roots -Filter Rscript.exe -Recurse -ErrorAction SilentlyContinue |
               Sort-Object FullName -Descending | Select-Object -First 1 -ExpandProperty FullName
}

if (-not $rscript) {
    Write-Warning "Rscript.exe not found under $($roots -join ', '). Locate it and add its folder to PATH by hand."
}
else {
    $RBin = Split-Path $rscript
    $env:Path = "$RBin;" + $env:Path
    $userPath = [Environment]::GetEnvironmentVariable("Path", "User")
    [Environment]::SetEnvironmentVariable("Path", "$RBin;$userPath", "User")
    Write-Host "R found at $RBin and added to PATH for this user." -ForegroundColor Green

    # Which R this is decides what the VM is worth as a test bench: an aarch64
    # build is a THIRD system, not the x64 runner CI uses. Print it, so the
    # answer is on the screen instead of assumed.
    & $rscript -e "cat(R.version.string, '/', R.version`$platform, '\n')"
}

Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "  net use Z: \\Mac\Home /persistent:yes"
Write-Host "  cd Z:\Documents\Progetti\RICERCA\semseeker\SEMseeker"
Write-Host "  Rscript dev\windows-vm-setup.R"
Write-Host ""
Write-Host "The drive mapping is not cosmetic: cmd.exe refuses a UNC path as a" -ForegroundColor Yellow
Write-Host "working directory, so running from \\Mac\Home directly will fail." -ForegroundColor Yellow
