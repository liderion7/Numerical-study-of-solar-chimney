# run_n128_sweep.ps1
# -----------------------------------------------------------------------------
# Batch driver for an angle sweep: run_cd_case.jl -> task_c.jl -> task_d.jl.
#
# RA decides which task line the products belong to:
#   RA = 1000 -> low-Ra validation   (results root: ..\validation_lowRa\data)
#   RA = 1e7  -> production Task C/D (results root: .\data)
# Switch RA (and $resultsRoot) below to run the other task line.
#
# All paths are relative to this script (main_task_cd/), so the whole repository
# can be moved without editing them. The Julia path is machine-specific.
# -----------------------------------------------------------------------------
$ErrorActionPreference = 'Continue'
$julia = 'C:\Users\18497\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe'
$here = $PSScriptRoot
$angles = 0, 15, 30, 45, 60, 75

$ra = '1000'
$resultsRoot = if ($ra -eq '1000') {
    Join-Path $here '..\validation_lowRa\data'
} else {
    Join-Path $here 'data'
}
$env:RESULTS_ROOT = $resultsRoot

foreach ($a in $angles) {
    Write-Output "=== START theta=$a N=128 ==="
    $env:THETA = "$a"
    $env:T_END = '550'
    $env:RA = $ra
    $env:PR = '0.71'
    $env:N_WALL = '128'
    $env:SAMPLE_DT = '0.2'
    $env:MONITOR_DT = '5.0'

    & $julia (Join-Path $here 'run_cd_case.jl')
    if ($LASTEXITCODE -ne 0) {
        Write-Output "CFD FAILED for theta=$a (exit $LASTEXITCODE)"
        continue
    }

    $case = Join-Path $resultsRoot "${a}degree"
    $env:CASE_DIR = $case
    & $julia (Join-Path $here 'task_c.jl')
    if ($LASTEXITCODE -ne 0) {
        Write-Output "task_c FAILED for theta=$a (exit $LASTEXITCODE)"
        continue
    }

    $clean = Join-Path $case 'task_c\task_c_clean_vbar.csv'
    if (Test-Path $clean) {
        & $julia (Join-Path $here 'task_d.jl')
        if ($LASTEXITCODE -ne 0) {
            Write-Output "task_d FAILED for theta=$a (exit $LASTEXITCODE)"
        } else {
            Write-Output "task_d OK for theta=$a"
        }
    } else {
        Write-Output "Task C not ready for theta=$a; skipping task_d"
    }
    Write-Output "=== DONE theta=$a ==="
}

Write-Output 'ALL_N128_SWEEP_DONE'
