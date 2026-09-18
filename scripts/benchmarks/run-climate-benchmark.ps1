[CmdletBinding()]
param(
    [ValidateRange(64, 2048)][int]$Resolution = 128,
    [ValidateRange(64, 2048)][int]$WindMapResolution = 2048,
    [ValidateRange(8, 1024)][int]$ClimateResolution = 64,
    [ValidateRange(0, 2048)][int]$HydrologyResolution = 0,
    [ValidateRange(1, 100)][int]$CouplingIterations = 8,
    [long]$Seed = 20260906,
    [string]$Maps = "koppen,temp,precip,precip-tif,jan-s-wind-speed,jul-s-wind-speed,jan-u-wind-speed,jul-u-wind-speed,jan-s-wind-lic,jan-s-wind-part,jul-s-wind-lic,jul-s-wind-part,jan-u-wind-lic,jan-u-wind-part,jul-u-wind-lic,jul-u-wind-part,jan-sst,jul-sst,jan-ocean-current,jul-ocean-current",
    [ValidateSet("calibrated", "radiative")][string]$TemperatureMode = "calibrated",
    [ValidateSet("baseline", "transient")][string]$ClimateMode = "baseline",
    [ValidateSet("ffsl", "mpdata")][string]$Transport = "ffsl",
    [ValidateSet("mixed-layer", "legacy")][string]$TropicalClosure = "legacy",
    [ValidateSet(1, 2, 3, 4, 6, 8, 12, 24)][int]$OceanStepHours = 24,
    [ValidateSet("reuse", "recompute")][string]$OceanCirculationCache = "reuse",
    [ValidateSet("krylov", "relaxation")][string]$OceanCirculationSolver = "krylov",
    [ValidateSet("gather", "legacy")][string]$OceanHeatKernel = "gather",
    [ValidateRange(0, 32)][int]$OceanHeatWorkers = 0,
    [ValidateSet("legacy", "face-drag", "face-forcing")][string]$OceanCoastalScheme = "face-forcing",
    [string]$OceanCapture = "",
    [ValidateSet("zonal", "jacobi", "legacy")][string]$AtmosphereSolver = "zonal",
    [string]$AtmosphereCapture = "",
    [string]$Executable = "out/build/x64-Debug/bin/UndiscoveredWorlds.exe"
)

$ErrorActionPreference = "Stop"
if ($Resolution % 2 -ne 0 -or $ClimateResolution % 4 -ne 0) {
    throw "Resolution must be even; ClimateResolution must be a multiple of four."
}
if ($WindMapResolution % 2 -ne 0) { throw "WindMapResolution must be even." }
if ($PSBoundParameters.ContainsKey("HydrologyResolution") -and ($HydrologyResolution -lt 8 -or $HydrologyResolution % 4 -ne 0)) {
    throw "HydrologyResolution must be a multiple of four and at least eight."
}
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "../.."))
Push-Location -LiteralPath $repoRoot
try {
    if (-not (Test-Path -LiteralPath $Executable)) { throw "Build the executable first: $Executable" }
    $logRoot = Join-Path $repoRoot "runs/logs/climate"
    New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
    $stamp = Get-Date -Format "yyyyMMdd-HHmmss-fff"
    $pendingLog = Join-Path $logRoot "pending-$stamp.log"
    $information = "${Resolution}x$($Resolution / 2); seed $Seed; climate $ClimateResolution; coupling limit $CouplingIterations; temperature $TemperatureMode; mode $ClimateMode; transport $Transport; tropical closure $TropicalClosure; maps $Maps"
    $arguments = @("--earth-climate-benchmark", "--seed", "$Seed", "--resolution", "$Resolution",
        "--climate-resolution", "$ClimateResolution", "--climate-temperature-mode", $TemperatureMode,
        "--climate-coupling-iterations", "$CouplingIterations",
        "--climate-mode", $ClimateMode, "--climate-transport", $Transport,
        "--climate-tropical-closure", $TropicalClosure,
        "--map", $Maps, "--no-rivers", "--no-lakes", "--no-deltas", "--benchmark-info", $information)
    # Optional arguments allow replay with older archived executables.
    if ($PSBoundParameters.ContainsKey("WindMapResolution")) {
        $arguments += @("--wind-map-resolution", "$WindMapResolution")
        $information += "; wind map width $WindMapResolution"
    }
    if ($PSBoundParameters.ContainsKey("HydrologyResolution")) {
        $arguments += @("--climate-hydrology-resolution", "$HydrologyResolution")
        $information += "; hydrology $HydrologyResolution"
    }
    if ($PSBoundParameters.ContainsKey("OceanStepHours")) {
        $arguments += @("--climate-ocean-step-hours", "$OceanStepHours")
        $information += "; ocean step $OceanStepHours hours"
    }
    if ($PSBoundParameters.ContainsKey("OceanCirculationCache")) {
        $arguments += @("--climate-ocean-circulation-cache", $OceanCirculationCache)
        $information += "; ocean circulation $OceanCirculationCache"
    }
    if ($PSBoundParameters.ContainsKey("OceanCirculationSolver")) {
        $arguments += @("--climate-ocean-circulation-solver", $OceanCirculationSolver)
        $information += "; ocean solver $OceanCirculationSolver"
    }
    if ($PSBoundParameters.ContainsKey("AtmosphereSolver")) {
        $arguments += @("--climate-atmosphere-solver", $AtmosphereSolver)
        $information += "; atmosphere solver $AtmosphereSolver"
    }
    if ($PSBoundParameters.ContainsKey("OceanHeatKernel")) {
        $arguments += @("--climate-ocean-heat-kernel", $OceanHeatKernel)
        $information += "; ocean heat kernel $OceanHeatKernel"
    }
    if ($PSBoundParameters.ContainsKey("OceanHeatWorkers")) {
        $arguments += @("--climate-ocean-heat-workers", "$OceanHeatWorkers")
        $information += "; ocean heat workers $OceanHeatWorkers"
    }
    if ($PSBoundParameters.ContainsKey("OceanCoastalScheme")) {
        $arguments += @("--climate-ocean-coastal-scheme", $OceanCoastalScheme)
        $information += "; ocean coastal scheme $OceanCoastalScheme"
    }
    if ($OceanCapture) {
        $arguments += @("--climate-ocean-capture", $OceanCapture)
        $information += "; ocean capture $OceanCapture"
    }
    if ($AtmosphereCapture) {
        $arguments += @("--climate-atmosphere-capture", $AtmosphereCapture)
        $information += "; atmosphere capture $AtmosphereCapture"
    }
    $arguments[[Array]::IndexOf($arguments, "--benchmark-info") + 1] = $information
    Write-Host "Benchmark: $information"
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    # Native stderr can contain diagnostic warnings on an otherwise successful
    # run (especially under Windows PowerShell 5). The exit code is authoritative.
    $logWriter = [System.IO.StreamWriter]::new($pendingLog, $false, [System.Text.UTF8Encoding]::new($false))
    $logWriter.AutoFlush = $true
    $logWriter.WriteLine("Executable SHA256: $((Get-FileHash -LiteralPath $Executable -Algorithm SHA256).Hash)")
    $sourceHead = & git rev-parse HEAD
    $sourceDirty = [bool](& git status --porcelain --untracked-files=normal)
    $logWriter.WriteLine("Source HEAD: $sourceHead; working tree modified: $sourceDirty")
    try {
        $ErrorActionPreference = "Continue"
        & $Executable @arguments 2>&1 | ForEach-Object {
            $line = $_.ToString()
            $logWriter.WriteLine($line)
            Write-Host $line
        }
        $benchmarkExit = $LASTEXITCODE
    }
    finally {
        $logWriter.Dispose()
        $ErrorActionPreference = "Stop"
    }
    $timer.Stop()
    if ($benchmarkExit -ne 0) { throw "Benchmark exited $benchmarkExit; log: $pendingLog" }
    # Read the ID printed by this process, rather than guessing the registry's next ID.
    $idLine = Select-String -LiteralPath $pendingLog -Pattern 'Climate benchmark run ID: (\d+)' | Select-Object -Last 1
    if (-not $idLine) { throw "No completed run ID in $pendingLog" }
    $runId = $idLine.Matches[0].Groups[1].Value
    $finalLog = Join-Path $logRoot "$runId.log"
    if (Test-Path -LiteralPath $finalLog) { throw "Refusing to replace existing log: $finalLog" }
    Move-Item -LiteralPath $pendingLog -Destination $finalLog
    $elapsedSeconds = $timer.Elapsed.TotalSeconds.ToString([System.Globalization.CultureInfo]::InvariantCulture)
    & uv run --offline python scripts/benchmarks/write-benchmark-report.py --run-id $runId --elapsed-seconds $elapsedSeconds
    if ($LASTEXITCODE -ne 0) { throw "Run $runId finished, but report generation failed." }
    Write-Host "Run $runId completed in $($timer.Elapsed). Log: $finalLog"
}
finally {
    Pop-Location
}
