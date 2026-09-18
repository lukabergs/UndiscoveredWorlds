param(
    [string]$BuildDir = 'out/build/x64-Debug',
    [string]$Output = 'out/benchmarks/geology',
    [int[]]$Seeds = @(12345, 24680, 424242),
    [ValidateSet(128, 256, 512)][int[]]$Widths = @(128, 256),
    [ValidateSet(128, 256, 512, 1024, 2048)][int]$ReferenceWidth = 512
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
Push-Location $repo
try {
    $pipeline = Join-Path $BuildDir 'bin/tectonic_pipeline.exe'
    if (!(Test-Path -LiteralPath $pipeline)) { throw "Build tectonic_pipeline first: $pipeline" }
    $bundles = @()
    foreach ($width in $Widths) {
        foreach ($seed in $Seeds) {
            $bundle = Join-Path $Output "seed-$seed-$width"
            if (Test-Path -LiteralPath $bundle) { throw "Output already exists; choose another -Output: $bundle" }
            Write-Host "Tectonics seed=$seed grid=$width x $($width / 2)"
            & $pipeline inspect --final --seed $seed --width $width --height ($width / 2) `
                --plates 10 --cycles 2 --cycle-step-limit 600 --erosion-period 60 `
                --aggregation-overlap-abs ([Math]::Max(64, [int][Math]::Floor($width * $width / 2000))) `
                --aggregation-overlap-rel 0.20 --folding-ratio 0.08 --bundle-output $bundle `
                --output (Join-Path $bundle 'summary.json')
            if ($LASTEXITCODE -ne 0) { throw "Tectonic export failed for seed $seed, width $width" }
            $bundles += $bundle
        }
    }
    & uv run scripts/benchmarks/benchmark-geology.py --bundle @bundles `
        --width $ReferenceWidth --output (Join-Path $Output 'report.json')
    if ($LASTEXITCODE -ne 0) { throw 'Geology benchmark failed; inspect report.json for invariant violations' }
}
finally { Pop-Location }
