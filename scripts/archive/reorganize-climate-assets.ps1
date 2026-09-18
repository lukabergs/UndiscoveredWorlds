[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))
$extraRoot = [System.IO.Path]::GetFullPath((Join-Path $repoRoot "extra"))
$extraPrefix = $extraRoot.TrimEnd([System.IO.Path]::DirectorySeparatorChar) +
    [System.IO.Path]::DirectorySeparatorChar

function Move-Safe {
    param(
        [Parameter(Mandatory = $true)][string]$Source,
        [Parameter(Mandatory = $true)][string]$Destination
    )

    $sourcePath = [System.IO.Path]::GetFullPath((Join-Path $repoRoot $Source))
    $destinationPath = [System.IO.Path]::GetFullPath((Join-Path $repoRoot $Destination))
    if (-not $sourcePath.StartsWith($extraPrefix, [System.StringComparison]::OrdinalIgnoreCase) -or
        -not $destinationPath.StartsWith($extraPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Asset move escaped the extra directory: $Source -> $Destination"
    }
    if (-not (Test-Path -LiteralPath $sourcePath)) {
        return
    }
    if (Test-Path -LiteralPath $destinationPath) {
        throw "Destination already exists; refusing to overwrite: $destinationPath"
    }

    $parent = Split-Path -Parent $destinationPath
    New-Item -ItemType Directory -Force -Path $parent | Out-Null
    Move-Item -LiteralPath $sourcePath -Destination $destinationPath
    Write-Host "Moved $Source -> $Destination"
}

Move-Safe "extra\img\earth\benchmark" "extra\climate\benchmarks\maps"
Move-Safe "extra\climate\benchmarks\maps\ref" "extra\reference\climate\previews"
Move-Safe "extra\validation\runs" "extra\climate\benchmarks\runs"

$validationRoot = Join-Path $extraRoot "validation"
if (Test-Path -LiteralPath $validationRoot) {
    Get-ChildItem -LiteralPath $validationRoot | ForEach-Object {
        Move-Safe ("extra\validation\" + $_.Name) ("extra\climate\benchmarks\work\" + $_.Name)
    }
}

Move-Safe "extra\climate.xlsx" "extra\climate\workbooks\climate.xlsx"
Move-Safe "extra\profiling.xlsx" "extra\climate\workbooks\profiling.xlsx"
Move-Safe "extra\climate_0.xlsx" "extra\climate\workbooks\archive\climate_0.xlsx"
Move-Safe "extra\climate_1.xlsx" "extra\climate\workbooks\archive\climate_1.xlsx"
Move-Safe "extra\climate.xlsx.inspect.ndjson" "extra\climate\workbooks\inspection\climate.xlsx.inspect.ndjson"

$sourceDirectories = @(
    "clara-a3-albedo-1991-2020",
    "era5-planette-2001-2020",
    "imerg-v07b-annual",
    "IMERGmonthlyClimatology2001to2022GeoTIFF",
    "ipcc-ar6-wgi-regions-v4",
    "nsidc-g02202-v6-2001-2020",
    "worldclim-2.1-10m"
)
foreach ($name in $sourceDirectories) {
    Move-Safe ("extra\reference\" + $name) ("extra\reference\climate\sources\" + $name)
}
Move-Safe "extra\MOD10CM_61-20260830_152128" "extra\reference\climate\sources\MOD10CM_61-20260830_152128"
Move-Safe "extra\wc2.1_5m_prec" "extra\reference\climate\sources\wc2.1_5m_prec"
Move-Safe "extra\wc2.1_5m_tavg" "extra\reference\climate\sources\wc2.1_5m_tavg"

$referenceRoot = Join-Path $extraRoot "reference"
if (Test-Path -LiteralPath $referenceRoot) {
    Get-ChildItem -LiteralPath $referenceRoot -File | ForEach-Object {
        Move-Safe ("extra\reference\" + $_.Name) ("extra\reference\climate\processed\" + $_.Name)
    }
}

Move-Safe "extra\img\earth\in" "extra\reference\earth\base-maps"
Move-Safe "extra\img\earth\out" "extra\archive\legacy-earth-output"
Move-Safe "extra\testing" "extra\archive\procedural-generation-tests"
Move-Safe "extra\ref\forhinhexes" "extra\reference\projects\forhinhexes"

$literature = @{
    "wind" = @(
        "41558_2020_848_Fig1_HTML.png",
        "GlobalWinds_BN.jpg"
    )
    "precipitation" = @(
        "CFS - 1-15March1993_Atmospheric_Precipitable_Water-small.gif",
        "HuffmanPatterson_Figure9.png",
        "rainydays_merra2_1980to2016.jpg",
        "redding-w.gif",
        "sf1.gif"
    )
    "terrain" = @(
        "DEM-0-in-grayscale-a-and-RGB-converted-image-including-relief-shading-of-the-same.webp",
        "france_rudolf_leuzinger.jpg",
        "Hypsography_USGS_1990.png",
        "lihs7oge2gg01.jpg",
        "Picture1.jpg"
    )
    "minerals" = @(
        "EMIT-First-Global-Maps-of-Surface-Minerals-in-Arid-Regions-scaled.jpg",
        "Mineral Systems CONUS AK HI PR.jpg",
        "Mineral_9.png",
        "minerals-15-00980-g001.png"
    )
    "climate" = @(
        "Biomes.png",
        "Sin título-1.png"
    )
    "papers" = @(
        "Short_communication_Analytical_models_for_2D_lands.pdf"
    )
}
foreach ($category in $literature.Keys) {
    foreach ($name in $literature[$category]) {
        Move-Safe ("extra\" + $name) ("extra\reference\literature\" + $category + "\" + $name)
    }
}

Get-ChildItem -LiteralPath $extraRoot -File |
    Where-Object { $_.Name -like "Sin t*tulo-1.png" } |
    ForEach-Object {
        Move-Safe ("extra\" + $_.Name) ("extra\reference\literature\climate\" + $_.Name)
    }

Move-Safe "extra\IPCC-WGI-reference-regions-v4.geojson.txt" `
    "extra\reference\climate\sources\ipcc-ar6-wgi-regions-v4\IPCC-WGI-reference-regions-v4.geojson.txt"
Move-Safe "extra\comandos.txt" "extra\archive\legacy-notes\comandos.txt"
Move-Safe "extra\seed.txt" "extra\archive\legacy-notes\seed.txt"

Write-Host "Climate asset reorganization complete. No files were deleted."
