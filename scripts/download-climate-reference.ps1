param(
    [string]$Destination = "extra/reference/climate/sources/worldclim-2.1-10m",
    [ValidateSet("core", "wind", "imerg", "all")]
    [string]$VariableSet = "core"
)

$ErrorActionPreference = "Stop"

$variables = if ($VariableSet -eq "all") {
    @("tavg", "prec", "srad", "wind", "vapr")
} elseif ($VariableSet -eq "wind") {
    @("wind")
} elseif ($VariableSet -eq "imerg") {
    @()
} else {
    @("tavg", "prec")
}

$destinationPath = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\$Destination"))
$archivePath = Join-Path $destinationPath "archives"
$gridPath = Join-Path $destinationPath "grids"

[System.IO.Directory]::CreateDirectory($archivePath) | Out-Null
[System.IO.Directory]::CreateDirectory($gridPath) | Out-Null

foreach ($variable in $variables) {
    $fileName = "wc2.1_10m_$variable.zip"
    $url = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/$fileName"
    $zipPath = Join-Path $archivePath $fileName

    if (-not (Test-Path -LiteralPath $zipPath)) {
        Write-Host "Downloading $fileName"
        Invoke-WebRequest -Uri $url -OutFile $zipPath
    }

    Write-Host "Extracting $fileName"
    Expand-Archive -LiteralPath $zipPath -DestinationPath $gridPath -Force
}

$hashes = foreach ($archive in Get-ChildItem -LiteralPath $archivePath -Filter "*.zip") {
    $hash = Get-FileHash -LiteralPath $archive.FullName -Algorithm SHA256
    [PSCustomObject]@{
        file = $archive.Name
        bytes = $archive.Length
        sha256 = $hash.Hash.ToLowerInvariant()
    }
}

$receipt = [ordered]@{
    dataset = "worldclim-2.1-1970-2000-10m"
    downloaded_utc = [DateTime]::UtcNow.ToString("yyyy-MM-ddTHH:mm:ssZ")
    source = "https://www.worldclim.org/data/worldclim21.html"
    archives = @($hashes)
}

$receipt | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $destinationPath "download-receipt.json") -Encoding utf8
Write-Host "WorldClim reference data is ready at $destinationPath"

if ($VariableSet -eq "imerg" -or $VariableSet -eq "all") {
    $imergPath = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\extra\reference\climate\sources\imerg-v07b-annual"))
    $imergGridPath = Join-Path $imergPath "grids"
    $imergArchive = Join-Path $imergPath "IMERG-Final.CLIM.200006-202305.V07B.tif.zip"
    $imergUrl = "https://gpm.nasa.gov/sites/default/files/data/climatologies/2024/IMERG-Final.CLIM.200006-202305.V07B.tif.zip"

    [System.IO.Directory]::CreateDirectory($imergGridPath) | Out-Null

    if (-not (Test-Path -LiteralPath $imergArchive)) {
        Write-Host "Downloading IMERG annual precipitation climatology"
        Invoke-WebRequest -Uri $imergUrl -OutFile $imergArchive
    }

    Expand-Archive -LiteralPath $imergArchive -DestinationPath $imergGridPath -Force
    $imergHash = Get-FileHash -LiteralPath $imergArchive -Algorithm SHA256
    $imergReceipt = [ordered]@{
        dataset = "imerg-final-v07b-grand-average-2000-2023"
        downloaded_utc = [DateTime]::UtcNow.ToString("yyyy-MM-ddTHH:mm:ssZ")
        source = "https://gpm.nasa.gov/data/imerg/precipitation-climatology"
        archive = [ordered]@{
            file = [System.IO.Path]::GetFileName($imergArchive)
            bytes = (Get-Item -LiteralPath $imergArchive).Length
            sha256 = $imergHash.Hash.ToLowerInvariant()
        }
    }

    $imergReceipt | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $imergPath "download-receipt.json") -Encoding utf8
    Write-Host "IMERG reference data is ready at $imergPath"
}
