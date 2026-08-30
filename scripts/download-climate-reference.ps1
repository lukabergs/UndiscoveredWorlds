param(
    [string]$Destination = "extra/reference/worldclim-2.1-10m",
    [ValidateSet("core", "all")]
    [string]$VariableSet = "core"
)

$ErrorActionPreference = "Stop"

$variables = if ($VariableSet -eq "all") {
    @("tavg", "prec", "srad", "wind", "vapr")
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
