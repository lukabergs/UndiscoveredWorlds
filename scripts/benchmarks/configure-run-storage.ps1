[CmdletBinding()]
param(
    [string]$StorageRoot = "D:/dev/undiscovered-worlds",
    [string[]]$RelativePaths = @(),
    [switch]$EnsureArchive,
    [switch]$Configure
)

$ErrorActionPreference = "Stop"
$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot "../..")).TrimEnd('\', '/')
$configPath = Join-Path $repoRoot "runs/storage.json"
if (Test-Path -LiteralPath $configPath) {
    $config = Get-Content -LiteralPath $configPath -Raw | ConvertFrom-Json
    $StorageRoot = $config.storage_root
}
$storageBase = [IO.Path]::GetFullPath($StorageRoot).TrimEnd('\', '/')
if ($storageBase -eq $repoRoot -or $storageBase.StartsWith($repoRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
    throw "Storage must be outside the code workspace."
}
if (-not (Test-Path -LiteralPath ([IO.Path]::GetPathRoot($storageBase)))) { throw "Storage drive is unavailable: $storageBase" }
$receiptRoot = Join-Path $repoRoot "out/storage-migration"
New-Item -ItemType Directory -Force -Path $receiptRoot | Out-Null

function Assert-Within([string]$Path, [string]$Base) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not $full.StartsWith($Base + '\', [StringComparison]::OrdinalIgnoreCase)) { throw "Path escapes intended root: $full" }
    return $full
}

foreach ($relative in $RelativePaths) {
    if ([IO.Path]::IsPathRooted($relative) -or $relative -match '(^|[\\/])\.\.([\\/]|$)') { throw "Expected a contained relative path." }
    if ($relative -notmatch '^(runs[\\/](diagnostics|maps|fields|logs|checkpoints|reports|work)|out[\\/](diagnostics|scratch))([\\/]|$)') {
        throw "Only explicitly scoped generated-data paths may be relocated: $relative"
    }
    $source = Assert-Within (Join-Path $repoRoot $relative) $repoRoot
    $destination = Assert-Within (Join-Path $storageBase $relative) $storageBase
    if ($EnsureArchive -and $relative -notmatch '^runs[\\/]reports[\\/][^\\/]+$') { throw "Archive must be a direct report subdirectory." }
    $item = Get-Item -LiteralPath $source -ErrorAction SilentlyContinue
    if ($item -and ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
        if ($item.LinkType -ne 'Junction' -or [IO.Path]::GetFullPath($item.Target) -ne $destination) { throw "Unexpected existing link: $source" }
        Write-Output "Already routed: $relative"
        continue
    }
    if ($item -and -not $item.PSIsContainer) { throw "Expected directory: $source" }
    $files = if ($item) { @(Get-ChildItem -LiteralPath $source -File -Recurse -Force) } else { @() }
    if ($item -and @(Get-ChildItem -LiteralPath $source -Recurse -Force -Attributes ReparsePoint).Count) { throw "Nested links require separate review: $source" }
    $bytes = ($files | Measure-Object -Property Length -Sum).Sum
    if ($EnsureArchive -and $bytes -gt 32MB) { throw "Existing archive exceeds automatic setup limit; migrate it between runs: $source" }
    $receiptName = ($relative -replace '[\\/]', '__') + '.json'
    $receiptPath = Join-Path $receiptRoot $receiptName
    if ((Test-Path -LiteralPath $destination) -and -not (Test-Path -LiteralPath $receiptPath)) { throw "Destination already exists without this migration's receipt: $destination" }
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    $receipt = [ordered]@{source=$source; destination=$destination; bytes=$bytes; files=$files.Count; status='copying'; sha256=@{}}
    $receipt | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $receiptPath -Encoding utf8
    Write-Output "Copying $relative ($([math]::Round($bytes / 1GB, 3)) GiB, $($files.Count) files)"
    if ($item) {
        foreach ($directory in Get-ChildItem -LiteralPath $source -Directory -Recurse -Force) {
            $suffix = [IO.Path]::GetRelativePath($source, $directory.FullName)
            New-Item -ItemType Directory -Force -Path (Assert-Within (Join-Path $destination $suffix) $storageBase) | Out-Null
        }
        foreach ($file in $files) {
            $suffix = [IO.Path]::GetRelativePath($source, $file.FullName)
            $target = Assert-Within (Join-Path $destination $suffix) $storageBase
            Copy-Item -LiteralPath $file.FullName -Destination $target -Force
            $before = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash
            $after = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
            if ($before -ne $after -or (Get-Item -LiteralPath $target).Length -ne $file.Length) { throw "Copy verification failed: $file" }
            $receipt.sha256[$suffix] = $before.ToLowerInvariant()
        }
        # No writers may be active. Detect changed file sets or timestamps before cutover.
        $latest = @(Get-ChildItem -LiteralPath $source -File -Recurse -Force)
        if ($latest.Count -ne $files.Count) { throw "Source changed while copying: $source" }
        foreach ($file in $files) {
            $now = Get-Item -LiteralPath $file.FullName
            if ($now.Length -ne $file.Length -or $now.LastWriteTimeUtc -ne $file.LastWriteTimeUtc) { throw "Active writer detected: $file" }
        }
    }
    $receipt.status = 'verified'
    $receipt | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $receiptPath -Encoding utf8
    $backup = Assert-Within ($source + '.storage-migration-backup') $repoRoot
    if (Test-Path -LiteralPath $backup) { throw "Unexpected backup exists: $backup" }
    if ($item) { Move-Item -LiteralPath $source -Destination $backup }
    try {
        New-Item -ItemType Junction -Path $source -Target $destination | Out-Null
        $link = Get-Item -LiteralPath $source
        if ($link.LinkType -ne 'Junction' -or [IO.Path]::GetFullPath($link.Target) -ne $destination) { throw "Junction verification failed." }
    } catch {
        if (-not (Test-Path -LiteralPath $source) -and $item) { Move-Item -LiteralPath $backup -Destination $source }
        throw
    }
    # Both targets are checked absolute paths within the authorized workspace/storage roots.
    if ($item) { Remove-Item -LiteralPath (Assert-Within $backup $repoRoot) -Recurse -Force }
    $receipt.status = 'complete'; $receipt.completed_utc = [DateTime]::UtcNow.ToString('o')
    $receipt | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $receiptPath -Encoding utf8
    Write-Output "Verified and routed: $relative"
}
if ($Configure) {
    [ordered]@{version=1; storage_root=$storageBase; archive_parent='runs/reports';
        keep_on_primary=@('src','scripts','refs','out/build','runs/registry','runs/metrics','runs/work/climate/seed_20260906');
        method='Directory junctions preserve logical paths; each new experiment archive is routed by configure-run-storage.ps1 -EnsureArchive.'} |
        ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $configPath -Encoding utf8
    Write-Output "Default storage configured: $storageBase"
}
