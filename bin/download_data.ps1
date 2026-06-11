# Navigate to project root directory
$PSScriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location "$PSScriptRoot\.."

$DataDir = "data"
if (!(Test-Path $DataDir)) {
    New-Item -ItemType Directory -Path $DataDir | Out-Null
}

$Samples = @(
    "SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513",
    "SRR1039516", "SRR1039517", "SRR1039520", "SRR1039521"
)

# Exact byte sizes from EBI database to verify completion
$ExpectedSizes = @{
    "SRR1039508_1" = 1243786980
    "SRR1039508_2" = 1230997455
    "SRR1039509_1" = 1144872640
    "SRR1039509_2" = 1152074727
    "SRR1039512_1" = 1592766115
    "SRR1039512_2" = 1588864003
    "SRR1039513_1" = 918567614
    "SRR1039513_2" = 921283600
    "SRR1039516_1" = 1478023204
    "SRR1039516_2" = 1481099287
    "SRR1039517_1" = 2034416969
    "SRR1039517_2" = 2029298187
    "SRR1039520_1" = 1162698060
    "SRR1039520_2" = 1164935922
    "SRR1039521_1" = 1289159510
    "SRR1039521_2" = 1290887405
}

function Get-EbiUrl($srr, $read_num) {
    $first6 = $srr.Substring(0, 6)
    $last_char = $srr.Substring($srr.Length - 1, 1)
    $suffix = "00$last_char"
    return "https://ftp.sra.ebi.ac.uk/vol1/fastq/$first6/$suffix/$srr/${srr}_$read_num.fastq.gz"
}

$MaxJobs = 4
$RunningJobs = @()

Write-Host "Starting robust parallel downloads (max concurrent: $MaxJobs)..."

foreach ($srr in $Samples) {
    foreach ($read in 1..2) {
        $file = "data/${srr}_$read.fastq.gz"
        $expected = $ExpectedSizes["${srr}_$read"]

        if (Test-Path $file) {
            $actual = (Get-Item $file).Length
            if ($actual -eq $expected) {
                Write-Host "$file is complete. Skipping."
                continue
            } elseif ($actual -gt $expected) {
                Write-Warning "$file is larger than expected ($actual > $expected). Deleting to redownload."
                Remove-Item $file -Force -ErrorAction SilentlyContinue
            } else {
                Write-Host "$file is incomplete ($actual < $expected). Resuming download..."
            }
        }

        $url = Get-EbiUrl $srr $read

        # Manage concurrency limit
        while ($RunningJobs.Count -ge $MaxJobs) {
            $RunningJobs = $RunningJobs | Where-Object { $_.HasExited -eq $false }
            if ($RunningJobs.Count -ge $MaxJobs) {
                Start-Sleep -Seconds 2
            }
        }

        Write-Host "Starting parallel download for $file..."
        
        # Self-healing loop command block
        $cmdBlock = @"
`$url = '$url'
`$file = '$file'
do {
    Write-Host "Downloading `$file..."
    curl.exe --speed-limit 1000 --speed-time 30 -L -C - -o `$file `$url
    `$ec = `$LASTEXITCODE
    if (`$ec -ne 0) {
        Write-Warning "Download of `$file failed/interrupted. Retrying in 5 seconds..."
        Start-Sleep -Seconds 5
    }
} while (`$ec -ne 0)
Write-Host "Finished downloading `$file successfully!"
"@

        # Run the self-healing loop in a background powershell process
        $proc = Start-Process powershell.exe -ArgumentList "-NoProfile", "-Command", $cmdBlock -PassThru -NoNewWindow
        $RunningJobs += $proc
    }
}

# Wait for all remaining background downloads to finish
$Remaining = $RunningJobs | Where-Object { $_.HasExited -eq $false }
if ($Remaining) {
    Write-Host "Waiting for remaining downloads to finish..."
    $Remaining | Wait-Process
}

Write-Host "All downloads completed successfully!"
