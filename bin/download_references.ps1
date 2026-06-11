l# Navigate to project root directory
$PSScriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location "$PSScriptRoot\.."

$RefsDir = "refs"
if (!(Test-Path $RefsDir)) {
    New-Item -ItemType Directory -Path $RefsDir | Out-Null
}

$GenomeUrl = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
$GtfUrl = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"

$GenomeFileGz = "$RefsDir/GRCh38.primary_assembly.genome.fa.gz"
$GtfFileGz = "$RefsDir/gencode.v46.annotation.gtf.gz"

$GenomeFile = "$RefsDir/GRCh38.primary_assembly.genome.fa"
$GtfFile = "$RefsDir/gencode.v46.annotation.gtf"

function Un-Gzip($GzipFile, $OutputFile) {
    try {
        $inputFS = [System.IO.File]::OpenRead($GzipFile)
        $outputFS = [System.IO.File]::Create($OutputFile)
        $gzipStream = New-Object System.IO.Compression.GzipStream($inputFS, [System.IO.Compression.CompressionMode]::Decompress)
        $gzipStream.CopyTo($outputFS)
        $gzipStream.Close()
        $inputFS.Close()
        $outputFS.Close()
        return $true
    } catch {
        Write-Error "Failed to decompress $GzipFile : $_"
        return $false
    }
}

function Download-With-Resume($url, $file) {
    do {
        Write-Host "Downloading $file..."
        curl.exe --speed-limit 1000 --speed-time 30 -L -C - -o $file $url
        $ec = $LASTEXITCODE
        if ($ec -ne 0) {
            Write-Warning "Download of $file failed or timed out. Retrying in 5 seconds..."
            Start-Sleep -Seconds 5
        }
    } while ($ec -ne 0)
}

Write-Host "=========================================================="
Write-Host "Downloading GRCh38 Primary Assembly Genome & GTF Reference..."
Write-Host "=========================================================="

if (!(Test-Path $GenomeFile) -or ((Get-Item $GenomeFile).Length -lt 100MB)) {
    Download-With-Resume -url $GenomeUrl -file $GenomeFileGz
} else {
    Write-Host "Genome file already exists or is already extracted."
}

if (!(Test-Path $GtfFile) -or ((Get-Item $GtfFile).Length -lt 10MB)) {
    Download-With-Resume -url $GtfUrl -file $GtfFileGz
} else {
    Write-Host "GTF Annotation already exists or is already extracted."
}

# Extract using native PowerShell GZip decompression
if (Test-Path $GenomeFileGz) {
    Write-Host "Extracting Genome FASTA..."
    $success = Un-Gzip -GzipFile $GenomeFileGz -OutputFile $GenomeFile
    if ($success -and (Test-Path $GenomeFile) -and ((Get-Item $GenomeFile).Length -gt 0)) {
        Remove-Item $GenomeFileGz
        Write-Host "Genome FASTA extracted successfully."
    } else {
        Write-Error "Extraction of Genome FASTA failed."
    }
}

if (Test-Path $GtfFileGz) {
    Write-Host "Extracting GTF Annotation..."
    $success = Un-Gzip -GzipFile $GtfFileGz -OutputFile $GtfFile
    if ($success -and (Test-Path $GtfFile) -and ((Get-Item $GtfFile).Length -gt 0)) {
        Remove-Item $GtfFileGz
        Write-Host "GTF Annotation extracted successfully."
    } else {
        Write-Error "Extraction of GTF Annotation failed."
    }
}

Write-Host "References task execution finished!"
