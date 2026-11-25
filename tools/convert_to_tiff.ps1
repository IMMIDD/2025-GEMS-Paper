<#
.SYNOPSIS
Converts all PNG figures of a folder to TIFF figures with the same names.

.DESCRIPTION
The TIFF files are written with LZW compression, without transparency (on white), and tagged with
300 dpi. The pixels of the figures are not changed. The PNG files are kept, existing TIFF files of
the same name are overwritten. Subfolders are not searched.

The script warns about figures that exceed the limits of PLOS: 2250 x 2625 px (7.5 x 8.75 in at
300 dpi) and 10 MB. Runs on Windows only.

.PARAMETER Folder
Folder with the PNG files, e.g. a result folder.

.PARAMETER Out
Output folder. Default: the same folder.

.EXAMPLE
powershell -ExecutionPolicy Bypass -File tools\convert_to_tiff.ps1 results\2026-11-24_21-55-04_000_demographic_maps
#>
param(
    [Parameter(Mandatory = $true)][string]$Folder,
    [string]$Out)

$ErrorActionPreference = "Stop"

$Dpi = 300

# limits of PLOS
$MaxWidth = 2250
$MaxHeight = 2625
$MaxFileSize = 10MB

Add-Type -AssemblyName System.Drawing

# saves the figure `png` as TIFF file `tif`; returns its size in pixels
function Save-Tiff($png, $tif) {
    $image = New-Object System.Drawing.Bitmap($png)
    try {
        $figure = New-Object System.Drawing.Bitmap($image.Width, $image.Height, [System.Drawing.Imaging.PixelFormat]::Format24bppRgb)
        try {
            $figure.SetResolution($Dpi, $Dpi)
            $graphics = [System.Drawing.Graphics]::FromImage($figure)
            $graphics.Clear([System.Drawing.Color]::White)
            $graphics.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::NearestNeighbor # copy the pixels as they are
            $graphics.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::Half
            $area = New-Object System.Drawing.Rectangle(0, 0, $image.Width, $image.Height)
            $graphics.DrawImage($image, $area, $area, [System.Drawing.GraphicsUnit]::Pixel)
            $graphics.Dispose()

            $codec = [System.Drawing.Imaging.ImageCodecInfo]::GetImageEncoders() | Where-Object { $_.MimeType -eq "image/tiff" }
            $options = New-Object System.Drawing.Imaging.EncoderParameters(1)
            $options.Param[0] = New-Object System.Drawing.Imaging.EncoderParameter(
                [System.Drawing.Imaging.Encoder]::Compression, [long][System.Drawing.Imaging.EncoderValue]::CompressionLZW)
            $figure.Save($tif, $codec, $options)
        }
        finally { $figure.Dispose() }
        return @($image.Width, $image.Height)
    }
    finally { $image.Dispose() }
}


if (-not (Test-Path -LiteralPath $Folder -PathType Container)) { throw "The folder '$Folder' does not exist." }
$Folder = (Resolve-Path -LiteralPath $Folder).ProviderPath
if (-not $Out) { $Out = $Folder }
$Out = (New-Item -ItemType Directory -Force $Out).FullName

$files = @(Get-ChildItem -LiteralPath $Folder -File | Where-Object { $_.Extension -eq ".png" })
if ($files.Count -eq 0) { Write-Warning "There are no PNG files in $Folder." }

$failed = 0
foreach ($file in $files) {
    $tif = Join-Path $Out ($file.BaseName + ".tif")
    try { $size = Save-Tiff $file.FullName $tif }
    catch {
        $failed++
        Write-Warning "$($file.Name) could not be converted: $($_.Exception.Message)"
        continue
    }
    $length = (Get-Item -LiteralPath $tif).Length
    Write-Host ("{0}: {1} x {2} px, {3:0.0} MB" -f (Split-Path -Leaf $tif), $size[0], $size[1], ($length / 1MB))
    if ($size[0] -gt $MaxWidth -or $size[1] -gt $MaxHeight) {
        Write-Warning "$(Split-Path -Leaf $tif) is larger than $MaxWidth x $MaxHeight px"
    }
    if ($length -gt $MaxFileSize) { Write-Warning "$(Split-Path -Leaf $tif) is larger than 10 MB" }
}

Write-Host "$($files.Count - $failed) of $($files.Count) figures saved in $Out"
if ($failed -gt 0) { exit 1 }
