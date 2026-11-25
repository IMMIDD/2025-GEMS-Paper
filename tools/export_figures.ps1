<#
.SYNOPSIS
Exports the drawn figures of the paper from PowerPoint in their print sizes (300 dpi).

.DESCRIPTION
Needs Windows with PowerPoint. Slides are 7.5 x 8.75 in (19.05 x 22.23 cm), a full page of PLOS.
For every slide that is not empty, the script

- exports the slide through PowerPoint (2250 px wide)
- keeps its left part up to the width in the name of the slide layout, e.g. "Full width (19.05 cm)"
  (2250 px) or "Text width (13.21 cm)" (1560 px)
- crops the empty space above and below (2 pt of white space remain)
- saves a PNG without transparency, tagged with 300 dpi (and a TIFF with -Tiff)

The notes of a slide control its export:

- first line: name of the file, e.g. "Fig_02_layered_architecture" (default: the slide number, e.g.
  "Slide03"), or "skip" for slides that are not exported
- second line "panel": the slide is a panel that Julia combines with a plot. With -UpdatePanels, it is
  also saved in the folder "panels" next to the deck, where Julia reads it

The script warns about objects outside the figure, fonts other than Times New Roman, Arial, and Symbol,
and text smaller than 8 or larger than 12 pt. Text in pictures and charts is not checked. Save the deck
first: unsaved changes are not exported.

.PARAMETER Deck
PowerPoint file (.pptx). Default: figures/drawn_figures.pptx.

.PARAMETER Out
Output folder. Default: a new folder results/<timestamp>_drawn_figures.

.PARAMETER UpdatePanels
Also replaces the panels that are part of the repository (folder "panels" next to the deck).

.PARAMETER Tiff
Also writes TIFF files (LZW compression), as required by PLOS for the submission.

.EXAMPLE
powershell -ExecutionPolicy Bypass -File tools\export_figures.ps1

.EXAMPLE
powershell -ExecutionPolicy Bypass -File tools\export_figures.ps1 -UpdatePanels -Tiff
#>
param(
    [string]$Deck,
    [string]$Out,
    [switch]$UpdatePanels,
    [switch]$Tiff)

$ErrorActionPreference = "Stop"

# print sizes (slides in points, figures in pixels)
$Dpi = 300
$SlideWidth = 540
$SlideHeight = 630
$Margin = 9 # px (2 pt)

# rules of PLOS
$Fonts = @("Times New Roman", "Arial", "Symbol")
$MinFontSize = 8
$MaxFontSize = 12
$MaxFileSize = 10MB

Add-Type -ReferencedAssemblies System.Drawing -TypeDefinition @"
using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

public static class FigureImage
{
    // first and last row with a pixel that is not white, in the left columns up to 'width'
    // (null if there is none)
    public static int[] ContentRows(Bitmap image, int width)
    {
        BitmapData data = image.LockBits(new Rectangle(0, 0, image.Width, image.Height),
            ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
        int first = -1, last = -1;
        byte[] row = new byte[3 * width];
        for (int y = 0; y < image.Height; y++)
        {
            Marshal.Copy(IntPtr.Add(data.Scan0, y * data.Stride), row, 0, row.Length);
            for (int i = 0; i < row.Length; i++)
            {
                if (row[i] < 250)
                {
                    if (first < 0) first = y;
                    last = y;
                    break;
                }
            }
        }
        image.UnlockBits(data);
        return first < 0 ? null : new int[] { first, last };
    }
}
"@

# all shapes of a slide, with the members of groups instead of the groups
function Get-FigureShapes($shapes) {
    foreach ($shape in $shapes) {
        if ($shape.Visible -eq 0) { continue }
        if ($shape.Type -eq 6) { foreach ($item in $shape.GroupItems) { if ($item.Visible -ne 0) { $item } } }
        else { $shape }
    }
}

# text frames of a shape (the cells of a table)
function Get-TextFrames($shape) {
    if ($shape.HasTable) {
        foreach ($row in $shape.Table.Rows) { foreach ($cell in $row.Cells) { $cell.Shape.TextFrame2 } }
    }
    elseif ($shape.HasTextFrame) { $shape.TextFrame2 }
}

# problems of a slide: objects outside the figure (`width` in points), fonts, font sizes
function Get-FigureProblems($slide, $width) {
    $sizes = @()
    foreach ($shape in (Get-FigureShapes $slide.Shapes)) {
        # bounding box of the (rotated) shape
        $angle = $shape.Rotation * [Math]::PI / 180
        $w = [Math]::Abs($shape.Width * [Math]::Cos($angle)) + [Math]::Abs($shape.Height * [Math]::Sin($angle))
        $h = [Math]::Abs($shape.Width * [Math]::Sin($angle)) + [Math]::Abs($shape.Height * [Math]::Cos($angle))
        $x = $shape.Left + $shape.Width / 2
        $y = $shape.Top + $shape.Height / 2
        if ($x - $w / 2 -lt -0.5 -or $x + $w / 2 -gt $width + 0.5 -or $y - $h / 2 -lt -0.5 -or $y + $h / 2 -gt $SlideHeight + 0.5) {
            "'$($shape.Name)' reaches outside the figure and is cut off"
        }

        foreach ($frame in (Get-TextFrames $shape)) {
            if (-not $frame.HasText) { continue }
            if ($frame.AutoSize -eq 2) { "'$($shape.Name)' shrinks its text automatically (AutoFit)" }
            foreach ($run in $frame.TextRange.Runs()) {
                $text = ($run.Text -replace "\s+", " ").Trim()
                if ($text -eq "") { continue }
                if ($text.Length -gt 30) { $text = $text.Substring(0, 30) + "..." }
                $size = [Math]::Round($run.Font.Size, 1)
                if ($size -lt $MinFontSize -or $size -gt $MaxFontSize) { $sizes += [pscustomobject]@{ Size = $size; Text = $text } }
                if ($Fonts -notcontains $run.Font.Name) { "font $($run.Font.Name): '$text'" }
            }
        }
    }
    # one line for all font sizes
    if ($sizes.Count -gt 0) {
        $list = $sizes | Group-Object Size | Sort-Object { $_.Group[0].Size } |
            ForEach-Object { "$($_.Group[0].Size) pt ($($_.Count)x, e.g. '$($_.Group[0].Text)')" }
        "text outside $MinFontSize-$MaxFontSize pt: " + ($list -join ", ")
    }
}

# lines of the notes of a slide
function Get-NoteLines($slide) {
    foreach ($shape in $slide.NotesPage.Shapes) {
        if ($shape.Type -eq 14 -and $shape.PlaceholderFormat.Type -eq 2) {
            return @($shape.TextFrame.TextRange.Text -split "[\r\n\v]" | ForEach-Object { $_.Trim() })
        }
    }
}

# width of the figure of a slide (points), from the name of its layout
function Get-FigureWidth($slide) {
    if ($slide.CustomLayout.Name -notmatch '\((\d+\.?\d*) cm\)') {
        throw "Slide $($slide.SlideIndex): the layout '$($slide.CustomLayout.Name)' does not name its width, e.g. '(13.21 cm)'."
    }
    return [double]::Parse($Matches[1], [Globalization.CultureInfo]::InvariantCulture) / 2.54 * 72
}

# crops the exported slide to the figure (`width` in pixels) and saves it; returns its size in pixels
function Save-Figure($png, $width, $path) {
    $image = New-Object System.Drawing.Bitmap($png)
    try {
        $rows = [FigureImage]::ContentRows($image, $width)
        if ($null -eq $rows) { return $null }
        $top = [Math]::Max($rows[0] - $Margin, 0)
        $height = [Math]::Min($rows[1] + $Margin, $image.Height - 1) - $top + 1
        $figure = New-Object System.Drawing.Bitmap($width, $height, [System.Drawing.Imaging.PixelFormat]::Format24bppRgb)
        $figure.SetResolution($Dpi, $Dpi)
        $graphics = [System.Drawing.Graphics]::FromImage($figure)
        $graphics.Clear([System.Drawing.Color]::White)
        $graphics.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::NearestNeighbor # copy the pixels as they are
        $graphics.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::Half
        $graphics.DrawImage($image, (New-Object System.Drawing.Rectangle(0, 0, $width, $height)),
            (New-Object System.Drawing.Rectangle(0, $top, $width, $height)), [System.Drawing.GraphicsUnit]::Pixel)
        $graphics.Dispose()

        $figure.Save("$path.png", [System.Drawing.Imaging.ImageFormat]::Png)
        if ($Tiff) {
            $codec = [System.Drawing.Imaging.ImageCodecInfo]::GetImageEncoders() | Where-Object { $_.MimeType -eq "image/tiff" }
            $options = New-Object System.Drawing.Imaging.EncoderParameters(1)
            $options.Param[0] = New-Object System.Drawing.Imaging.EncoderParameter(
                [System.Drawing.Imaging.Encoder]::Compression, [long][System.Drawing.Imaging.EncoderValue]::CompressionLZW)
            $figure.Save("$path.tif", $codec, $options)
        }
        $figure.Dispose()
        return @($width, $height)
    }
    finally { $image.Dispose() }
}


$root = Split-Path $PSScriptRoot
if (-not $Deck) { $Deck = Join-Path $root "figures\drawn_figures.pptx" }
$Deck = (Resolve-Path $Deck).Path
if (-not $Out) { $Out = Join-Path $root ("results\{0}_drawn_figures" -f (Get-Date -Format "yyyy-MM-dd_HH-mm-ss_fff")) }
$Out = (New-Item -ItemType Directory -Force $Out).FullName

# work on a copy: the deck may be open in PowerPoint
$temp = Join-Path ([IO.Path]::GetTempPath()) ("figures_" + [Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory $temp | Out-Null
$copy = Join-Path $temp "deck.pptx"
Copy-Item $Deck $copy

$powerpoint = New-Object -ComObject PowerPoint.Application
$presentation = $null
try {
    try { $presentation = $powerpoint.Presentations.Open($copy, -1, 0, 0) } # read-only, without a window
    catch { throw "PowerPoint cannot open the deck ($($_.Exception.Message)). Save it, close PowerPoint, and try again." }
    $page = $presentation.PageSetup
    if ([Math]::Abs($page.SlideWidth - $SlideWidth) -gt 0.1 -or [Math]::Abs($page.SlideHeight - $SlideHeight) -gt 0.1) {
        throw "The slides must be 19.05 x 22.23 cm (7.5 x 8.75 in)."
    }

    $names = @{}
    foreach ($slide in $presentation.Slides) {
        $notes = @(Get-NoteLines $slide)
        if ($notes[0] -eq "skip") {
            Write-Host "Slide $($slide.SlideIndex): skipped (notes)"
            continue
        }
        $name = if ($notes[0] -match '^Fig_\w+$') { $notes[0] } else { "Slide{0:00}" -f $slide.SlideIndex }
        $isPanel = $notes.Count -gt 1 -and $notes[1] -eq "panel"
        if ($names.ContainsKey($name)) { throw "Slides $($names[$name]) and $($slide.SlideIndex) have the same figure name ($name)." }
        $names[$name] = $slide.SlideIndex

        $width = Get-FigureWidth $slide
        $png = Join-Path $temp "$name.png"
        $slide.Export($png, "PNG", [int]($SlideWidth / 72 * $Dpi), [int]($SlideHeight / 72 * $Dpi))
        $size = Save-Figure $png ([int][Math]::Round($width / 72 * $Dpi)) (Join-Path $Out $name)
        if ($null -eq $size) {
            Write-Host "Slide $($slide.SlideIndex): skipped (empty)"
            continue
        }
        Write-Host "Slide $($slide.SlideIndex): $name, $($size[0]) x $($size[1]) px"
        if ($isPanel -and $UpdatePanels) {
            $panels = (New-Item -ItemType Directory -Force (Join-Path (Split-Path $Deck) "panels")).FullName
            Copy-Item (Join-Path $Out "$name.png") $panels -Force
            Write-Host "    updated panels\$name.png"
        }
        foreach ($problem in (Get-FigureProblems $slide $width)) { Write-Warning "$name - $problem" }
        foreach ($file in (Get-ChildItem (Join-Path $Out "$name.*") -Include "$name.png", "$name.tif")) {
            if ($file.Length -gt $MaxFileSize) { Write-Warning "$name - $($file.Name) is larger than 10 MB" }
        }
    }
    Write-Host "Saved in $Out"
}
finally {
    if ($null -ne $presentation) { $presentation.Close() }
    if ($powerpoint.Presentations.Count -eq 0) { $powerpoint.Quit() } # keep PowerPoint open if it is in use
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($powerpoint)
    Remove-Item -Recurse -Force $temp
}
