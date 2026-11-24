# Extracts the tables of the BA commuter statistics (30.06.2019) in this folder
# into CSV files, because Julia cannot read the binary Excel format (.xlsb).
# The rows and values are copied unchanged ("*" = anonymized, "-" = none):
#   gemband_dlk_0.xlsb, sheet "Gemeindedaten"    -> ba_commuters_by_area_2019.csv
#   krpend_XX_0.xlsb,   sheet "Auspendler Kreise" -> ba_commuter_flows_by_county_2019.csv
# Requires Microsoft Excel (Windows). Opens the files read-only.
# Usage (from the repository root): powershell -File inputs/commuters/extract_ba_commuters.ps1

$ErrorActionPreference = "Stop"
$dir = $PSScriptRoot
$utf8 = New-Object System.Text.UTF8Encoding($false)

$source = "# Source: Statistik der Bundesagentur fuer Arbeit, Tabellen, {0}, Nuernberg, Stichtag 30. Juni 2019 ({1}); https://statistik.arbeitsagentur.de/"
$license = "# License: (c) Statistik der Bundesagentur fuer Arbeit. Storing, copying and redistributing (also in excerpts) with source note; contents must not be altered."

function Csv-Field($v) {
    $s = if ($null -eq $v) { "" } else { ([string]$v).Trim() -replace "\s*\r?\n\s*", " " }
    if ($s -match '[",]') { '"' + $s.Replace('"', '""') + '"' } else { $s }
}

# all values of a sheet as a 2D array (1-based), read at once
function Sheet-Values($wb, $name) { return ,$wb.Worksheets.Item($name).UsedRange.Value2 }

$xl = New-Object -ComObject Excel.Application
$xl.Visible = $false; $xl.DisplayAlerts = $false; $xl.AutomationSecurity = 3
try {
    # --- employees by place of residence and place of work, commuters: Germany, states, districts, municipalities
    $wb = $xl.Workbooks.Open((Join-Path $dir "gemband_dlk_0.xlsb"), 0, $true)
    $v = Sheet-Values $wb "Gemeindedaten"
    $lines = @(
        "# Employees subject to social insurance (SvB) on 30.06.2019 by place of residence and place of work, and commuters, per area (Germany, states, administrative districts, districts, municipalities).",
        ($source -f "Gemeindedaten der sozialversicherungspflichtig Beschaeftigten nach Wohn- und Arbeitsort", "gemband_dlk_0.xlsb, sheet Gemeindedaten"),
        $license,
        "# Changes: rows of the sheet copied unchanged (empty rows and notes dropped); columns renamed. Codes: 8 digits = municipality, 5 = district (urban districts only as municipality), 3 = administrative district, 2 = state. ""*"" = anonymized (values 1 or 2), empty = none.",
        "ags,name,residence,residence_men,residence_women,residence_german,residence_foreign,residence_under_20,residence_20_to_24,residence_55_plus,workplace,residence_is_workplace,in_commuters,out_commuters,establishments")
    for ($r = 1; $r -le $v.GetLength(0); $r++) {
        $code = ([string]$v[$r, 1]).Trim()
        if ($code -notmatch '^\d+$') { continue }
        $lines += (@(1..15 | ForEach-Object { Csv-Field $v[$r, $_] }) -join ",")
    }
    $wb.Close($false)
    [System.IO.File]::WriteAllLines((Join-Path $dir "ba_commuters_by_area_2019.csv"), [string[]]$lines, $utf8)
    Write-Host "ba_commuters_by_area_2019.csv: $($lines.Count - 5) rows"

    # --- out-commuters per district of residence and district (or aggregate) of work
    $lines = @(
        "# Out-commuters (employees subject to social insurance) on 30.06.2019 per district of residence (home) and district or aggregate of work (destination). Destination codes: 5 digits = district; 3 digits = administrative district (ending in 0: the whole state), as total or as ""Uebrige Kreise"" (districts not listed); 2 digits = state, as total or as ""Uebrige Regierungsbezirke"" (administrative districts not listed); ZZ = ""Uebrige Bundeslaender"" (states not listed) or total in Germany; Z = total. Small flows are only contained in these aggregates.",
        ($source -f "Pendlerverflechtungen der sozialversicherungspflichtig Beschaeftigten nach Kreisen", "krpend_01_0.xlsb to krpend_16_0.xlsb, sheet Auspendler Kreise"),
        $license,
        "# Changes: rows of the 16 sheets copied unchanged (empty rows and notes dropped); the district of residence repeated in each row; columns renamed. ""*"" = anonymized (values 1 or 2), ""-"" = none.",
        "home,home_name,destination,destination_name,total,men,women,german,foreign")
    foreach ($s in 1..16) {
        $f = "krpend_{0:D2}_0.xlsb" -f $s
        $wb = $xl.Workbooks.Open((Join-Path $dir $f), 0, $true)
        $v = Sheet-Values $wb "Auspendler Kreise"
        $origin = $null; $originName = $null; $n = 0
        for ($r = 1; $r -le $v.GetLength(0); $r++) {
            $a = ([string]$v[$r, 1]).Trim()
            if ($a -match '^\d{5}$') { $origin = $a; $originName = Csv-Field $v[$r, 2]; continue }
            $c = ([string]$v[$r, 3]).Trim()
            if ($null -eq $origin -or $c -eq "") { continue }
            $lines += ((@($origin, $originName) + @(3..9 | ForEach-Object { Csv-Field $v[$r, $_] })) -join ",")
            $n++
        }
        $wb.Close($false)
        Write-Host "$f : $n rows"
    }
    [System.IO.File]::WriteAllLines((Join-Path $dir "ba_commuter_flows_by_county_2019.csv"), [string[]]$lines, $utf8)
    Write-Host "ba_commuter_flows_by_county_2019.csv: $($lines.Count - 5) rows"
} finally {
    $xl.Quit()
    [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($xl)
}
