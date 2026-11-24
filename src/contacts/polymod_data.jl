###
### CONTACT DATA: POLYMOD SOCIAL CONTACT SURVEY (download and loading)
###

export obtain_polymod
export load_polymod


"""
    POLYMOD_DIR

Folder of the POLYMOD raw files (not part of the repository, downloaded by
`obtain_polymod`).
"""
const POLYMOD_DIR = joinpath("data", "polymod")

"""
    POLYMOD_URL

Version of the POLYMOD data used in the paper: Mossong J, Hens N, Jit M, et al.
(2020). POLYMOD social contact data (Version 2) [Data set]. Zenodo.
https://doi.org/10.5281/zenodo.3874557. License: CC BY-NC 4.0. Survey described
in Mossong J, et al. (2008). Social contacts and mixing patterns relevant to the
spread of infectious diseases. PLoS Med 5(3):e74.
https://doi.org/10.1371/journal.pmed.0050074.
"""
const POLYMOD_URL = "https://zenodo.org/records/3874557/files/"

"""
    POLYMOD_FILES

Files of `POLYMOD_URL` read by `load_polymod` (and the data dictionary, for
reference) with their SHA-256 checksums.
"""
const POLYMOD_FILES = [
    "2008_Mossong_POLYMOD_participant_common.csv" =>
        "cccca443c96bd0f2d5768a5b4ab489d85265a4cb9e9cbc51f32b38d45dc6ae1a",
    "2008_Mossong_POLYMOD_participant_extra.csv" =>
        "9a0973e84183accd9c99b0f5566018a1cf9b9b499b9874e4fd1ea6ec220ddf6f",
    "2008_Mossong_POLYMOD_hh_common.csv" =>
        "d3799cb076826a95e00a508169df59a9c448970536e87c3d4af1d818cdf56759",
    "2008_Mossong_POLYMOD_sday.csv" =>
        "76cc79379568a42d799a35e0ac0eac07174b85308d40513a7020e7496b15277a",
    "2008_Mossong_POLYMOD_contact_common.csv" =>
        "79a3f0db3fd5dca1c698a26260fb263e061711e82d99a3d26dd999191b0445a6",
    "2008_Mossong_POLYMOD_dictionary.xls" =>
        "c7a0358c6b8075009416d76616333214d00afa25fed6f01e5ccff16d1743ff74"]

"""
    POLYMOD_LOCATIONS

Location flags of a contact in the contact file (several can be set).
"""
const POLYMOD_LOCATIONS = [:cnt_home, :cnt_work, :cnt_school,
    :cnt_transport, :cnt_leisure, :cnt_otherplace]


"""
    obtain_polymod(dir = POLYMOD_DIR)

Downloads the files of `POLYMOD_FILES` from `POLYMOD_URL` into `dir` if they do
not exist, and checks that they are the version used in the paper. Returns `dir`.
"""
function obtain_polymod(dir = POLYMOD_DIR)
    mkpath(dir)
    for (name, sha) in POLYMOD_FILES
        file = joinpath(dir, name)
        if !isfile(file)
            printinfo("Downloading $name from $POLYMOD_URL")
            # download to a temporary file first, so an interrupted download leaves no file
            tmp = file * ".part"
            Downloads.download(POLYMOD_URL * name, tmp)
            mv(tmp, file, force = true)
        end
        checksum = bytes2hex(open(SHA.sha256, file))
        if checksum != sha
            throw(ErrorException("$file is not the version of the POLYMOD data used in " *
                "the paper ($POLYMOD_URL). Delete it to download it again."))
        end
    end
    return dir
end


"""
    load_polymod(; dir = POLYMOD_DIR)

Loads the POLYMOD survey (downloaded if needed, see `obtain_polymod`) and returns
a named tuple:

- `participants`: one row per participant, from the participant, household and
  survey day files. Columns used later: `part_id`, `part_age`, `country`
  (`"BE"`, `"DE"`, `"FI"`, `"GB"`, `"IT"`, `"LU"`, `"NL"`, `"PL"`; the dictionary
  says `"UK"`, the data uses `"GB"`), `type`
  (1: adult, 2: teenager or older child, 3: parent on behalf of a child),
  `part_occupation` (1: working, 2: retired, 3: at home, 4: unemployed,
  5: in full-time education, 6: other; per the dictionary, of the person filling
  in the diary, which is the parent for `type` 3), `child_care_detail` (1: day care at home,
  2: day care outside home, 3: pre-school, 4: school), `dayofweek`
  (0: Sunday to 6: Saturday).
- `contacts`: one row per person contacted on the survey day. Columns used later:
  `part_id`, `cnt_age_exact`, `cnt_age_est_min`, `cnt_age_est_max`, and the
  location flags `POLYMOD_LOCATIONS` (`Bool`).

Codes as in the data dictionary (`2008_Mossong_POLYMOD_dictionary.xls`). Raw
values, no derived columns. The loaded data reproduces Table 1 of Mossong et al.
(2008) exactly (number of participants, mean and standard deviation of contacts by
age, sex, household size, day of the week and country; checked 22.11.2026).
"""
function load_polymod(; dir = POLYMOD_DIR)
    obtain_polymod(dir)
    read_part(name) = CSV.read(joinpath(dir, "2008_Mossong_POLYMOD_$name.csv"), DataFrame,
        missingstring = ["NA", ""])

    common = read_part("participant_common")
    extra = read_part("participant_extra")
    hh = read_part("hh_common")
    sday = read_part("sday")
    contacts = read_part("contact_common")

    # one survey day per participant, otherwise the joins would duplicate participants
    for (name, df) in ["participant_common" => common, "participant_extra" => extra,
            "sday" => sday]
        allunique(df.part_id) ||
            throw(ErrorException("Participant ids in $name are not unique."))
    end

    participants = leftjoin(common, select(hh, :hh_id, :country, :hh_size), on = :hh_id)
    participants = leftjoin(participants, extra, on = :part_id)
    participants = leftjoin(participants, select(sday, :part_id, :dayofweek), on = :part_id)
    sort!(participants, :part_id)

    # every contact must belong to a participant
    issubset(Set(contacts.part_id), Set(participants.part_id)) ||
        throw(ErrorException("Contacts refer to unknown participants."))
    for col in POLYMOD_LOCATIONS
        eltype(contacts[!, col]) <: Union{Missing, Bool} ||
            throw(ErrorException("Location flag $col is not read as Bool."))
    end

    return (participants = participants, contacts = contacts)
end
