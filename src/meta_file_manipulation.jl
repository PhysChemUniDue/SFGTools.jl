# Maniplate entries in all metadata files to match new format.
# pwd() has to be ./2011-11-11/name/... while executing

for (root, dirs, files) in walkdir(".")
    meta_files = files[contains.(files, "data.txt")]
    for file in meta_files
        f = open(file, "r+")
        content = read(f)
        close(f)
    end
end

replace(str, r"asdf (?<x>\d*\.?\d*)", s"fdsa\t\1")
