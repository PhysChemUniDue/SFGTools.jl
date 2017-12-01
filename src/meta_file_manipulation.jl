# Maniplate entries in all metadata files to match new format.
# pwd() has to be ./2011-11-11/name/... while executing

timestamp = string(Int64(Dates.value(now())))

# Loop through all files and search for metadata specifically
for (root, dirs, files) in walkdir(".")

    meta_files = files[contains.(files, "data.txt")]
    for file in meta_files

        path = joinpath(root, file)

        # Read File
        open(path, "r+") do f
            global content = readstring(f)
        end

        # Expressions to be replaced
        regex1 = r"Micos_Positions_X_Y_Z (?<x>\d*\.?\d*) (?<y>\d*\.?\d*) (?<z>\d*\.?\d*)"

        if ismatch(regex1, content)

            # Move the file to a backup destination
            backup_dest = joinpath(root, "../meta_backup", timestamp)
            if !isdir(backup_dest)
                mkpath(backup_dest)
            end
            mv(path, joinpath(backup_dest, replace(file, "data", "backup")))

            sub1 = Base.SubstitutionString("micos_position_x\t\\g<x>\r\nmicos_position_y\t\\g<y>\r\nmicos_position_z\t\\g<z>")
            repl = replace(content, regex1, sub1)

            # Write to file
            open(path, "w") do f
                print(f, repl)
            end

            println("Changed $path")
        end
    end
end
