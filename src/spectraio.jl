import FileIO.load
using AndorSIF
using CSV
using Dates
using DelimitedFiles
using JSON: json
import LightXML

"""
    list_spectra(; exact="", inexact="", date::Tuple{Int64,Int64,Int64}, group=false)

Return a `DataFrame` with matching spectra.
First `grab(datafolder)` the spectra to genereate a `.spectralist` file inside your current working
directory. Note that you have to have the `.spectralist` file in your current working directory
for `list_spectra` to work.

Filter the dataframe:
* `exact="name"`: Match the exact name of the spectrum
* `inexact="part"`: Check if the name contains part of the string
* `date=(2011, 01, 31)`: Filter for specific date
* `group=true`: Groups spectra by name

The result with some useful information is stored in a DataFrame. No spectra are
loaded hereby to save time. Load the actual spectra via `load_spectra(id)`
where `id` is `list_spectra().id` See the DataFrames package for more info.
"""
function list_spectra(; exact=""::AbstractString,
                        inexact=""::AbstractString,
                        date=(0,0,0)::Tuple{Int64,Int64,Int64},
                        group=false)

  dir = "./"

  spectrafile = joinpath(dir, ".spectralist")
  if isfile(spectrafile)
    dat = readdlm(spectrafile; comments=false)
  else
    error("The file $spectrafile does not exist.")
  end

  df = CSV.File(spectrafile) |> DataFrame

  # Filter the dataframe
  if exact != ""
      df = df[lowercase.(df.name) .== lowercase(exact), :]
  end

  if inexact != ""
      df = df[occursin.(lowercase(inexact), lowercase.(df.name)), :]
  end

  if !all(iszero.(date))
      df = df[Dates.yearmonthday.(df.date) .== [date], :]
  end

  sort!(df)

  if group
      df = combine(df -> DataFrame(
            N = length(df.id),
            sizes = [unique(df.sizes)],
            dates = [unique(floor.(df.date, Dates.Day(1)))],
            id = [df.id]
        ),
        groupby(df, :name)
      )
  end

  return df

end


"""
    grab(dir="./"; getall=false, singledata="")

Make a data file that contains information where to find spectra connected to
an ID. All spectra are internally referenced by this id. `dir` is the main data
directory (`dir/2011-01-31/Name/...`).

By default this function writes only newly found spectra to the .spectralist file.
If you want to rewrite the whole file pass `getall=true` as a keyword argument
to the function.

Returns the number of added spectra and the number of spectra in total.

Optional: specify a single data folder to grab with singledata="path/to/datafolder"
"""
function grab(dir="./"; getall=false, singledata="")


    #grab single data folder in $singledata?
    if singledata != ""

        idexisting = Int64[]
        idlist = Int64[]
        namelist = String[]
        dirlist = String[]
        datelist = DateTime[]
        sizelist = []
        numlist = Int64[]

        for (root, dirs, files) in walkdir(singledata)
            for dir in dirs
                if dir == "raw"

                    # search for xml files first
                    format = :xml
                    mlist = searchdir(joinpath(root, dir), "data.xml")

                    # if still empty throw an error
                    isempty(mlist) && error("No meta files found in $(root)/$(dir)")

                    mdict = get_metadata(joinpath(root, dir, mlist[1]))

                    # Get ID
                    id = Int64(Dates.value(DateTime(mdict["timestamp"])))
                    if !any(idexisting .== id)
                        push!(idlist, id)
                        push!(namelist, splitdir(splitdir(root)[1])[2])
                        push!(dirlist, joinpath(abspath(root), dir))
                        push!(datelist, DateTime(mdict["timestamp"]))
                        if format == :xml
                            sif_files = searchdir(joinpath(root, dir), "data.sif")
                            filestats = stat(joinpath(root, dir, sif_files[1]))
                            specsize = filestats.size
                            push!(sizelist, specsize)
                        elseif format == :txt
                            push!(sizelist, Int64[N_PIXEL/mdict["x_binning"], N_PIXEL/mdict["y_binning"], length(mlist)])
                        end
                        push!(numlist, parse(Int64, splitdir(splitdir(dirlist[end])[1])[2]))
                    end
                end
            end
        end

        # # create new SPECTRALIST file
        # df = DataFrame(id=idlist, name=namelist, path=dirlist, date=datelist, sizes=sizelist, number=numlist)
        # CSV.write(".spectralist", df; append=false)

    # grab all spectra files
    else

        # Check if the .spectralist file exists.
        # If it exist read its contents.
        if !isfile(".spectralist") || getall
            idexisting = Int64[]
        else
            df = CSV.read(".spectralist", DataFrame)
            idexisting = convert.(Int64, df.id)
        end

        idlist = Int64[]
        namelist = String[]
        dirlist = String[]
        datelist = DateTime[]
        sizelist = []
        numlist = Int64[]

        for (root, dirs, files) in walkdir(dir)
            for dir in dirs
                if dir == "raw"

                    # search for xml files first
                    format = :xml
                    mlist = searchdir(joinpath(root, dir), "data.xml")
                    # if the list is empty search for txt files
                    if isempty(mlist)
                        mlist = searchdir(joinpath(root, dir), "data.txt")
                        format = :txt
                    end

                    # if still empty throw an error
                    isempty(mlist) && error("No meta files found in $(root)/$(dir)")

                    mdict = get_metadata(joinpath(root, dir, mlist[1]))

                    # Get ID
                    id = Int64(Dates.value(DateTime(mdict["timestamp"])))
                    if !any(idexisting .== id)
                        push!(idlist, id)
                        push!(namelist, splitdir(splitdir(root)[1])[2])
                        push!(dirlist, joinpath(abspath(root), dir))
                        push!(datelist, DateTime(mdict["timestamp"]))
                        if format == :xml
                            sif_files = searchdir(joinpath(root, dir), "data.sif")
                            filestats = stat(joinpath(root, dir, sif_files[1]))
                            specsize = filestats.size
                            push!(sizelist, specsize)
                        elseif format == :txt
                            push!(sizelist, Int64[N_PIXEL/mdict["x_binning"], N_PIXEL/mdict["y_binning"], length(mlist)])
                        end
                        push!(numlist, parse(Int64, splitdir(splitdir(dirlist[end])[1])[2]))
                    end
                end
            end
        end
    end

    df = DataFrame(id=idlist, name=namelist, path=dirlist, date=datelist, sizes=sizelist, number=numlist)

    # open(".spectralist", "a") do f
    #     # writedlm(f, [idlist namelist dirlist datelist sizelist numlist])
    #     CSV.write(f, df)
    # end

    if !isfile(".spectralist") || getall || singledata != ""
        CSV.write(".spectralist", df; append=false)
    else
        CSV.write(".spectralist", df; append=true)
    end

    # println("Collected $(length(idlist)) spectra. ($(length(idexisting) + length(idlist)) overall)")

    length(idlist), length(idexisting) + length(idlist)


end

"""
    load_spectra(id::Int64)

Load the spectra specified by `id`.
"""
function load_spectra(id::Int64; format=:sif)
  # Load grabbed data
  dir = getdir(id)

  # Load the spectrum files
  s = read_as_3D(dir; format=format)

  sfspectrum = SFSpectrum[SFSpectrum(id, s)]
end

"""
    load_spectra(id::Array{Int64})

Load the spectra specified by `id`.
"""
function load_spectra(id::AbstractArray{Int64}; format=:sif)
  sfspectra = Array{SFSpectrum,1}(undef, length(id))
  for i in 1:length(id)
    sfspectrum = load_spectra(id[i], format=format)
    sfspectra[i] = sfspectrum[1]
  end
  return sfspectra
end

"""
Get an attribute of a spectrum or an array of spectra. To get available
attribute list these via `get_metadata(spectrum)`.
"""
function get_attribute(s::Array{SFSpectrum}, attr::AbstractString)
    values = []
    for i in s
        val = get_attribute(i, attr)
        push!(values, val)
    end
    return values
end

function get_attribute(id::Int64, attr::AbstractString)

    """
    Flattens an array of dicts into a dict of arrays.
    """
    function flattendictarray(dictarray)
    typeof(dictarray) <: Dict && (return dictarray)
    d = Dict()
        for k in keys(dictarray[1])
            d[k] = [x[k] for x in dictarray]
        end
    d
    end

    dict = get_metadata(id)
    # check if we deal with an old txt metadatafile or a new xml metadata file
    # (the old one should have an x_binning key somewhere)
    haskey(dict, "x_binning") ? (format = :txt) : (format = :xml)
    if format == :txt
        # metafile is in txt format. just return the value for the key
        attr in keys(dict) ? val = dict[attr] : val = nothing
        return val
    elseif format == :xml
        # format is xml lookup the dictionary to get to the right entry
        translate = Dict(
            "name" => splitpath(getdir(id))[end-2],
            "comment" => dict["comment"],
            "grating" => "n/a", #dict["spectra pro 300"]["grating"],
            "spectrometer_wavelength" => dict["spectra pro 300"]["wavelength set"],
            "mirror_position" => "n/a", #dict["spectra pro 300"]["mirror position"],
            "pump_dl_position" => dict["ls-180"]["pos"],
            "ccd_exposure_time" => dict["ixon"]["exposure_time"],
            "ccd_readout_speed" => "n/a",
            "x_binning" => dict["ixon"]["horizontal_bin"],
            "y_binning" => dict["ixon"]["vertical_bin"],
            "ccd_temperature" => dict["ixon"]["temperature"],
            "timestamp" => dict["timestamp"],
            "twin1_wavelength" => dict["twins"]["twin1 wavelength"],
            "twin2_wavelength" => dict["twins"]["twin2 wavelength"],
            "twin1_shutter" => dict["twins"]["twin1 shutter"],
            "twin2_shutter" => dict["twins"]["twin2 shutter"],
            "pump_ir_wavelength" => dict["ekspla laser"]["ekspla wavelength"],
            "micos_position_x" => dict["smc stages"]["xpos"],
            "micos_position_y" => dict["smc stages"]["ypos"],
            "micos_position_z" => dict["smc stages"]["zpos"],
            "probe_ir_power" => dict["ir power meters"]["probe ir power"] |> flattendictarray,
            "pump_ir_power" => dict["ir power meters"]["pump ir power"] |> flattendictarray,
            "vis_power" => dict["pm100usb"]["sample power"],
        )
        return translate[attr]
    end
end

get_attribute(s::SFSpectrum, attr::AbstractString) = get_attribute(s.id, attr)


"""
Check if a spectrum has the attribute `attr`.
"""
function is_attribute(id::Int64, attr::AbstractString)
    dict = get_metadata(id)
    attr in keys(dict)
end

function is_attribute(s::Array{SFSpectrum}, attr::AbstractString)
    values = []
    for i in s
        val = is_attribute(i, attr)
        push!(values, val)
    end
    return values
end

is_attribute(s::SFSpectrum, attr::AbstractString) = is_attribute(s.id, attr)

"""
Read tiff files and put them in a 3D matrix
"""
function read_as_3D(path::AbstractString; format=:sif)
    filelist = String[]
    format == :sif && (filelist = searchdir(path, ".sif"))
    # if there are no sif files check for tiff files instead
    isempty(filelist) && (format = :tiff)
    format == :tiff && (filelist = searchdir(path, ".tiff"))
    if isempty(filelist)
        error("Could not find a file with format sif or tiff in $directory.")
    end

    if format == :sif
        # transpose the image that's loaded from the sif file
        I = load(joinpath(path, filelist[1]))
        C = Array{Float64,3}(undef, size(I,1), size(I,2), length(filelist))
        C[:,:,1] .= I[:,:,1]
    elseif format == :tiff
        I = load(joinpath(path, filelist[1]))
        C = Array{eltype(I),3}(undef, size(I,1), size(I,2), length(filelist))
        C[:,:,1] .= I
    end
    @inbounds for i = 2:length(filelist)
        if format == :sif
            C[:,:,i] .= load(joinpath(path, filelist[i]))[:,:,1]
        elseif format == :tiff
            C[:,:,i] .= load(joinpath(path, filelist[i]))
        end
    end
    if format == :sif
        F = C
    elseif format == :tiff
        F = (C .|> Float64) .* typemax(UInt16) .|> round
    end
    return F
end

"""
Get the directory of a spectrum with a given id. If the ID does not exist return an empty string.
"""
function getdir(id::Int64)
    df = CSV.read(".spectralist", DataFrame)
    idx = findall((in)(id), df.id)
    isempty(idx) && error("Could not find spectrum with id $id.")
    dir = df.path[idx[1]]
end


"""
Reads an xml metadata file
"""
function read_xml(path::String)

    function loopelements(element, dict)
        # println("Looping")
        for e in LightXML.child_elements(element)
            # println("Current Element Type: $(LightXML.name(e))")
            if haskey(lvtypes, LightXML.name(e))
                # Go deeper
                if any(LightXML.name(e) .== ["Cluster", "Array"])
                    # Module Data
                    clustername = LightXML.find_element(e, "Name") |> LightXML.content
                    # Get rid of the Module Data in the Name of the Cluster,
                    # strip whitespace, make lowercase
                    clustername_clean =
                        replace(clustername, "Module Data" => "") |> strip |> lowercase
                    # println("Found Cluster: $clustername")
                    dict[clustername_clean] = Dict()
                    loopelements(e, dict[clustername_clean])
                else
                    n, v = parseelement(e)
                    # println("Found Element $n with value $v")
                    global dict[n] = v
                end
            elseif any(LightXML.name(e) .== ["NumElts", "Name", "Dimsize"])
                # Ignore entry
                # println("Ignoring")
                continue
            else
                @warn "Could not find type for $(LightXML.name(e))"
            end
        end
    end

    function parseelement(element)
        typ = lvtypes[LightXML.name(element)]
        n = LightXML.find_element(element, "Name") |> LightXML.content |> strip |> lowercase
        v_str = LightXML.find_element(element, "Val") |> LightXML.content
        if typ == String
            global v = v_str
        elseif typ == Bool
            if v_str == "1"
                v = true
            elseif v_str == "0"
                v = false
            end
        else
            v = parse(typ, v_str)
        end
        n, v
    end

    lvtypes = Dict(
        "Cluster" => Dict,
        "Array"   => Array,
        "Refnum"  => String,
        "DBL"     => Float32,
        "I32"     => Int32,
        "I64"     => Int64,
        "U8"      => UInt8,
        "String"  => String,
        "Boolean" => Bool,
    )

    dict = Dict()

    xdoc = LightXML.parse_file(path)
    xroot = LightXML.root(xdoc)

    loopelements(xroot, dict)

    LightXML.free(xdoc)

    dict
end


"""
Get all available metadata for a spectrum. Takes a SFSpectrum struct or an id as
an Int64.
"""
get_metadata(s::SFSpectrum) = get_metadata(s.id)

function get_metadata(id::Int64)
    dir = getdir(id)
    mdict = get_metadata(dir)
end

function get_metadata(path::AbstractString)
    path
    path == "" && (return Dict())

    mdict = Dict{String, Any}()
    mfiles = Array{String}[]

    if splitext(path)[2] != ".txt" && splitext(path)[2] != ".xml"
        # The path is a folder
        @assert isdir(path)
        # Search for xml files by default
        mfiles = searchdir(path, "data.xml")
        format = :xml
        # If no xml files were found search for txt files
        if length(mfiles) == 0
            mfiles = searchdir(path, "data.txt")
            format = :txt
        end
        paths = Array{String}(undef, size(mfiles, 1))
        [paths[i] = joinpath(path, mfiles[i]) for i = 1:length(mfiles)]
    else
        # Single file
        @assert isfile(path)
        paths = [path]
        mfiles = [path]
        # Get the file extenstion of the input file
        format_str = splitext(path)[2]
        format = replace(format_str, "." => "") |> lowercase |> Symbol
    end

    if format == :xml
        # get the first metadatafile to extract the keys and initialize arrays
        # with the values of the right type
        dict_first = read_xml(paths[1])
        dictkeys = keys(dict_first) #the values of these keys can be dictionaries
        if length(paths) > 1
            for dictkey in dictkeys, k in keys(dict_first[dictkey])
                # So we loop through these dictionaries if dictkey is indeed a
                # dictonary. Otherwise we just put it in on the first layer of
                # the dictionary
                if typeof(dict_first[dictkey]) <: Dict
                    val = dict_first[dictkey][k]
                    valtype = typeof(val)
                    !haskey(mdict, dictkey) &&
                        (mdict[dictkey] = Dict{String, Array{Any}}())
                    # Put the values
                    mdict[dictkey][k] = valtype[val]
                else
                    val = dict_first[dictkey]
                    valtype = typeof(val)
                    # Put the values
                    mdict[dictkey] = valtype[val]
                    # continue # otherwise we loop through all the keys
                end
            end
            for i = 2:length(paths), dictkey in dictkeys, k in keys(dict_first[dictkey])
                dict = read_xml(paths[i])
                if typeof(dict[dictkey]) <: Dict
                    dictkey, k
                    val = dict[dictkey][k]
                    push!(mdict[dictkey][k], val)
                else
                    val = dict[dictkey]
                    push!(mdict[dictkey], val)
                    # continue # otherwise we loop through all the keys
                end
            end
        elseif length(paths) == 1
            mdict = dict_first
        end

        # next we read the metadata from the .sif file
        # First search for the sif files in the folder
        if length(paths) == 1
            # We deal with a single file which has a .xml extension that we need
            # to get rid of and change it to .sif.
            sif_files = [splitext(paths[1])[1] * ".sif"]
        else
            # We deal with a path and perhaps multiple files
            sif_files = joinpath.(path, searchdir(path, "data.sif"))
        end
        if length(sif_files) == 1
            # load the first metadata
            ixon_meta = load(sif_files[1]).ixon
            mdict["ixon"] = ixon_meta
        elseif length(sif_files) > 1
            # load the first metadata
            ixon_meta = load(sif_files[1]).ixon
            mdict["ixon"] = Dict{String, Any}()
            for key in keys(ixon_meta)
                val = ixon_meta[key]
                typ = typeof(val)
                # make an array of the proper type for value val
                mdict["ixon"][key] = typ[val]
            end
            # and loop through the remaining spectra in the folder
            for i = 2:length(sif_files), key in keys(ixon_meta)
                dict = load(joinpath(path, sif_files[i])).ixon
                push!(mdict["ixon"][key], dict[key])
            end
        else
            @warn "No corresponding sif file found in $path"
        end

    elseif format == :txt
        values = Array{Any,2}
        data = readdlm(paths[1], '\t'; comments=false)
        keyz = Array{String}
        keyz = data[:,1]
        value = data[:,2]
        if length(mfiles) == 1
          for (i, key) in enumerate(keyz)
              mdict[key] = value[i]
          end
        elseif length(mfiles) > 1
          valuemat = Array{Any,2}(undef, length(keyz), length(mfiles))
          valuemat[:,1] = data[:,2]
          for i = 2:length(mfiles)
              data = readdlm(paths[i], '\t'; comments=false)
              valuemat[:,i] = data[:,2]
          end
          for (i, key) in enumerate(keyz)
              mdict[key] = valuemat[i,:]
          end
        end
    end

    # We need the comment only once
    # length(mfiles) > 1 && (mdict["comment"] = mdict["comment"][1])
    !haskey(mdict, "comment") && (mdict["comment"] = "")

    return mdict
  end


function searchdir(directory::AbstractString, key::AbstractString)
    filter!(x->occursin(key, x), readdir(directory))
end


"""
`savejson(path::String, spectra::Array{SFSpectrum{T,1},1})`
Saves the array `spectra` to `path` in JSON format.

This is primarily to read the data set by the Matlab fit program.
"""
function savejson(path::String, spectra::Array{SFSpectrum{T,1},1}) where T
    splitext(path)[end] != ".json" && (path = path * ".json")

    for i = 1:length(spectra)
        length(size(spectra[i])) > 1 && error("Please make the Spectra 1D")
    end

    ids = Int64[]
    wavelengths = Array{Float64}[]
    wavenumbers = Array{Float64}[]
    signals = Array{Float64}[]
    names = String[]

    for i = 1:length(spectra)
        λ = get_ir_wavelength(spectra[i])
        ν = 1 ./ λ .* 1e7

        push!(ids, spectra[i].id)
        push!(wavelengths, λ)
        push!(wavenumbers, ν)
        push!(signals, spectra[i].s)
        push!(names, get_attribute(spectra[i], "name"))
    end

    dict = Dict(
        "ids"         => ids,
        "wavelengths" => wavelengths,
        "wavenumbers" => wavenumbers,
        "signals"     => signals,
        "names"       => names
    )

    jsonstring = json(dict)

    open(path, "w+") do io
        write(io, jsonstring)
    end

end
using HDF5

""" Save the data in a HDF5 File. \n
Condition: \n
\t Use Template provided by Marvin or Nelli.\n
Register a new Sample name or choose a sample name from the SFGDashboard Sample Dropdown, e.g. "CaAra" or "CaAra d4".\n
Example Input:\n
    \tsample = "CaAra" \n
    \tsample_prep = "20200504"\n
    \tdate = "20200508"\n
    \tdirectory = "date" * "/" * "CaAra_20200504_1_DL-Scan_r--pumped_ppp"\n
    save_data(sample,32,"ppp","DL d-",date)\n
    save_data(sample::String, surface_density::Int, polarisation_comb::String, scan::String, date::String)
"""
function save_data(sample::String, surface_density::Int, polarisation_comb::String, scan::String, date::String)

    # sample surface density
    surface_density = "$surface_density mNm⁻¹"
        
    # scan type (delay or wavenumber scan)
    if occursin("dl",lowercase(scan)) == true
        scan_type = "delay_scan"
        pump_resonance= split(scan)[2]*"pumped"
    elseif occursin("wl",lowercase(scan)) == true || occursin("wn",lowercase(scan)) == true
        scan_type = "wavenumber_scan"
    else 
        error("""Scan type cannot be $scan !!
        Example for delay scan with d- pumped scan ="dl d-" 
        Example for wavelength scan scan = "wl".
        """)
    end

    # generate filename
    filename = sample *"_"* "$surface_density"*"mNm" *"_"* "$polarisation_comb" *"_"* scan *"_"* sample_prep *".h5"
    
    # calculate pump wavenumber (Ekspla) for delay scan 
    ekspla_wavelength = get_metadata(raw[1])["ekspla laser"]["ekspla wavelength"][1]
    ekspla_wavenumber = 10^7 / ekspla_wavelength
    
    h5open(filename, "w") do fid
        g0 = g_create(fid, sample)
        g1 = g_create(g0, surface_density)
        g2 = g_create(g1, polarisation_comb)
        g3 = g_create(g2, scan_type)
        
        if scan_type == "delay_scan"
            g4 = g_create(g3, pump_resonance)
            g5 = g_create(g4, date)
            
            g6 = g_create(g5, "Data")
            g6["sig_matrix"] = sig03
            g6["ref_matrix"] = ref03
            g6["pump_wavenumber"] = ekspla_wavenumber
            g6["wavenumber"] = ν
            g6["dltime"] = dltime_sorted
            
            for i in 1:length(mode_name)
                g6["sig_mean_$(mode_name[i])"] = mean(sig03[:,pixel[i]], dims=2)[:,1]
            end
            
            for i in 1:length(mode_name)
                g6["ref_mean_$(mode_name[i])"] = mean(ref03[:,pixel[i]], dims=2)[:,1]
            end
            
            g6["comment"] = get_metadata(raw[1])["comment"][1]
            g6["folder_name"] = directory
            
        else scan_type == "wavenumber_scan"
            g4 = g_create(g3, date)
            
            g5 = g_create(g4, "Data")
            g5["sig_matrix"] = sig03
            g5["ref_matrix"] = ref03
            g5["pump_wavenumber"] = ekspla_WN # Ekspla pump wavenumber for wavenumber scan
            g5["wavenumber"] = ν
            g5["dltime"] = dltime_sorted
            
            for i in 1:length(mode_name)
                g5["sig_mean_$(mode_name[i])"] = mean(sig03[:,pixel[i]], dims=2)[:,1]
            end
            
            for i in 1:length(mode_name)
                g5["ref_mean_$(mode_name[i])"] = mean(ref03[:,pixel[i]], dims=2)[:,1]
            end
            
            g5["comment"] = get_metadata(raw[1])["comment"][1]
            g5["folder_name"] = directory    
        end
    end
end

"""
Save Spectra in a Matlab file.
"""
# function save_mat(filename::AbstractString, s::SFSpectrum)
#     save_mat(filename, SFSpectrum[s])
# end
#
# function save_mat(filename::AbstractString, s::Array{SFSpectrum})
#     # Check if filename has a proper extension
#     if splitext(filename)[end] != ".mat"
#         filename *= ".mat"
#     end
#
#     f = matopen(filename, "w")
#
#     length(s)
#
#     for i = 1:length(s)
#
#         # Put Stuff in dict
#         name = ""
#         try
#             name = get_attribute(s[i], "name")
#         catch
#             name = string(s[i].id)
#         end
#
#         data = Dict()
#         try
#             data = Dict(
#                 "name" => get_attribute(s[i], "name"),
#                 "signal" => s[i].s,
#                 "wavelength" => get_ir_wavelength(s[i]),
#                 "wavenumber" => get_ir_wavenumber(s[i]),
#             )
#         catch
#             data = Dict(
#                 "name" => name,
#                 "signal" => s[i].s
#             )
#         end
#         data
#
#         # Write to file
#         write(f, "data$i", data)
#     end
#     close(f)
# end
