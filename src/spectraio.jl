import FileIO.load
using CSV
asdf
"""
    list_spectra(; exact="", inexact="", date::Tuple{Int64,Int64,Int64}, group=false)

List available spectra.
First `grab(datafolder)` the spectra to genereate a `.spectralist` file. You can
optionally pass the path to the `.spectralist` file in the arguments.

Filter the dataframe:
* `exact="name"`: Match the exact name of the spectrum
* `inexact="part"`: Check if the name contains part of the string
* `date=(2011, 01, 31)`: Filter for specific date
* `group=true`: Groups spectra by name

The result with some useful information is stored in a DataFrame. No spectra are
loaded hereby to save time. Load the actual spectra via `load_spectra(id)`
where `id` is `list_spectra()[:id]` See the DataFrames package for more info.
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

  df = CSV.read(spectrafile; allowmissing=:none)

  # Filter the dataframe
  if exact != ""
      df = df[lowercase.(df[:name]) .== lowercase(exact), :]
  end

  if inexact != ""
      df = df[contains.(lowercase.(df[:name]), lowercase(inexact)), :]
  end

  if !all(iszero.(date))
      df = df[Dates.yearmonthday.(df[:date]) .== [date], :]
  end

  if group
      df = by(df, :name, df -> DataFrame(
        N = length(df[:id]),
        sizes = [unique(df[:sizes])],
        dates = [unique(floor.(df[:date], Dates.Day(1)))],
        id = [df[:id]]))
  end

  sort!(df)

  return df

end

# For compatibility with older versions
function list_spectra(str::String)
    list_spectra()
end

"""
    grab(dir="./"; getall=false)

Make a data file that contains information where to find spectra connected to
an ID. All spectra are internally referenced by this id. `dir` is the main data
directory (`dir/2011-01-31/Name/...`).

By default this function writes only newly found spectra to the .spectralist file.
If you want to rewrite the whole file pass `getall=true` as a keyword argument
to the function.
"""
function grab(dir="./"; getall=false)

    # Check if the .spectralist file exists. If not create it.
    # If it exist read its contents.
    if !isfile(".spectralist") || getall
        f = open(".spectralist", "w") do f
        end
        idexisting = Int64[]
    else
        content = readdlm(".spectralist")
        idexisting = convert.(Int64, content[:,1])
    end

    idlist = Int64[]
    namelist = String[]
    dirlist = String[]
    datelist = DateTime[]
    sizelist = Array{Int64}[]
    numlist = Int64[]

    for (root, dirs, files) in walkdir(dir)
        for dir in dirs
            if dir == "raw"
                mlist = searchdir(joinpath(root, dir), "data.txt")
                mdict = get_metadata(joinpath(root, dir, mlist[1]))

                # Get ID
                id = Int64(Dates.value(DateTime(mdict["timestamp"])))
                if !any(idexisting .== id)
                    push!(idlist, id)
                    push!(namelist, splitdir(splitdir(root)[1])[2])
                    push!(dirlist, joinpath(abspath(root), dir))
                    push!(datelist, DateTime(mdict["timestamp"]))
                    push!(sizelist, Int64[N_PIXEL/mdict["x_binning"], N_PIXEL/mdict["y_binning"], length(mlist)])
                    push!(numlist, parse(Int64, splitdir(splitdir(dirlist[end])[1])[2]))
                end
            end
        end
    end

    df = DataFrame(id=idlist, name=namelist, path=dirlist, date=datelist, sizes=sizelist, number=numlist)

    open(".spectralist", "a") do f
        # writedlm(f, [idlist namelist dirlist datelist sizelist numlist])
        CSV.write(f, df)
    end

    println("Collected $(length(idlist)) spectra. ($(length(idexisting) + length(idlist)) overall)")


end

"""
    load_spectra(id::Int64)

Load the spectra specified by `id`.
"""
function load_spectra(id::Int64)
  # Load grabbed data
  dir = getdir(id)

  # Load the spectrum files
  s = read_as_3D(dir, Float64)

  sfspectrum = SFSpectrum[SFSpectrum(id, s)]
end

"""
    load_spectra(id::Array{Int64})

Load the spectra specified by `id`.
"""
function load_spectra(id::AbstractArray{Int64})
  sfspectra = Array{SFSpectrum,1}(length(id))
  for i in 1:length(id)
    sfspectrum = load_spectra(id[i])
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
    dict = get_metadata(id)
    attr in keys(dict) ? val = dict[attr] : val = nothing
    val
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
function read_as_3D(path::AbstractString, astype=Float64)
    filelist = searchdir(path, ".tiff")
    if isempty(filelist)
        println("Could not find a .tiff file in $directory.")
        return
    end
    I = load(joinpath(path, filelist[1]))
    C = Array{UInt16,3}(size(I,1), size(I,2), length(filelist))
    for (i, file) in enumerate(filelist)
        I = load(joinpath(path, file))
        C[:,:,i] = reinterpret(UInt16, I)
        C = convert(Array{astype}, C)
    end
    return C
end

"""
Get the directory of a spectrum with a given id. If the ID does not exist return an empty string.
"""
function getdir(id::Int64)
    df = CSV.read(".spectralist"; allowmissing=:none)
    idx = findin(df[:id], id)
    isempty(idx) && error("Could not find spectrum with id $id.")
    dir = df[:path][idx[1]]
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
  
  path == "" && (return Dict())

  if splitext(path)[2] != ".txt"
      mfile = searchdir(path, "data.txt")[1]
      path = joinpath(path, mfile)
  end
  data = readdlm(path, '\t'; comments=false)
  keys = data[:,1]
  values = data[:,2]

  mdict = Dict{String, Any}()
  for (i, key) in enumerate(keys)
      mdict[key] = values[i]
  end

  return mdict
end


function searchdir(directory::AbstractString, key::AbstractString)
    filter!(x->contains(x, key), readdir(directory))
end


"""
Save Spectra in a Matlab file.
"""
function save_mat(filename::AbstractString, s::SFSpectrum)
    save_mat(filename, SFSpectrum[s])
end

function save_mat(filename::AbstractString, s::Array{SFSpectrum})
    # Check if filename has a proper extension
    if splitext(filename)[end] != ".mat"
        filename *= ".mat"
    end

    f = matopen(filename, "w")

    length(s)
    
    for i = 1:length(s)

        # Put Stuff in dict
        name = ""
        try
            name = get_attribute(s[i], "name")
        catch
            name = string(s[i].id)
        end

        data = Dict()
        try 
            data = Dict(
                "name" => get_attribute(s[i], "name"),
                "signal" => s[i].s,
                "wavelength" => get_ir_wavelength(s[i]),
                "wavenumber" => get_ir_wavenumber(s[i]),
            )
        catch
            data = Dict(
                "name" => name,
                "signal" => s[i].s
            )
        end
        data

        # Write to file
        write(f, "data$i", data)
    end
    close(f)
end