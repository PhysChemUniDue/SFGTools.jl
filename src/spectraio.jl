import FileIO.load

"""
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

  const N_PIXEL = 512
  dir = "./"
  
  spectrafile = joinpath(dir, ".spectralist")
  if isfile(spectrafile)
    dat = readdlm(spectrafile; comments=false)
  else
    error("The file $spectrafile does not exist.")
  end

  idlist = convert(Array{Int64}, dat[:,1])
  namelist = convert(Array{AbstractString}, dat[:,2])
  dirlist = convert(Array{AbstractString}, dat[:,3])
  datelist = DateTime[]
  sizelist = Array{Int64,1}[]
  numlist = Int64[]

  for dir in dirlist
      mfilelist = searchdir(dir, "data.txt")
      mfile = mfilelist[1]
      mdict = get_metadata(joinpath(dir, mfile))

      sz = Int64[N_PIXEL/mdict["x_binning"], N_PIXEL/mdict["y_binning"], length(mfilelist)]
      d = DateTime(mdict["timestamp"])

      push!(sizelist, sz)
      push!(datelist, d)
      push!(numlist, parse(Int64, splitdir(splitdir(dir)[1])[2]))
  end

  df = DataFrame(
    id = idlist,
    name = namelist,
    number = numlist,
    date = datelist,
    sizes = sizelist,
  )

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
Make a data file that contains information where to find spectra connected to
an ID. All spectra are internally referenced by this id. `dir` is the main data
directory (`dir/2011-01-31/Name/...`).

By default this function writes only newly found spectra to the .spectralist file.
If you want to rewrite the whole file pass `getall=true` as a keyword argument
to the function.

```julia
julia> grab(directory; getall=true)
```
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
                end
            end
        end
    end

    println("Collected $(length(idlist)) spectra. ($(length(idexisting) + length(idlist)) overall)")
    # writedlm(".spectralist", [idlist namelist dirlist])
    open(".spectralist", "a") do f
        writedlm(f, [idlist namelist dirlist])
    end
end

"""
Load the spectra specified by `id`. `id` can be an Int64 or Array{Int64}.
"""
function load_spectra(id::Int64)
  # Load grabbed data
  dir = getdir(id)

  # Load the spectrum files
  s = read_as_3D(dir, Float64)

  sfspectrum = SFSpectrum(id, s)
end

function load_spectra(id::AbstractArray{Int64})
  sfspectra = SFSpectrum[]
  for i in id
    sfspectrum = load_spectra(i)
    push!(sfspectra, sfspectrum)
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
function read_as_3D(directory::AbstractString, astype=Float64)
    path = joinpath(directory)
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


function getdir(id::Int64)
    data = readdlm(".spectralist"; comments=false)
    idx = findin(data, id)
    dir = data[idx,:][3]
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
    save_mat(filename, [s])
end

function save_mat(filename::AbstractString, s::Array{SFSpectrum})
    # Check if filename has a proper extension
    if splitext(filename)[2] != ".mat"
        filename *= ".mat"
    end

    f = matopen(filename, "w")
    
    for i = 1:length(s)

        # Put Stuff in dict
        try
            name = get_attribute(s[i], "name")
        catch
            name = get_attribute(s[i], "timestamp")
        end

        data = Dict(
            "name" => get_attribute(s[i], "name"),
            "signal" => mean(s[i]),
            "wavelength" => get_ir_wavelength(s[i]),
            "wavenumber" => get_ir_wavenumber(s[i]),
        )

        # Write to file
        write(f, "data$i", data)
    end
    close(f)
end