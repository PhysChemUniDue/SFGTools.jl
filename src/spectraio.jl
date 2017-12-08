
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
function list_spectra(data_directory="./"::AbstractString;
                        exact=""::AbstractString,
                        inexact=""::AbstractString,
                        date=(0,0,0)::Tuple{Int64,Int64,Int64},
                        group=false)

  const N_PIXEL = 512
  dat = readdlm(".spectralist"; comments=false)

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
      # Strangely this gives an array of type DataArrays.DataArray{Any,1}
      # We need it to be a Bool array
      any_array =  Dates.yearmonthday.(df[:date]) .== [date]
      bool_array = convert(DataArrays.DataArray{Bool,1}, any_array)
      df = df[bool_array, :]
  end

  if group
      df = by(df, :name, df -> DataFrame(
        N = length(df[:id]),
        sizes = [unique(df[:sizes])],
        dates = [unique(round.(df[:date], Dates.Day(1)))],
        id = [df[:id]]))
  end

  return df

end

"""
Make a data file that contains information where to find spectra connected to
an ID. All spectra are internally referenced by this id. `dir` is the main data
directory (`dir/2011-01-31/Name/...`).
"""
function grab(dir="./")
    idlist = Int64[]
    namelist = String[]
    dirlist = String[]
    for (root, dirs, files) in walkdir(dir)
        for dir in dirs
            if dir == "raw"
                mlist = searchdir(joinpath(root, dir), "data.txt")
                mdict = get_metadata(joinpath(root, dir, mlist[1]))
                push!(idlist, Int64(Dates.value(DateTime(mdict["timestamp"]))))
                push!(namelist, splitdir(splitdir(root)[1])[2])
                push!(dirlist, joinpath(root, dir))
            end
        end
    end
    println("Collected $(length(idlist)) spectra")
    writedlm(".spectralist", [idlist namelist dirlist])
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

function load_spectra(id::AbstractArray)
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
    val = dict[attr]
end

get_attribute(s::SFSpectrum, attr::AbstractString) = get_attribute(s.id, attr)

function read_as_3D(directory::AbstractString, astype=Float64)
    path = joinpath(directory)
    filelist = searchdir(path, ".tiff")
    if isempty(filelist)
        println("Could not find a .tiff file in $directory.")
        return
    end
    I = FileIO.load(joinpath(path, filelist[1]))
    C = Array{UInt16,3}(size(I,1), size(I,2), length(filelist))
    for (i, file) in enumerate(filelist)
        I = FileIO.load(joinpath(path, file))
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
  data = readdlm(path; comments=false)
  keys = data[:,1]
  values = data[:,2]

  mdict = Dict{String, Any}()
  for (i, key) in enumerate(keys)
      mdict[key] = values[i]
  end

  return mdict
end


function searchdir(directory, key)
    filter!(x->contains(x, key), readdir(directory))
end
