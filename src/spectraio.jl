using FileIO

"""
List available spectra.
The working directory has to be like the following:
`pdw()/2011-11-11/SpectrumName/...``

The result with some useful information is stored in a DataFrame. No spectra are
loaded hereby to save time. Load the actual spectra via `load_spectra(ids)`.
"""
function list_data(data_directory="./"; filterkey="")
  # Lets put the result in a DataFrame for nice displaying and
  # nice functionality like sorting

  # Initialize things that end up in the data frame
  ids = Int64[]        #defined by first spectrums timestamp
  names = String[]
  numbers = Int64[]    #subfolder number
  dates = DateTime[]
  sizes = Array{Int64,1}[]
  directories = String[]

  # List only directories
  dir_list_dates = filter!(x -> isdir(joinpath(data_directory, x)),
                           readdir(data_directory))

  for date_dir in dir_list_dates
    # List those directories with that at least contain a folder named "1".
    # Data found is the directory
    # with the name that was given to the spectrum.
    name_folders = filter!(x -> ispath(joinpath(data_directory, date_dir, x, "1")),
                         readdir(joinpath(data_directory, date_dir)))

    for name_folder in name_folders
      # Loop through the numberd folders
      num_folders = readdir(joinpath(data_directory, date_dir, name_folder))

      for num_folder in num_folders
        # Search for txt files which hold metadata information
        path = joinpath(data_directory, date_dir, name_folder, num_folder)
        meta_file_list = searchdir(joinpath(path, "raw"), ".txt")

        # Check if theres metadata in the folder. If not continue with next
        # folder.
        if isempty(meta_file_list)
            println("Didn't find metadata in $path.")
            continue
        end

        meta_dict = read_metadata(joinpath(path, "raw", meta_file_list[1]))
        size_spectra = Int64[512/meta_dict["x_binning"], 512/meta_dict["y_binning"], length(meta_file_list)]
        date = DateTime(meta_dict["timestamp"])

        push!(ids, Int64(Dates.value(date)))
        push!(names, name_folder)
        push!(numbers, parse(Int64, num_folder))
        push!(dates, date)
        push!(sizes, size_spectra)
        push!(directories, path)
      end
    end
  end

  df = DataFrame(
    id = ids,
    name = names,
    number = numbers,
    date = dates,
    siz = sizes,
    dir = directories,
  )

  # Filter the dataframe
  if filterkey != ""
      df = df[lowercase.(df[:name]) .== lowercase(filterkey), :]
  end

  return df

end


function load_spectrum(df::DataFrame, id::Int64)
  # Get the row we are interested in
  subdf = df[df[:id] .== id, :]

  # Get the stuff that was already in the data frame
  name = subdf[:name][1]
  directory = subdf[:dir][1]

  # Load the spectrum files
  spectrum = read_as_3D(directory, Float64)

  sfspectrum = SFSpectrum(id, name, directory, spectrum)
end

function load_spectrum(df::DataFrame, id::AbstractArray)
  sfspectra = SFSpectrum[]
  for i in id
    # sfspectrum = load_spectrum(df, i)
    push!(sfspectra, sfspectrum)
  end
end


function read_as_3D(directory::AbstractString, astype=Float64)
    path = joinpath(directory, "raw")
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


function read_metadata(meta_file)
  data = readdlm(meta_file, '\t')
  keys = data[:,1]
  values = data[:,2]

  meta_dict = Dict{String, Any}()
  for (i, key) in enumerate(keys)
      meta_dict[key] = values[i]
  end

  return meta_dict
end


function searchdir(directory, key)
    filter!(x->contains(x, key), readdir(directory))
end
