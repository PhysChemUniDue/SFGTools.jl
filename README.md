# SFGTools Overview

## Installation
From the REPL do
```julia
julia> Pkg.clone("https://github.com/PhysChemUniDue/SFGTools.jl")
```

## Usage

*EXAMPLE NOTEBOOK PROVIDED*

### Load Spectra
Use the module via
```julia
using SFGTools
```
Generate a `.spectralist` file in your current director by executing
```julia
julia> grab("./mydatadir")
```
where `mydatadir` has the folders with the specific dates as subfolders. The `.spectralist` file contains information of all spectra located in the subfolders of `mydatadir`. These are mainly the ID of the spectrum and its location in the file structure.

To get an overview over all available spectra do
```julia
julia> df = list_spectra()
```
This returns a dataframe with some basic information about the spectra. You can also get a filtered dataframe by applying filter:

* ```exact="name"```        → Filters for the exact name of the spectrum
* ```inexact="partofname"```  → Filters for matching parts in the spectrums name
* ```d=(2011, 01, 31)```      → Filters for a specific date

Example:
```julia
julia> df = list_spectra(inexact="Au", d=(2011, 11, 1))
```
To get all the spectra whose name contains "Au" from the 1st of November 2011.

To load the spectra you want you have to get their IDs. You can do this manually or get them from the filtered dataframe.
```julia
julia> id = df[:id]
julia> data = load(id)
```

### SFSpectrum Objects
The loaded `data` is an `SFSpectrum` object or an array of `SFSpectrum` objects depending on if `id` was a single value or an array. These have two fields:

* `id`: The unique ID of the spectrum (this is built from the timestamp of the first spectrum in the series)
* `s`: The spectrum as a 3D Matrix (the first dimension corresponds to the horizontal on the CCD, the second dimension to the vertical axis of the CCD and the third dimension to the number of the spectrum in the series)

The `id` enables us to load metadata of the spectrum easily by executing
```julia
julia> get_metadata(data[1])
```
or getting attributes of every spectrum in the series:
```julia
julia> get_attribute(data, "ccd_temperature")
```

You can get get the IR wavelengths/wavenumbers corresponding to the pixes via
```julia
julia> get_ir_wavelength(data[1])
julia> get_ir_wavenumber(data[1])
```

### Processing
See:

* `rm_background!`
* `rm_events!`
