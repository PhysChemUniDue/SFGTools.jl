__precompile__()
module SFGTools

using DataFrames
# using MAT

include("sfspectrum.jl")
include("spectraio.jl")
include("processing.jl")
include("get_polarisationcomb.jl")
include("nr_delay.jl")

const N_PIXEL = 512
const VIS_WAVELENGTH = 513.3  # After Ethalon - October 2022

export  SFSpectrum,
        average,
        grab,
        list_spectra,
        load_spectra,
        fieldcorrection!,
        get_attribute,
        get_metadata,
        is_attribute,
        rm_events!,
        get_wavelength,
        get_wavenumber,
        get_ir_wavelength,
        get_ir_wavenumber,
        get_pump_delay,
        dlpos2t,
        get_variables,
        pixelshift,
        pixelshift!,
        findpixelshift,
        savejson,
        save_data,
        failed_sifs,
        check_sif_files,
        get_pol_comb,
        get_comment,
        get_timestamp,
        nr_delay
        # save_mat
end
