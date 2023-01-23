__precompile__()
module SFGTools

using DataFrames
# using MAT

include("sfspectrum.jl")
include("spectraio.jl")
include("processing.jl")
<<<<<<< Updated upstream
=======
include("get_polarisationcomb.jl")
include("nr_delay.jl")
incldue("save2dashboard.jl")
>>>>>>> Stashed changes

const N_PIXEL = 512
const VIS_WAVELENGTH = 512.4  # According to service protocol of May 2018

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
<<<<<<< Updated upstream
        save_data
=======
        save_data,
        failed_sifs,
        check_sif_files,
        get_pol_comb,
        get_comment,
        get_timestamp,
        nr_delay,
        save_dl_scan
>>>>>>> Stashed changes
        # save_mat
end
