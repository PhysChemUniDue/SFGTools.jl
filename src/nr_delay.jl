const length_delaystage = 2*0.05
const max_steps_delaystage = 51492.8
const c_light = 299792458
const sec_per_step = length_delaystage/max_steps_delaystage / c_light

"""
    Get the time delay of the visible in respect to IR pulse in ps. The time zero position must be indicated (t0) or be in the comments.

julia>  nr_delay(s::SFSpectrum;t0= nothing)
0.94 
"""
function nr_delay(s::SFSpectrum;t0 = nothing,sigdigits=2)
    comment = get_comment(s)
    r_time_zero_probe = r"time\s*zero\s*probe\s*=\s*(\d{5})" 
    matched = match(r_time_zero_probe,comment)
    if matched === nothing && t0 === nothing
        error("No time zero in comments. Add the kwarg t0 to use the function.")
    end
    if matched === nothing 
        time_zero_probe = t0
    else
        time_zero_probe = tryparse(Int,matched.captures[1])
    end
    xpos = get_metadata(s)["smc stages"]["xpos"] |> Float64
    if t0 !== nothing
        delay = (t0 - xpos) * sec_per_step *1e12 # time delay in ps
            return round(delay,sigdigits=sigdigits)
    else
        delay = (time_zero_probe - xpos) * sec_per_step *1e12 # time delay in ps
            return round(delay,sigdigits=sigdigits)
    end
end


