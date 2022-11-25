const length_delaystage = 2*0.05
const max_steps_delaystage = 51492.8
const c_light = 299792458
const sec_per_step = length_delaystage/max_steps_delaystage / c_light

"""
    Get the time delay of the visible in respect to IR pulse in ps. The time zero position musst be in the comments.

julia>  nr_delay(s::SFSpectrum)
0.94 
"""
function nr_delay(s::SFSpectrum;kwargs...)
    comment = get_comment(s)
    r_time_zero_probe = r"time\s*zero\s*probe\s*=\s*(\d{5})" 
    time_zero_probe = tryparse(Int,match(r_time_zero_probe,comment).captures[1])
        if isempty(time_zero_probe)
            delay = "No time zero in comments"
            return delay
        end
    xpos = get_metadata(s)["smc stages"]["xpos"] |> Float64
    delay = (time_zero_probe - xpos) * sec_per_step *1e12 # time delay in ps
    return round(delay,sigdigits=2)
end


