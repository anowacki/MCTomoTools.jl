# Functions to examine quality of chains

acceptance_ratios(chain::Chain) = acceptance_ratios(chain.samples)

function acceptance_ratios(samples::AbstractVector{<:RawSample})
    total = Dict{Int,Int}(i=>0 for i in keys(STEP_DESCRIPTION))
    accepted = Dict{Int,Int}(i=>0 for i in keys(STEP_DESCRIPTION))
    for sample in samples
        i = Int(sample.step)
        accepted[i] += Bool(sample.accepted) ? 1 : 0
        total[i] += 1
    end
    Dict(STEP_DESCRIPTION[key] => accepted[key]/total[key] for key in keys(total))
end
