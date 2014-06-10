module OBC

using Distributions

export
    Sampler, 
    MHRecord, 
    sample, 
    #propose, 
    #energy, 
    #reject, 
    plot_all, 
    plot_trace, 
    logsum
    #set_energy_limits,

abstract Sampler

include("mh.jl")

type BinaryClassifier{T<:Sampler}
    cls1 :: T
    cls2 :: T
    mcmc1 :: MHRecord
    mcmc2 :: MHRecord
end

include("cls_utils.jl")
include("utils.jl")

end # module
