module OBC

using Distributions
using Iterators

export
    Sampler, 
    MHRecord, 
    AMWGRecord, 
    sample, 
    #propose, 
    #energy, 
    #reject, 
    plot_all, 
    plot_trace, 
    logsum,
    get_bbox,
    gen_grid

    #set_energy_limits,

abstract Sampler
abstract MCMC

include("mh.jl")
include("amwg.jl")

type BinaryClassifier{T<:Sampler,S<:MCMC}
    cls1 :: T
    cls2 :: T
    mcmc1 :: S
    mcmc2 :: S
end

include("cls_utils.jl")
include("utils.jl")

end # module
