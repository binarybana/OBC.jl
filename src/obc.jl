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
include("cls_utils.jl")
include("utils.jl")

end # module
