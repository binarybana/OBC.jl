module OBC

using Distributions
using Iterators
using SAMC

export
    plot_all, 
    plot_trace, 
    logsum,
    get_bbox,
    gen_grid,
    gelman_rubin, 
    get_bbox, 
    gen_grid, 
    gen_unit_grid

type BinaryClassifier{T<:Sampler,S<:MCMC}
    cls1 :: T
    cls2 :: T
    mcmc1 :: S
    mcmc2 :: S
end

include("cls_utils.jl")
include("utils.jl")

end # module
