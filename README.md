# OBC

An optimal Bayesian classification library and runtime for RNA-Seq data.
module MCMC

```julia
export Sampler, MHRecord, sample

abstract Sampler # Must support methods:
#propose
#energy
#reject

#include("samc.jl")
include("mh.jl")
#include("mpm.jl")
#include("test.jl")

end #module
```
