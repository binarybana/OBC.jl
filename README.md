# OBC

An optimal Bayesian classification library and runtime for RNA-Seq data.

## Installation Instructions
 - [Install the latest release of Julia](http://julialang.org/downloads/)
 - Start a Julia terminal and run the following commands:

```julia
Pkg.update()
Pkg.clone("git://github.com/binarybana/OBC.jl.git")
```

## Usage

You are now ready to use the OBC Julia library. The core operations look 
something like the following,
```julia
using OBC
data1,data2 = ... # your datasets as integer valued matrices (samples x genes)
d1,d2 = ... # the normalization factors for each dataset (float arrays)
cls = MPM.mpm_classifier(data1, data2, d1=d1, d2=d2)
MPM.sample(cls, 10000)
bemc = MPM.bee_e_mc(cls, (mean(d1),mean(d2)))
```

For a full example (with code to generate synthetic data) see the ``run.jl`` 
runner script.

