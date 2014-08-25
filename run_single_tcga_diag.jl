reload("src/OBC.jl")
using HDF5, JLD
using Stats
using Distributions

loc = expanduser("~/largeresearch/seq-data/rnaseq.jld")
@load loc luad_norm luad_norm_sampids luad_norm_geneids

blas_set_num_threads(1)

d = 10.0
iters = 10000
burn = 2000
thin = 50
D = 2
rounds = 1
numsamps = 100

srand(123)

samps = sample(1:size(luad_norm,2), numsamps)
genes = sample(1:size(luad_norm,1), D*rounds)

i=1
data = luad_norm[genes[(i-1)*D+1:i*D],samps]'

######################################################################
# Fit Model
######################################################################

#samplers = MPM.AMWGRecord[]
#for i=1:5
    #push!(samplers, MPM.mpm_model(data,burn=burn,thin=thin,d=d,usepriors=true))
    #sample(samplers[end],iters)
#end

dnostics = OBC.gelman_rubin(samplers)

#for i=1:5
    #sample(samplers[i],iters)
#end

dnostics2 = OBC.gelman_rubin(samplers)
