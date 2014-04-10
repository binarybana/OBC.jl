using OBC, MPM
using HDF5, JLD
using PyPlot
using Distributions
using DataFrames

blas_set_num_threads(1)

kmax = 1
iters = 10000

# Load Karen's RNA-Seq counts here
loc = "data/karen-clean1.jld"
@load loc data geneids sampleids

# Read in sample key
key = readtable("data/karen-rnaseq-key.csv")
key = key[sortperm(key[:Sample]),:]

# Read in gene list (start with Karen's augmented with DE first?)
targets = readcsv("data/karen_refseq_rn5.csv", UTF8String)

sel = Bool[x in targets for x=geneids]
finalgenes = geneids[sel]
data = float(data[:,sel])

numgenes = sum(sel)

# combinatorial gene selection
comparisons = [
    (["Saline", "fish", "pectin"], ["Saline", "fish", "cellulose"]),
    (["Saline", "corn", "pectin"], ["Saline", "corn", "cellulose"]),
    (["Saline", "fish", "pectin"], ["Saline", "corn", "cellulose"])]
# FIXME Check on these ^^

numcomparisons = length(comparisons)

type MCMCResult
    gene_names :: Vector{ASCIIString}
    gene_inds :: Vector{Int}
    comparison :: Any

    bee :: Float64
    bee2 :: Float64
    brmse :: Float64
    acceptance :: (Float64, Float64)
end

num_gene_lists = binomial(numgenes,2) + binomial(numgenes,3)
println("Attempting $(binomial(numgenes,2)) + $(binomial(numgenes,3)) = $(numcomparisons*num_gene_lists) MCMC runs...")

results = MCMCResult[]

data1 = data2 = 0

for compare = comparisons
    groups = collect(groupby(key, [:treatment,:oil,:fiber]))
    # Function to get rows specified by the selector c
    getrows = c -> groups[findfirst(y->vec(DataArray(y[1,2:end])) == c, groups)].rows

    samps1 = getrows(compare[1])
    samps2 = getrows(compare[2])

    for D = 2:3
        for genes = combinations(1:numgenes, D)

            t1 = time()
            data1 = data[samps1,genes]
            data2 = data[samps2,genes]

            d = 100.0
            cov = eye(D,D) 
            mu = ones(D,kmax)
            prior = MPM.MPMPrior(D=2, kmax=kmax, kappa=10., S=eye(2)*10)
            pmoves = MPM.MPMPropMoves()
            pmoves.lammove = 0.005
            pmoves.mumove = 0.6
            pmoves.priorkappa = 100.0

            ######################################################################
            # Class 1
            ######################################################################

            start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
                cov, #sigma :: Matrix{Float64}
                ones(kmax)/kmax, #w :: Vector{Float64}
                clamp(log(data1'/d),-3.0,Inf), #lam :: Matrix{Float64}
                1) #k :: Int
            
            obj_a = MPM.MPMCls(prior, data1, deepcopy(start), pmoves, d)
            obj_a.usepriors = true

            mymh_a = MPM.MHRecord(obj_a,burn=0,thin=50)
            sample(mymh_a,iters)

            ######################################################################
            # Class 2
            ######################################################################

            start = MPM.MPMParams(mu, #mu :: Matrix{Float64}
                cov, #sigma :: Matrix{Float64}
                ones(kmax)/kmax, #w :: Vector{Float64}
                clamp(log(data2'/d),-3.0,Inf), #lam :: Matrix{Float64}
                1) #k :: Int

            obj_b = MPM.MPMCls(prior, data2, deepcopy(start), pmoves, d)
            mymh_b = MPM.MHRecord(obj_b,burn=0,thin=50)
            sample(mymh_b,iters)

            gext, (n1, n2, xstride, ystride, grid) = MPM.get_grid(vcat(data1,data2), positive=true)

            bes = Float64[]
            bee2s = Float64[]
            for i=1:5
                be,bee2 = MPM.bee_moments(grid, mymh_a.db, mymh_b.db, 20, xstride*ystride)
                be,bee2
                push!(bes, be)
                push!(bee2s, bee2)
            end
            brmse = median(bee2s) - median(bes)^2

            acc = mymh_a.count_accept / mymh_a.count_total, mymh_b.count_accept / mymh_b.count_total
            push!(results, MCMCResult(
                geneids[genes], 
                genes, 
                compare, 
                median(bes), 
                median(bee2s), 
                brmse, 
                acc))

            println("Comparison took $(time()-t1) seconds")
        end
    end
end
