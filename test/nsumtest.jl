using Distributions
#using Base.Collections
using Iterators
reload("src/nsum.jl")

#@show closest_grid_loc([1,2,3], 30)
#@show closest_grid_loc([1,15,25], 30)
#@show closest_grid_loc([60,59,61,32], 30)

function gen_grid(mins, maxs, N=30)
    D = length(mins)
    stepsizes = ceil(float(maxs .- mins) ./N)
    ranges = [mins[i]:stepsizes[i]:maxs[i] for i=1:D]
    grid = hcat(map(collect, Iterators.product(ranges...))...)'
    return map(length,ranges), stepsizes, grid
end

function gen_unit_grid(mins, maxs)
    D = length(mins)
    ranges = [mins[i]:maxs[i] for i=1:D]
    hcat(map(collect, Iterators.product(ranges...))...)'
end

rate = 1000
iters = 0
evals = 0
function f(x, vals)
    global iters,evals
    iters += 1
    evals += length(vals)
    pdfs = pdf(Poisson(rate), floor(vec(x)))
    vals[:] = vec(prod(reshape(pdfs,size(x)...), 1))
end

function ftree(x)
    #println("Computing fvals for $(x[1])")
    global iters,evals
    iters += 1
    evals += length(x[1])
    for (pt,vals)=zip(x...)
        assert(length(vals) == 0)
        pdfs = pdf(Poisson(rate), pt)
        push!(vals, prod(pdfs))
    end
end

function ftree_skinny(x)
    #println("Computing fvals for $(x[1])")
    global iters,evals
    iters += 1
    evals += length(x[1])
    for (pt,vals)=zip(x...)
        assert(length(vals) == 0)
        val = 1.0
        val *= pdf(Poisson(0.01), pt[1])
        for v in pt[2:end]
            val *= pdf(Poisson(rate), v)
        end
        push!(vals, val)
    end
end

feasy(x) = (v = zeros(size(x)...); f(x,v); sum(v))

rtest = NSum.Region(3, 512)
NSum.subdivide!(rtest)


upp = 2000
D = 3
#println("####################")
#ndivs=10
#xs=linspace(0,upp,ndivs)
#spacing=upp/ndivs
#est = spacing*sum(map(x->pdf(Poisson(rate), ifloor(x)), xs))
#println("Estimate using $ndivs ndivs is $est")

println("####################")
@time tot1,r = NSum.nsum(1, ftree_skinny, zeros(D).+upp, abstol=0.03, maxevals=10)
@show tot1
println("With $iters iters and $evals fun evals")
iters = 0
evals = 0
@time tot2,r = NSum.nsum(1, ftree_skinny, zeros(D).+upp, abstol=0.03, maxevals=100)
@show tot2
#iters = 0
#evals = 0
#@time tot,r = IntSum.nsum(ftree, zeros(D).+upp)
#@show tot
println("With $iters iters and $evals fun evals")

#allnodes = collect(r)

#using PyPlot
#close("all")
#figure()
#for n in allnodes
    #if length(n.vals)!=0
        #plot(n.mins[1], n.mins[2], "g.")
    #end
#end
#lens, steps, grid = gen_grid(zeros(D), zeros(D).+upp.-1, 50)

#vals = zeros(size(grid,1))
#f(grid',vals)

#imshow(reshape(vals,lens...), extent=[0,upp,0,upp], origin="lower")

## NOTE TO SELF: trying to get above code to work to verify that my new 
# find uneven branch code with midpoint is working


#using Cubature
#println("####################")
#@time tot1 = hcubature_v(f, zeros(D), zeros(D).+upp, abstol=0.01)
#iters = 0
#evals = 0
#@time tot2 = hcubature_v(f, zeros(D), zeros(D).+upp, abstol=0.01)
#@show tot1, tot2
#println("With $iters iters and $evals fun evals")

#iters = 0
#evals = 0
#lens, steps, grid = gen_grid(zeros(D).+100, zeros(D).+140, 20)
#vals = Array(Float64, size(grid,1))
#f(grid',vals)
#@show steps
#@show sum(vals)
