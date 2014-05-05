using Distributions
reload("src/nsum.jl")

#@show closest_grid_loc([1,2,3], 30)
#@show closest_grid_loc([1,15,25], 30)
#@show closest_grid_loc([60,59,61,32], 30)

const rate = 430
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

feasy(x) = (v = zeros(size(x)...); f(x,v); sum(v))

rtest = NSum.Region(3, 512)
NSum.subdivide!(rtest)


upp = 2000
D = 3
println("####################")
@time tot1,r = NSum.nsum(ftree, zeros(D).+upp)
iters = 0
evals = 0
@time tot2,r = NSum.nsum(ftree, zeros(D).+upp)
@show tot1, tot2
#iters = 0
#evals = 0
#@time tot,r = IntSum.isum(ftree, zeros(D).+upp)
#@show tot
println("With $iters iters and $evals fun evals")

allnodes = collect(r)

using PyPlot
close("all")
figure()
for n in allnodes
    if length(n.vals)!=0
        plot(n.mins[1], n.mins[2], "g.")
    end
end
grid = NSum.gen_unit_grid(zeros(D), zeros(D).+upp.-1)

vals = zeros(size(grid,1))
f(grid',vals)

imshow(reshape(vals,upp,upp), origin="lower")


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
