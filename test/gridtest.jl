using Grid
using Cubature
using Distributions
using Iterators

iters = 0
function f(x)
    global iters
    iters += 1
    prod(pdf(Poisson(3.0), floor(x)))
end

function gen_grid(mins, maxs, N=30)
    D = length(mins)
    stepsizes = ceil(float(maxs .- mins) ./N)
    ranges = [mins[i]:stepsizes[i]:maxs[i] for i=1:D]
    grid = hcat(map(collect, collect(product(ranges...)))...)'
    return map(length,ranges), stepsizes, grid
end

const tol = 0.001
const upp = 100

println("####################")
@time val = hquadrature(f, [0], [upp], abstol=tol)
println(val)
println("1D used $iters iters")

println("####################")
iters = 0
@time val = hcubature(f, zeros(2), zeros(2).+upp, abstol=tol)
println(val)
println("2D used $iters iters")

println("####################")
iters = 0
@time val = hcubature(f, zeros(3), zeros(3).+upp, abstol=tol)
println(val)
println("3D used $iters iters")

println("####################")
println("Now using Grid...")
println("####################")

valsi = nothing
iiters = 0
function fi(x)
    global iters, iiters
    v,g,h = valgradhess(valsi,x...)
    if abs(prod(g)) > 1e-15 || abs(det(h)) > 1e-15
        f(x)
    else
        iiters += 1
        valsi[x...]
    end
end

N = 20

println("####################")
_,_,pts = gen_grid([0],[upp], N)
@time vals = [f(pt) for pt in pts]
valsi = InterpGrid(vals, BCnearest, InterpQuadratic)
iters = 0
iiters = 0
@time val = hquadrature(fi, [0], [upp], abstol=tol)
println(val)
println("1D used $iters actual iters, $iiters interp iters")

println("####################")
lens,_,pts = gen_grid(zeros(2), zeros(2).+upp, N)
@time vals = [f(pts[i,:]) for i=1:size(pts,1)]
vals = reshape(vals, lens...)
valsi = InterpGrid(vals, BCnearest, InterpQuadratic)
iters = 0
iiters = 0
@time val = hcubature(fi, zeros(2), zeros(2).+upp, abstol=tol)
println(val)
println("2D used $iters actual iters, $iiters interp iters")

#println("####################")
#lens,_,pts = gen_grid(zeros(3), zeros(3).+upp, N)
#@time vals = [f(pts[i,:]) for i=1:size(pts,1)]
#vals = reshape(vals, lens...)
#valsi = InterpGrid(vals, BCnearest, InterpQuadratic)
#iters = 0
#iiters = 0
#@time val = hcubature(fi, zeros(3), zeros(3).+upp, abstol=tol)
#println(val)
#println("3D used $iters actual iters, $iiters interp iters")

hess = reshape([det(valgradhess(valsi, pts[i,:]...)[3]) for i=1:size(pts,1)], lens...);
grad = reshape([prod(valgrad(valsi, pts[i,:]...)[2]) for i=1:size(pts,1)], lens...);
