module IntSum

using Iterators
using Grid
using Cubature
using Base.Collections
using Distributions
import Base: ndims, show, length, collect

function gen_grid(mins, maxs, N=30)
    D = length(mins)
    stepsizes = ceil(float(maxs .- mins) ./N)
    ranges = [mins[i]:stepsizes[i]:maxs[i] for i=1:D]
    grid = hcat(map(collect, collect(product(ranges...)))...)'
    return map(length,ranges), stepsizes, grid
end

function gen_unit_grid(mins, maxs)
    D = length(mins)
    ranges = [mins[i]:maxs[i] for i=1:D]
    hcat(map(collect, collect(product(ranges...)))...)'
end

function closest_grid_loc(pt, N)
    divs = div(pt,N)
    rems = pt%N
    smaller = rems .< div(N,2)
    N .* (divs .+ [smaller[i]?0:1 for i=1:length(pt)])
end

#function intsum_(f, mins, maxs, vals)
    ## Evaluate f on mins, maxs, using unit grid and preallocated vals vector
    #grid = gen_unit_grid(mins, maxs)
    #resize!(vals, size(grid,1))
    #f(grid', vals)
    #return sum(vals)
#end

function intsum(f::Function, mins, maxs)
    # Pseudo code:
    # take a dirty first grid
    # get the max point
    # in all directions around this ('grid' centered) point:
    #   compute exact integral at integer points
    #   return this summation
    #
    #   eventually:
    #   keep expanding in all directions where sum is above a threshold
    #   stop once all boundary regions are below this threshold
    
    assert(length(mins) == length(maxs))
    D = length(mins)
    maxN = iround(1000^(1/D))
    lens, steps, grid = gen_grid(mins, maxs, maxN)
    vals = Array(Float64, size(grid,1))
    f(grid',vals)
    #@show steps
    #@show sum(vals)
    if all(steps .== 1)
        return sum(vals)
    end
    
    maxind = indmax(vals)
    maxpt = grid[maxind,:]
    #@show maxpt
    maxsnap = int(closest_grid_loc(maxpt, maxN))

    pq = PriorityQueue(Base.Order.Reverse)
    #memory = (NTuple{D,Int}=>Float64)[]
    memory = Set{Vector{Int}}()
    #memory = Set{NTuple{D,Int}}()

    cubeitr = product(repeated(0:1,D)...)
    minlist = {[bool(i)?v-maxN:v for (i,v) in zip(itr,maxsnap)] for itr=cubeitr}
    maxlist = {[bool(i)?v-1:v+maxN-1 for (i,v) in zip(itr,maxsnap)] for itr=cubeitr}
    tot = 0.0
    maxval = 0.0
    for (mins,maxs) in zip(minlist,maxlist)
        if any(mins.<0)
            continue
        end
        #@show mins, maxs
        grid = gen_unit_grid(mins, maxs)
        resize!(vals, size(grid,1))
        f(grid', vals)
        subtot = sum(vals)
        tot += subtot
        maxval = max(maxval, subtot)
        #memory[mins] = subtot
        push!(memory, mins)
        pq[mins] = subtot
    end

    #println(pq)
    while length(pq) > 0
        (mins,val) = peek(pq)
        dequeue!(pq)
        if val > 0.1*maxval
            minlist = {Int[bool(i)?v-maxN:v+maxN for (i,v) in zip(itr,mins)] for itr=cubeitr}
            maxlist = {Int[bool(i)?v-1:v+2*maxN-1    for (i,v) in zip(itr,mins)] for itr=cubeitr}
            for (mins,maxs) in zip(minlist,maxlist)
                if any(mins.<0) || mins in memory
                    continue
                end
                grid = gen_unit_grid(mins, maxs)
                resize!(vals, size(grid,1))
                f(grid', vals)
                subtot = sum(vals)
                #@show mins, maxs, subtot
                tot += subtot
                maxval = max(maxval, subtot)
                push!(memory, mins)
                pq[mins] = subtot
            end
        end
    end
    tot
end

type Region{K}
    #mins :: NTuple{K,Int}
    mins::Vector{Int}
    vals::Vector{Float64}
    #subs :: NTuple{K,Region{K}}
    subs::Vector{Region{K}}
    level::Int
    len::Int
end

Region(k, len) = Region{k}(zeros(k), [], [], 1, len)
Region(k, len, mins) = Region{k}(mins, [], [], 1, len)

ndims{K}(::Type{Region{K}}) = K
ndims(r::Region) = ndims(typeof(r))

graball(curr,state) = push!(state, curr)
collect(r::Region) = (allnodes = Region[]; dfs(graball, allnodes, r); allnodes)

function show(io::IO, r::Region)
    allnodes = collect(r)

    for i=1:maxdepth(r)
        print(io, "Level $i: ")
        levelnodes = filter(allnodes) do x
            x.level == i
        end
        for n in sort(levelnodes, by=x->x.mins, lt=lexless)
            print(io, n.mins)
            if length(n.vals)!=0
                @printf io ":%.2f, " n.vals[1]
            end
        end
        println(io, "")
    end
end

function subdivide!{K}(r::Region{K})
    if r.len == 1
        warn("Trying to subdivide a unit region!")
    elseif length(r.subs)>0
        for sub in r.subs
            subdivide!(sub)
        end
    else
        sizehint(r.subs, 2^ndims(r))
        for c in Counter(zeros(K).+2)
            newmin = (c .- 1)*div(r.len,2) .+ r.mins
            newr = Region{K}(newmin, sum(c)==K ? r.vals : [], 
                [], r.level+1, div(r.len,2))
            push!(r.subs, newr)
        end
        r.vals = []
    end
    nothing
end

function dfs(f::Function, state, r::Region, onlyleaves=false)
    !onlyleaves && f(r, state)
    tovisit = copy(r.subs)
    while length(tovisit) > 0
        curr = pop!(tovisit)
        if onlyleaves && length(curr.subs)==0 # only run on leaves
            f(curr, state)
        elseif !onlyleaves
            f(curr, state)
        end
        append!(tovisit, curr.subs)
    end
end

function maxdepth_(r::Region, state) 
    if r.level > state[1] 
        state[1] = r.level
    end
end

function minlen_(r::Region, state) 
    if r.len < state[1] 
        state[1] = r.len
    end
end

length_(r::Region, state) = state[1] += 1
length(r::Region) = (state = {0}; dfs(length_, state, r); state[1])

maxdepth(r::Region) = (state = {0}; dfs(maxdepth_, state, r); state[1])
minlen(r::Region) = (state = {typemax(Int)}; dfs(minlen_, state, r); state[1])

function treesum(r::Region)
    if length(r.vals)!=0
        return r.vals
    else
        s = zeros(1)
        dfs(treesum, s, r, true)
        return s
    end
end

function treesum(r::Region, state) 
    state[:] = state .+ r.vals .* r.len^ndims(r)
end

function addifdirty(r::Region, state) 
    if length(r.vals) == 0
        push!(state[1], r.mins)
        push!(state[2], r.vals)
    end
end

function find_uneven_branch(r::Region)
    if length(r.subs[1].vals) != 0 # we are right above leaf nodes
        if r.subs[1].len == 1 # but we can't go lower
            return(-Inf, r)
        else
            # Get a measure of discrepancy:
            vals = [x.vals[1] for x=r.subs]
            ind = indmax(vals)
            #return (std(vals), r.subs[ind])
            return (sum(vals), r)
        end
    else # we are still high, so we need to recurse downwards
        # return maximum discrepancy branch:
        maxpaths = [find_uneven_branch(x) for x in r.subs]
        ind = indmax([x[1] for x in maxpaths])
        return maxpaths[ind]
    end
end

function intsum2(f::Function, maxs)
    #pseudocode:
    #create octree using closest larger point to maxs
    #subdivide a few times, adding points to evaluate to dirty list
    #evaluate in batch
    #update region values
    #
    D = length(maxs)
    len = nextpow2(maximum(iround(maxs)))
    r = Region(D, len)
    subdivide!(r)
    subdivide!(r)
    subdivide!(r)
    #subdivide!(r)
    #subdivide!(r)
    
    dirtylist = ({}, {})
    dfs(addifdirty, dirtylist, r, true)
    f(dirtylist)
    tots = zero(1)
    for v in dirtylist[2]
        tots .+= v
    end
    println("current estimate: $tots")

    count = 0
    while abs(tots[1]-1) > 0.01 && count < 100 #maximum(tots[1:2])
        count += 1
        # find good one
        stdcurr, curr = find_uneven_branch(r)
        println("Chose a branch with std: $stdcurr, branch: $curr")
        # subdivide a few times
        subdivide!(curr)
        #subdivide!(curr)
        # compute again
        dirtylist = ({}, {})
        dfs(addifdirty, dirtylist, curr, true)
        f(dirtylist)
        fill!(tots, 0.0)
        dfs(treesum, tots, r, true)
        println("Tots estimate: $tots")
    end
    tots, r
end

#@show closest_grid_loc([1,2,3], 30)
#@show closest_grid_loc([1,15,25], 30)
#@show closest_grid_loc([60,59,61,32], 30)

const rate = 60
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
    println("Computing fvals for $(x[1])")
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

#rtest = Region(3, 512)
#subdivide!(rtest)


upp = 128
D = 2
println("####################")
#@time tot1,r = intsum2(ftree, zeros(D).+upp)
#iters = 0
#evals = 0
#@time tot2,r = intsum2(ftree, zeros(D).+upp)
#@show tot1, tot2
iters = 0
evals = 0
@time tot,r = intsum2(ftree, zeros(D).+upp)
@show tot
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
grid = gen_unit_grid(zeros(D), zeros(D).+upp.-1)

vals = zeros(size(grid,1))
f(grid',vals)

imshow(reshape(vals,upp,upp), origin="lower")


#println("####################")
#@time tot1 = intsum(f, zeros(D), zeros(D).+upp)
#iters = 0
#evals = 0
#@time tot2 = intsum(f, zeros(D), zeros(D).+upp)
#@show tot1, tot2
#println("With $iters iters and $evals fun evals")

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

end
