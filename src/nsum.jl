module NSum

using Iterators
using Grid
#using Cubature
#using Base.Collections
#using Distributions
import Base: ndims, show, length, collect

export nsum

type Region{K}
    mins::Vector{Int}
    vals::Vector{Float64}
    mvals::Vector{Vector{Float64}}
    subs::Vector{Region{K}}
    #mins :: NTuple{K,Int}
    #subs :: NTuple{K,Region{K}}
    level::Int
    len::Int
end

Region(k, len) = Region{k}(zeros(k), [], {Float64[]}, [], 1, len)
Region(k, len, mins) = Region{k}(mins, [], {Float64[]}, [], 1, len)

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
    if r.len == 2
        #warn("Trying to subdivide a 'unit' region!")
        return nothing
    elseif length(r.subs)>0
        for sub in r.subs
            subdivide!(sub)
        end
    else
        sizehint(r.subs, 2^ndims(r))
        for c in Counter(zeros(K).+2)
            newmin = (c .- 1)*div(r.len,2) .+ r.mins
            newr = nothing
            if sum(c) == K
                newr = Region{K}(newmin, r.vals, {Float64[]}, [], r.level+1, div(r.len,2))
            elseif sum(c) == K*2
                newr = Region{K}(newmin, r.mvals[1], {Float64[]}, [], r.level+1, div(r.len,2))
            else
                newr = Region{K}(newmin, [], {Float64[]}, [], r.level+1, div(r.len,2))
            end
            push!(r.subs, newr)
        end
        r.vals = Float64[]
        r.mvals = {Float64[]}
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

#function treesum{K<:Integer}(r::Region{K}, dims)
#function treesum(r, dims)
function treesum(r::Region, dims)
    s = zeros(dims)
    dfs(treesum_, s, r, true)
    return s
end

function treesum_(r::Region, state) 
    if r.len == 2
        state[:] = reduce(.+, state, r.mvals)
        state[:] = state .+ r.vals
    else
        volume = r.len^ndims(r)
        state[:] = state .+ (r.mvals[1] .+ r.vals)./2 .* volume
    end
end

function addifdirty(r::Region, state) 
    if length(r.vals) == 0
        push!(state[1], r.mins)
        push!(state[2], r.vals)
    end
    if length(r.mvals[1]) == 0
        if r.len == 2 # Add all points
            empty!(r.mvals)
            for c in Counter(zeros(ndims(r)).+2)
                sum(c) == ndims(r) && continue # r.mins already handled in r.val
                push!(state[1], r.mins .+ c .- 1)
                push!(r.mvals, Float64[])
                push!(state[2], r.mvals[end]) 
            end
        else # just midpoint
            push!(state[1], r.mins .+ r.len/2)
            push!(state[2], r.mvals[1])
        end
    end
end

function max_(r::Region, state) 
    # Only valid for leaf nodes
    if r.len == 2 # Ignore fully explored regions
        return
    end
    if length(state) == 0
        push!(state, r)
    else
        k = length(r.mins)
        tempr = max(r.vals[1], r.mvals[1][1]) * r.len^k
        tempstate = max(state[1].vals[1], state[1].mvals[1][1]) * state[1].len^k # FIXME CHECK
        if tempr > tempstate
            state[1] = r
        end
    end
end

function find_max(r::Region)
    s = {}
    dfs(max_, s, r, true)
    length(s) == 0 && return r # We've evaluated everything
    s[1]
end

function find_min_containing(r::Region, pt)
    while true
        for n in r.subs
            if all(0 .<= pt-n.mins .< n.len) 
                r = n
                if length(r.subs) == 0
                    return r
                end
                break # out to while loop
            end
        end
    end
end

function find_uneven_branch(r::Region)
    if length(r.vals) != 0 # we are a leaf node
        if r.len == 2 # but we can't go lower
            return(-Inf, [r])
        else
            # Get a measure of discrepancy:
            return (abs(r.vals[1] - r.mvals[1][1]), [r])
        end
    else # we are still high, so we need to recurse downwards
        # return maximum discrepancy branch:
        maxpaths = [find_uneven_branch(x) for x in r.subs]
        ind = indmax([x[1] for x in maxpaths])
        push!(maxpaths[ind][2], r)
        return maxpaths[ind]
    end
end

function nsum(fdim::Integer, f::Function, maxs; abstol=0.01, maxevals=30)
    #pseudocode:
    #create octree using closest larger point to maxs
    #subdivide a few times, adding points to evaluate to dirty list
    #evaluate in batch
    #update region values
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
    tots = zeros(fdim)
    for v in dirtylist[2]
        tots .+= v
    end
    #println("current estimate: $tots")

    count = 0
    findtype = 1
    while abs(tots[1]-1) > abstol && count < maxevals
        count += 1

        if findtype == 1
            #println("Getting max")
            #get max 
            maxr = find_max(r)
            r === maxr && break # Exhausted everything

            #get regions in other 2^D-1 directions from max
            #others = Array(Any, 2^D-1)
            #resize!(others, 0)
            #for c in Counter(zeros(D).+2)
                #newloc = (c .- 1) .* 2 .-1 .+ maxr.mins
                #any(newloc .< 0) && continue
                #push!(others, newloc)
            #end
            ##println(others)

            ##subdivide these other regions
            #rs = map(others) do x
                #find_min_containing(r, x)
            #end
            #map(subdivide!, rs)
            #map(subdivide!, rs)
            subdivide!(maxr)
            findtype $= 1
        else
            dis, maxrs = find_uneven_branch(r)
            dis === -Inf && break # Exhausted everything
            maxr = maxrs[3]
            #println("Found $dis discrepancy at pt $(maxr.mins) $(maxr.len)")
            subdivide!(maxr)
            findtype $= 1
        end

        dirtylist = ({}, {})
        dfs(addifdirty, dirtylist, r, true)
        
        if length(dirtylist[2])>0
            f(dirtylist)
        end
        fill!(tots, 0.0)
        dfs(treesum_, tots, r, true)
        #println("Tots estimate: $tots")
    end
    #println("Finished $count iters")
    tots, r
end

end
