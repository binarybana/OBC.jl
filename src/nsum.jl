module NSum

using Grid
#using Cubature
#using Base.Collections
#using Distributions
import Base: ndims, show, length, collect

export nsum

typealias Location Vector{Int}
typealias Dimension Vector{Int}
#mins :: NTuple{K,Int}

type EvalPoint
    location::Location
    vals::Vector{Float64}
end

abstract Region{K}

type RegionTree{K} <: Region{K}
    root::Location
    sub1::Region{K}
    sub2::Region{K}
    level::Int
    dims::Dimension
    cutdim::Int
end

type RegionLeaf{K} <: Region{K}
    root::Location
    evals::Vector{EvalPoint}
    level::Int
    dims::Dimension
end

function _split(root, dims, level, cut_plane=1) 
    newdims = copy(dims)
    K = length(dims)
    newdims[cut_plane]/=2
    newroot = copy(root)
    newroot[cut_plane]+=newdims[cut_plane]
    #@show cut_plane, root, newroot, dims, newdims
    leaf1 = RegionLeaf{K}(root, [], level+1, newdims)
    leaf2 = RegionLeaf{K}(newroot, [], level+1, newdims)
    # This is already handled by tree_from_leaf. FIXME Delete
    #for (i,c) in enumerate(Counter(zeros(K).+2))
        ## Now generate all points for leaf1 and leaf2
        #newpt = (c .- 1) .* div(newdims,2)
        #p1 = newpt .+ leaf1.root
        #p2 = newpt .+ leaf2.root

        #push!(leaf1.evals, EvalPoint(p1,[]))
        #push!(leaf2.evals, EvalPoint(p2,[]))
    #end
    return (leaf1, leaf2)
end
split_empty_tree(r::Region, cut_plane=1) = _split(r.root, r.dims, r.level, cut_plane)

function RegionTree(dims::Dimension)
    k = length(dims)
    leaf1, leaf2 = _split(zeros(k), dims, 1)
    RegionTree{k}(zeros(k), leaf1, leaf2, 1, dims, 1)
end

ndims{K}(::Type{RegionLeaf{K}}) = K
ndims{K}(::Type{RegionTree{K}}) = K
ndims(r::Region) = ndims(typeof(r))

graball(curr,state) = push!(state, curr)
collect(r::Region) = (allnodes = Region[]; dfs(graball, allnodes, r); allnodes)

function show(io::IO, r::Region)
    if isa(r, RegionTree)
        allnodes = collect(r)
    else
        allnodes = [r]
    end

    for i=r.level:maxdepth(r)
        print(io, "Level $i: ")
        levelnodes = filter(allnodes) do x
            x.level == i
        end
        for n in sort(levelnodes, by=x->x.root, lt=lexless)
            print(io, n.root)
            if isa(n, RegionLeaf) && !isempty(n.evals[1].vals)
                print(io, n.evals[1].vals)
            end
        end
        println(io, "")
    end
end

function tree_from_leaf{K}(r::RegionLeaf{K}, cut_plane)
    leaf1, leaf2 = _split(r.root, r.dims, r.level, cut_plane)
    # Copy evals from r into leaf1 and leaf2
    #leaf1_buf = EvalPoints[]
    #leaf2_buf = EvalPoints[]
    #for (i,c) in enumerate(Counter(zeros(K).+2))
        ## first save all pre-computed points
        #if c[cut_plane] == 1 # r.evals[i] point goes to leaf1
            #push!(leaf1_buf, r.evals[i])
        #else
            #push!(leaf2_buf, r.evals[i])
        #end
    #end
    newdims = leaf1.dims
    for (i,c) in enumerate(Counter(zeros(K).+2))
        # Now generate all points for leaf1 and leaf2
        newpt = (c .- 1) .* div(newdims,2)
        p1 = newpt .+ leaf1.root
        p2 = newpt .+ leaf2.root

        p1prev = p2prev = false
        # Have we already computed p1?
        for pt in r.evals
            if p1 == pt.location
                p1prev = true
                push!(leaf1.evals, pt)
            elseif p2 == pt.location
                p2prev = true
                push!(leaf2.evals, pt)
            end
        end
        p1prev || push!(leaf1.evals, EvalPoint(p1,[]))
        p2prev || push!(leaf2.evals, EvalPoint(p2,[]))
    end
    RegionTree{K}(r.root, leaf1, leaf2, r.level, r.dims, cut_plane)
end

function pick_cut_plane(r::RegionLeaf)
    K = ndims(r)
    assert(any(r.dims .> 2))
    candidates = find(x->x>2, r.dims) # only consider directions with enough space
    if rand()<0.2 || length(r.evals) == 0 || any(x->isempty(x.vals), r.evals) # no choice but to pick randomly
        return -Inf, candidates[rand(1:length(candidates))]
    else
        discrep = -Inf
        direc = -1
        #@show r.root, length(r.evals), r.evals
        for eval in r.evals[2:end-1]
            poss_direc = findfirst(eval.location .- r.evals[1].location)
            #@show poss_direc, eval
            r.dims[poss_direc] == 2 && continue
            grad = r.evals[1].vals[end] - eval.vals[end] 
            if grad > discrep
                discrep = grad
                direc = poss_direc
            end
        end
        return discrep, direc
    end
end

function cut!(r::RegionTree)
    if isa(r.sub1, RegionTree)
        cut!(r.sub1)
    elseif any(r.sub1.dims .> 2) #RegionLeaf, so we have to cut from above
        _, cut_plane = pick_cut_plane(r.sub1)
        r.sub1 = tree_from_leaf(r.sub1, cut_plane)
    end

    # argh.. I know about DRY, but how?
    if isa(r.sub2, RegionTree)
        cut!(r.sub2)
    elseif any(r.sub2.dims .> 2) #RegionLeaf, so we have to cut from above
        _, cut_plane = pick_cut_plane(r.sub2)
        r.sub2 = tree_from_leaf(r.sub2, cut_plane)
    end
    nothing
end

function dfs(f::Function, state, r::RegionLeaf, onlyleaves=false)
    f(r, state)
end

function dfs(f::Function, state, r::RegionTree, onlyleaves=false)
    !onlyleaves && f(r, state)
    tovisit = Region[r.sub1, r.sub2]
    while length(tovisit) > 0
        curr = pop!(tovisit)
        if !onlyleaves # run on all
            f(curr, state)
        elseif isa(curr,RegionLeaf) # or only on leaves
            f(curr, state)
        end
        isa(curr, RegionTree) && push!(tovisit, curr.sub1, curr.sub2)
    end
end

function recursive_dfs(f::Function, state, r::RegionTree, onlyleaves=false)
    assert(isa(state, Dict))
    if !haskey(state, :history)
        state[:history] = Region[r]
    else
        push!(state[:history],r)
    end
    !onlyleaves && f(r, state)
    recursive_dfs(f, state, r.sub1, onlyleaves)
    recursive_dfs(f, state, r.sub2, onlyleaves)
    pop!(state[:history]) # clean up after ourselves
end

function recursive_dfs(f::Function, state, r::RegionLeaf, onlyleaves=false)
    push!(state[:history],r)
    f(r, state)
    pop!(state[:history]) # clean up after ourselves
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

function treesum(r::Region, fdim)
    s = zeros(fdim)
    dfs(treesum_, s, r, true)
    return s
end

treesum_(r::RegionTree, state) = nothing
function treesum_(r::RegionLeaf, state) 
    cellvolume = prod(r.dims) / 2^ndims(r) 
    for eval in r.evals
        state[:] += eval.vals * cellvolume
    end
end

addifdirty(r::RegionTree, state) = nothing
function addifdirty(r::RegionLeaf, state) 
    assert(length(r.evals) > 0)
    for eval in r.evals
        if length(eval.vals)==0 # Don't want to add already computed ones
            push!(state, eval)
        end
    end
end

function max_(r::RegionLeaf, state) 
    if all(r.dims .== 2) # Ignore fully explored regions
        return
    end
    k = ndims(r)
    buf = zeros(state[:fdim]) 
    treesum_(r, buf)
    if !haskey(state,:maxhistory) || buf[end] > state[:maxvalue]
        state[:maxhistory] = copy(state[:history])
        state[:maxvalue] = buf[end]
    end
end

function find_max(r::Region, fdim)
    s = {:fdim=>fdim}
    recursive_dfs(max_, s, r, true)
    !haskey(s,:maxhistory) && return r # We've evaluated everything
    return s[:maxhistory]
end

function find_min_containing(r::Region, pt)
    @assert all(0 .<= pt .- r.root .< r.dims)
    while true
        r = pt[r.cutdim] .< r.dims[r.cutdim] ? r.sub1 : r.sub2
        isa(r, RegionLeaf) && return r
    end
end

function find_uneven_branch(r::RegionTree)
    # we are still high, so we need to recurse downwards
    # return maximum discrepancy branch:
    path1 = find_uneven_branch(r.sub1)
    path2 = find_uneven_branch(r.sub2)
    maxpath = path1[1] > path2[1] ? path1 : path2
    push!(maxpath[2], r)
    return maxpath 
end

function find_uneven_branch(r::RegionLeaf)
    if all(r.dims .== 2) # but we can't go lower
        return(-Inf, Region[r])
    else
        # Get a measure of discrepancy:
        discrep,direc = pick_cut_plane(r)
        return (discrep, Region[r])
    end
end

function nsum(fdim::Integer, f::Function, maxs; abstol=0.01, maxevals=30)
    #pseudocode:
    #create octree using closest larger point to maxs
    #subdivide a few times, adding points to evaluate to dirty list
    #evaluate in batch
    #update region values
    lens = map(nextpow2,iround(maxs))
    r = RegionTree(lens)
    for i=1:3
        cut!(r)
    end
    
    dirtylist = {}
    dfs(addifdirty, dirtylist, r, true)
    f(dirtylist)
    tots = zeros(fdim)
    dfs(treesum_, tots, r, true)
    #println("current estimate: $tots")

    count = 0
    findtype = 1
    totnum = fdim==1 ? 1 : fdim - 1
    while any(abs(tots[1:totnum].-1) .> abstol) && count < maxevals
        count += 1

        if findtype == 1
            maxr = find_max(r, fdim)
            r === maxr[end] && break # Exhausted everything
            maxr = maxr[end-1]
            #println("found max: $(maxr.root) w/ dims: $(maxr.dims)")

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
            cut!(maxr)
            cut!(maxr)
            findtype $= 1
        else
            dis, maxrs = find_uneven_branch(r)
            dis === -Inf && break # Exhausted everything
            maxr = maxrs[2] # 1 is the leaf 
            #println("found uneven point: $(maxr.root) w/ dims: $(maxr.dims)")
            cut!(maxr)
            cut!(maxr)
            findtype $= 1
        end

        dirtylist = {}
        dfs(addifdirty, dirtylist, r, true)
        
        if length(dirtylist)>0
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
