
type MHRecord
    obj :: Sampler
    mapvalue :: Sampler
    mapenergy :: Float64
    db :: Vector{Any}

    count_accept :: Int
    count_total :: Int
    iteration :: Int

    burn :: Int
    thin :: Int
    verbose :: Bool

    scheme_propose :: Dict{Int,Int}
    scheme_accept :: Dict{Int,Int}
end

MHRecord(obj::Sampler; burn=0, thin=1, verbose=false) = MHRecord(obj, 
                                obj, #mapvalue
                                Inf, #mapenergy
                                Any[], #db
                                0,0,1, #accept, total, iteration
                                burn,thin, #burn, thin
                                verbose, #verbose
                                (Int=>Int)[],
                                (Int=>Int)[])

import Distributions.sample

function sample(rec::MHRecord, iters::Int; verbose=false)
    oldenergy = energy(rec.obj)
    #dbsize = div(rec.iteration + iters - rec.burn, rec.thin)
    #dbsize = dbsize > 0 ? dbsize : 0
    #init_db(rec, dbsize) # FIXME
    println("Initial energy: $oldenergy")

    for current_iter = rec.iteration:(rec.iteration+iters-1)
        scheme = propose(rec.obj)
        rec.scheme_propose[scheme] = get(rec.scheme_propose, scheme, 0) + 1
        newenergy = energy(rec.obj)

        if newenergy < rec.mapenergy #I need to decide if I want this or not
            rec.mapenergy = newenergy
            rec.mapvalue = deepcopy(rec.obj)
        end

        ### Acceptance of new moves ###
        r = oldenergy - newenergy
        if r > 0.0 || rand() < exp(r) #Accept
            rec.count_accept += 1
            oldenergy = newenergy
            rec.scheme_accept[scheme] = get(rec.scheme_accept, scheme, 0) + 1
            if verbose
                println("A: old: $oldenergy, new: $newenergy, diff: $(oldenergy-newenergy)")
            end
        else #Reject
            reject(rec.obj)
            if verbose
                println("R: old: $oldenergy, new: $newenergy, diff: $(oldenergy-newenergy)")
            end
        end
        rec.count_total += 1

        if rec.iteration >= rec.burn && rec.iteration%rec.thin == 0
            push!(rec.db, deepcopy(rec.obj.curr))
        end
        
        if (rec.iteration) % 1000 == 0
            @printf "Iteration: %8d, best energy: %7f, current energy: %7f\n" rec.iteration rec.mapenergy oldenergy
        end
        rec.iteration += 1
    end
    #save_state_db(rec, temperature) #FIXME
    println("Accepted samples: $(rec.count_accept)")
    println("Total samples: $(rec.count_total)")
    println("Acceptance: $(rec.count_accept/rec.count_total)")
end

# To be overwritten
energy(x::Sampler) = error("Must be instantiated for your sampler object")
propose(x::Sampler) = error("Must be instantiated for your sampler object")
reject(x::Sampler) = error("Must be instantiated for your sampler object")
