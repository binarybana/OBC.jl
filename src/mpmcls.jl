# Everything dealing with TWO classifier objects

function mpm_classifier(data1, data2; burn=1000, thin=50, d1=100.0, d2=100.0, kappa=10.0, usepriors=true)
    assert(size(data1,2) == size(data1,2))
    D = size(data1, 2)
    cov = eye(D,D) .* 0.1
    mu = ones(D)
    prior = MPM.MPMPrior(D=D, kappa=kappa, S=eye(D)*0.05*10)

    # Class 1
    start = MPM.MPMParams(mu, cov, 
        clamp(log(data1./d1)',-8.0,Inf)) #lam 
    obj_a = MPM.MPMSampler(prior, data1, deepcopy(start), d1)
    obj_a.usepriors = usepriors
    #mymh_a = MPM.MHRecord(obj_a,burn=burn,thin=thin)
    mymh_a = MPM.AMWGRecord(obj_a, blocks, burn=burn,thin=thin)

    # Class 2
    start.lam = clamp(log(data2./d2)',-8.0,Inf) # lam
    obj_b = MPM.MPMSampler(prior, data2, deepcopy(start), d2)
    obj_b.usepriors = usepriors
    #mymh_b = MPM.MHRecord(obj_b,burn=burn,thin=thin)
    mymh_b = MPM.AMWGRecord(obj_b, blocks, burn=burn,thin=thin)

    OBC.BinaryClassifier(obj_a, obj_b, mymh_a, mymh_b)
end

function sample(cls::OBC.BinaryClassifier, iters=10000; verbose=false)
    t0 = time()
    OBC.sample(cls.mcmc1, iters, verbose=verbose)
    t1 = time()
    OBC.sample(cls.mcmc2, iters, verbose=verbose)
    t2 = time()
    return t1-t0, t2-t1
end

import Base: vcat 

function vcat(xs::OBC.BinaryClassifier...)
    x = xs[1]
    for y in xs[2:end]
        append!(x.mcmc1.db, y.mcmc1.db)
        append!(x.mcmc2.db, y.mcmc2.db)
    end
    return x
end

function gelman_rubin(samplers::Vector{OBC.BinaryClassifier})
    return [gelman_rubin([x.mcmc1 for x in samplers]), gelman_rubin([x.mcmc1 for x in samplers])]
end

function acceptance_rates(cls::OBC.BinaryClassifier)
    mc1,mc2 = cls.mcmc1, cls.mcmc2
    accs = Dict{String,Any}()
    accs["1"] = cls.mcmc1.block_accept
    accs["2"] = cls.mcmc2.block_accept
    return accs
end

function bee_moments(cls::OBC.BinaryClassifier; dmean=10.0, numlam=20, numpts=20)
    bee_moments(cls, (dmean,dmean), numlam=numlam, numpts=numpts)
end

function bee_moments(cls::OBC.BinaryClassifier, dmeans::NTuple{2,Float64}; numlam=20, numpts=20)
    # FIXME This is only valid for c=0.5
    dmean1,dmean2 = dmeans

    db1 = cls.mcmc1.db
    db2 = cls.mcmc2.db
    assert(length(db1) == length(db2))
    dblen = length(db1)
    dim = length(db1[1].mu)
    points = Array(Int, dim, numpts*2)
    lamsx1 = Array(Float64, dim, numpts)
    lamsx2 = Array(Float64, dim, numpts)
    lams1 = Array(Float64, dim, numlam)
    lams2 = Array(Float64, dim, numlam)
    llams1 = Array(Float64, dim, numlam)
    llams2 = Array(Float64, dim, numlam)
    rlams1 = Array(Float64, dim, numlam)
    rlams2 = Array(Float64, dim, numlam)

    beem1 = beem2 = 0.0

    for i in 1:dblen 
        dbpass = 0.0
        #dbpass = 0.0
        curr1 = db1[i]
        curr2 = db2[i]
        # Generate x values (from lams or from new lams?)
        rand!(MultivariateNormal(curr1.mu, curr1.sigma), lamsx1)
        rand!(MultivariateNormal(curr2.mu, curr2.sigma), lamsx2)

        if any(lamsx1 .> 35) || any(lamsx2 .> 35)
            #warn("Encountered lambda value greater than 35, skipping iteration")
            continue
        end

        for i=1:dim, j=1:numpts
            points[i,j] = rand(Poisson(dmean1*exp(lamsx1[i,j])))
            points[i,j+numpts] = rand(Poisson(dmean2*exp(lamsx2[i,j])))
        end

        g1 = calc_g(points', db1, numlam, dmean=dmean1)
        g2 = calc_g(points', db2, numlam, dmean=dmean2)
        g1 = clamp(g1, -10_000_000.0, Inf)
        g2 = clamp(g2, -10_000_000.0, Inf)

        # Now compute p(y|x,\theta)
        rand!(MultivariateNormal(curr1.mu, curr1.sigma), lams1)
        rand!(MultivariateNormal(curr2.mu, curr2.sigma), lams2)
        for k=1:numlam*dim
            rlams1[k] = dmean1 * exp(lams1[k])
            llams1[k] = log(dmean1) + lams1[k]
            rlams2[k] = dmean2 * exp(lams2[k])
            llams2[k] = log(dmean2) + lams2[k]
        end
        for j in 1:numpts*2 # each point
            accumlam1 = 0.0
            accumlam2 = 0.0
            #accumlam1 = -Inf
            #accumlam2 = -Inf
            for s in 1:numlam # each lambda
                accumD1 = 0.0
                accumD2 = 0.0
                for d in 1:dim # each dimension
                    accumD1 += points[d,j]*llams1[d,s] - rlams1[d,s] - lgamma(points[d,j]+1)
                    accumD2 += points[d,j]*llams2[d,s] - rlams2[d,s] - lgamma(points[d,j]+1)
                end
                accumlam1 += exp(accumD1)
                accumlam2 += exp(accumD2)
                #accumlam1 = accumlam1 == -Inf ? accumD1 : logsum(accumlam1, accumD1)
                #accumlam2 = accumlam2 == -Inf ? accumD2 : logsum(accumlam2, accumD2)
            end
            temp = g1[j] < g2[j] ? accumlam1 : accumlam2 #FIXME c!=0.5
            z = accumlam1 + accumlam2 #FIXME c!=0.5
            if z != 0.0
                dbpass += temp/z #FIXME c!=0.5
            end
            #z = logsum(accumlam1, accumlam2)
            #@show exp(temp-z)
            #dbpass = dbpass == -Inf ? temp-z : logsum(dbpass,temp-z) #FIXME c!=0.5
            #dbpass += exp(temp-z)
        end
        temp = dbpass/numpts/2 #FIXME c!=0.5
        beem1 += temp
        beem2 += temp^2
    end
    beem1 /= dblen
    beem2 /= dblen
    return beem1, beem2
end

function bee_e_eff(g0, g1, volume)
    # FIXME This is only valid for c=0.5
    # Calculate expectation of error given the two effective class conditional
    # densities (posterior predictive densities)
    exp(logsum(min(g0,g1) .+ log(volume)))*0.5
end

function bee_e_grid(points, db1, db2, numlam, volume; dmean=10.0)
    # FIXME This is only valid for c=0.5
    g1 = calc_g(points, db1, numlam, dmean=dmean)
    g2 = calc_g(points, db2, numlam, dmean=dmean)
    bee_e_eff(g1,g2,volume)
end

function bee_e_data_grid(data, db1, db2, numlam; dmean=10.0)
    mins,maxs = get_bbox(data)
    lens,steps,points = gen_unit_grid(mins,maxs)
    println("Num points: $(size(points)), maxs: $maxs")
    bee_e_grid(points, db1, db2, numlam, prod(steps))
end

function bee_e_data(data, db1, db2, numlam; dmean=10.0, maxtry=10)
    # FIXME This is only valid for c=0.5
    local volume, points, g1, g2
    g1sum = g2sum = 0.0
    factor = 0.0
    D = size(data,2)
    maxN = iround(20_000^(1/D))
    trycount = 1
    numpts = 0
    while max(abs(g1sum-1.0),abs(g2sum-1.0)) > 0.1
        mins, maxs = get_bbox(data, factor=factor)
        lens, steps, points = gen_grid(mins, maxs, maxN)
        volume = prod(steps)
        g1 = calc_g(points, db1, numlam, dmean=dmean)
        g2 = calc_g(points, db2, numlam, dmean=dmean)
        g1sum = exp(logsum(g1 .+ log(volume)))
        g2sum = exp(logsum(g2 .+ log(volume)))
        println("g1 sum: $g1sum")
        println("g2 sum: $g2sum")
        numpts += length(points)
        #println("steps: $steps")
        #println(maxs, " ", size(points), " ", factor)
        factor += 1.0
        if trycount > maxtry
            return -min(g1sum, g2sum)
        else
            trycount += 1
        end
    end
    println("bee_e_data used $trycount iterations, and $numpts * numclasses evaluations")
    bee_e_eff(g1,g2,volume)
end

function bee_e_mc_naive(cls::OBC.BinaryClassifier; dmean=10.0, numpts=100)
    #FIXME only c=0.5
    #Generate points from each class
    pts1 = gen_posterior_points(numpts, dmean, cls.mcmc1.db)
    pts2 = gen_posterior_points(numpts, dmean, cls.mcmc2.db)
    acc_numpts = size(pts1,2) + size(pts2,2) # requested != generated

    #Now classify and count up mistakes
    points = hcat(pts1,pts2)
    labels = [zeros(Float64,size(pts1,2)), ones(Float64,size(pts2,2))]
    errs = abs(predict(cls.mcmc1.db, cls.mcmc2.db, points', dmean=dmean) .- labels)
    return sum(errs)/acc_numpts
end

function bee_e_mc(cls::OBC.BinaryClassifier; dmean=10.0, numpts=100)
    bee_e_mc(cls, (dmean,dmean), numpts=numpts)
end

function bee_e_mc(cls::OBC.BinaryClassifier, dmeans::NTuple{2,Float64}; numpts=100)
    #FIXME only c=0.5
    #Generate points from each class
    dmean1,dmean2 = dmeans
    pts1 = gen_posterior_points(numpts, dmean1, cls.mcmc1.db)
    pts2 = gen_posterior_points(numpts, dmean2, cls.mcmc2.db)
    acc_numpts = size(pts1,2) + size(pts2,2) # requested != generated
    points = hcat(pts1,pts2)
    g0 = calc_g(points', cls.mcmc1.db, 20, dmean=dmean1)
    g1 = calc_g(points', cls.mcmc2.db, 20, dmean=dmean2)
    g0 = clamp(g0, -10_000_000.0, Inf)
    g1 = clamp(g1, -10_000_000.0, Inf)
    z = mapslices(logsum, hcat(g0,g1), 2)
    res = exp(logsum(min(g0,g1) .- z .- log(acc_numpts)))
    return res
end

bee_moments_2sample(cls::OBC.BinaryClassifier; dmean=10.0, numpts=100) = bee_moments_2sample(cls, (dmean,dmean), numpts=numpts)

function bee_moments_2sample(cls::OBC.BinaryClassifier, dmeans::NTuple{2,Float64}; numpts=100)
    #FIXME only c=0.5
    #Generate points from each class
    dmean1,dmean2 = dmeans
    pts1a = gen_posterior_points(numpts, dmean1, cls.mcmc1.db)
    pts1b = gen_posterior_points(numpts, dmean1, cls.mcmc1.db)
    pts2a = gen_posterior_points(numpts, dmean2, cls.mcmc2.db)
    pts2b = gen_posterior_points(numpts, dmean2, cls.mcmc2.db)

    @assert size(pts1a,2) == size(pts1b,2)
    @assert size(pts2a,2) == size(pts2b,2)

    acc_numpts = size(pts1a,2) + size(pts2a,2) # requested != generated
    pointsa = hcat(pts1a,pts2a)
    pointsb = hcat(pts1b,pts2b)
    g1a = calc_g(pointsa', cls.mcmc1.db, 20, dmean=dmean1)
    g1b = calc_g(pointsb', cls.mcmc1.db, 20, dmean=dmean1)
    g2a = calc_g(pointsa', cls.mcmc2.db, 20, dmean=dmean2)
    g2b = calc_g(pointsb', cls.mcmc2.db, 20, dmean=dmean2)

    g1a = clamp(g1a, -10_000_000.0, Inf)
    g1b = clamp(g1b, -10_000_000.0, Inf)
    g2a = clamp(g2a, -10_000_000.0, Inf)
    g2b = clamp(g2b, -10_000_000.0, Inf)

    #agree = ~((g1a .< g2a) $ (g1b .< g2b))
    agree(x,y) = ~(x $ y)
    labela = (g1a .< g2a)
    labelb = (g1b .< g2b)

    za = [logsum(g1a[i],g2a[i]) for i=1:acc_numpts]
    zb = [logsum(g1b[i],g2b[i]) for i=1:acc_numpts]

    lpyx = min(g1a,g2a) .- za
    lpyz = min(g1b,g2b) .- zb

    ## BEE calculation
    bee = exp(logsum(vcat(lpyx,lpyz)) - log(acc_numpts*2))

    ## BEE second moment
    bee2 = -Inf
    for i=1:acc_numpts, j=1:acc_numpts
        if agree(labela[i],labelb[j])
            bee2 = logsum(bee2, lpyx[i] + lpyz[j])
        end
    end
    bee2 = exp(bee2 - log(acc_numpts))

    return bee, bee2, bee2 - bee^2
end

function bee_e_nsum(cls::OBC.BinaryClassifier, numlam; dmean=10.0, abstol=0.03, maxevals=30)
    # FIXME This is currently only valid for c=0.5
    data = vcat(cls.cls1.data, cls.cls2.data)
    db1,db2 = cls.mcmc1.db, cls.mcmc2.db
    mins, maxs = get_bbox(data, factor=2)
    assert(length(db1) == length(db2))
    global iters = 0
    global evals = 0
    function error_1st(data)
        tabpoints = hcat([x.location for x in data]...)
        #points: dxn array to evaluate at
        #vals: 2xn values to store into
        numpts = length(data)
        iters += 1
        evals += numpts
        g1 = calc_g(tabpoints', db1, numlam, dmean=dmean)
        g2 = calc_g(tabpoints', db2, numlam, dmean=dmean)
        # For c != 0.5, we'll need to multiply by c depending on which class is
        # less
        errpts = exp(min(g1,g2))
        for (i,pt) in enumerate(data)
            resize!(pt.vals, 3)
            pt.vals[1] = exp(g1[i])
            pt.vals[2] = exp(g2[i])
            pt.vals[3] = errpts[i]
        end
    end
    tot1,r = NSum.nsum(3, error_1st, maxs, abstol=abstol, maxevals=maxevals)
    #println("maxs: $maxs")
    #println("NSum used $iters iterations, and $evals * numclasses evaluations")
    tot1[3:end] /= 2 # for c=0.5
    return tot1, r
end

function bee_e_cube(data, db1, db2, numlam; dmean=10.0, abstol=0.03, maxevals=0)
    # FIXME This is only valid for c=0.5
    # Get the first moments using histories and adaptive integration by
    # calling Cubature.hquadrature_v with a function that evaluates the error
    # at multiple points at a time
    mins, maxs = get_bbox(data, factor=2)
    assert(length(db1) == length(db2))
    global iters = 0
    global evals = 0
    function error_1st(points, vals)
        #points: dxn array to evaluate at
        #vals: 2xn values to store into
        numpts = size(points, 2)
        iters += 1
        evals += numpts
        g1 = calc_g(points', db1, numlam, dmean=dmean)
        g2 = calc_g(points', db2, numlam, dmean=dmean)
        vals[:] = exp(min(g1,g2))
    end
    val, err = hcubature_v(error_1st, mins, maxs, abstol=abstol, maxevals=maxevals)
    println("maxs: $maxs, mins: $mins")
    println("Cubature used $iters iterations, and $evals * numclasses evaluations")
    return val*0.5, err
end

function predict(db0, db1, points; dmean=10.0)
    g0 = calc_g(points, db0, 20, dmean=dmean)
    g1 = calc_g(points, db1, 20, dmean=dmean)
    return vec((g0 .- g1) .< 0) * 1.0
end

function error_points(cls::OBC.BinaryClassifier, points, labels; dmean=10.0)
    return sum(abs(predict(cls.mcmc1.db, cls.mcmc2.db, points, dmean=dmean) - labels))/size(labels,1)
end

