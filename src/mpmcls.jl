# Everything dealing with TWO classifier objects

function mpm_classifier(data1, data2; burn=1000, thin=50, d=100.0, usepriors=true)
    assert(size(data1,2) == size(data1,2))
    kmax = 1
    D = size(data1, 2)
    cov = eye(D,D) .* 0.1
    mu = ones(D,kmax)
    prior = MPM.MPMPrior(D=D, kmax=kmax, kappa=10., S=eye(D)*10)
    pmoves = MPM.MPMPropMoves()
    pmoves.lammove = 0.005
    pmoves.mumove = 0.6
    pmoves.priorkappa = 100.0

    # Class 1
    start = MPM.MPMParams(mu, cov, ones(kmax)/kmax, #w 
        clamp(log(data1'/d),-8.0,Inf), #lam 
        1) #k :: Int
    obj_a = MPM.MPMSampler(prior, data1, deepcopy(start), pmoves, d)
    obj_a.usepriors = usepriors
    mymh_a = MPM.MHRecord(obj_a,burn=burn,thin=thin)

    # Class 2
    start.lam = clamp(log(data2'/d),-8.0,Inf) # lam
    obj_b = MPM.MPMSampler(prior, data2, deepcopy(start), pmoves, d)
    obj_b.usepriors = usepriors
    mymh_b = MPM.MHRecord(obj_b,burn=burn,thin=thin)

    OBC.BinaryClassifier(obj_a, obj_b, mymh_a, mymh_b)
end

function sample(cls::OBC.BinaryClassifier, iters=10000)
    t0 = time()
    OBC.sample(cls.mcmc1, iters)
    t1 = time()
    OBC.sample(cls.mcmc2, iters)
    t2 = time()
    return t1-t0, t2-t1
end

function accuracies(cls::OBC.BinaryClassifier)
    acc1 = cls.mcmc1.count_accept / cls.mcmc1.count_total
    acc2 = cls.mcmc2.count_accept / cls.mcmc2.count_total
    acc1,acc2
end

function bee_moments(points,db1::Vector{Any},db2::Vector{Any},numlam,volume;dmean=10.0)
    # FIXME This is only valid for c=0.5
    g1 = calc_g(points, db1, numlam)
    g2 = calc_g(points, db2, numlam)
    bee_moments(points,db1,db2,g1,g2,numlam,volume,dmean=dmean)
end

function bee_moments(points,db1::Vector{Any}, db2::Vector{Any}, 
        g1::Vector{Float64},g2::Vector{Float64},numlam,volume;dmean=10.0)
    # FIXME This is only valid for c=0.5
    bee = exp(logsum(min(g1,g2) .+ log(volume)))*0.5
    numpts = size(points, 1)
    dims = size(points,2)
    dblen = length(db1)
    assert(length(db1) == length(db2))

    lams1 = Array(Float64, dims, numlam)
    lams2 = Array(Float64, dims, numlam)
    llams1 = Array(Float64, dims, numlam)
    llams2 = Array(Float64, dims, numlam)
    rlams1 = Array(Float64, dims, numlam)
    rlams2 = Array(Float64, dims, numlam)
    bee2 = dbpass = 0.0
    points = points'
    lgpoints = lgamma(points .+ 1)

    for i in 1:dblen 
        dbpass = 0.0
        curr1 = db1[i]
        curr2 = db2[i]
        assert(curr1.k==1 && curr2.k==1)
        rand!(MultivariateNormal(curr1.mu[:,1], curr1.sigma), lams1)
        rand!(MultivariateNormal(curr2.mu[:,1], curr2.sigma), lams2)
        for k=1:numlam*dims 
            rlams1[k] = dmean * exp(lams1[k])
            llams1[k] = log(dmean) + lams1[k]
            rlams2[k] = dmean * exp(lams2[k])
            llams2[k] = log(dmean) + lams2[k]
        end
        for j in 1:numpts # each point
            accumlam1 = 0.0
            accumlam2 = 0.0
            for s in 1:numlam # each lambda
                accumD1 = 0.0
                accumD2 = 0.0
                for d in 1:dims # each dimension
                    accumD1 += points[d,j]*llams1[d,s] - rlams1[d,s] - lgpoints[d,j]
                    accumD2 += points[d,j]*llams2[d,s] - rlams2[d,s] - lgpoints[d,j]
                end
                accumlam1 += exp(accumD1)
                accumlam2 += exp(accumD2)
            end
            temp = g1[j] < g2[j] ? accumlam1 : accumlam2
            #temp = accumlam1 # for testing purposes
            dbpass += temp
        end
        temp = 0.5*dbpass/numlam*volume
        bee2 += temp^2
    end
    bee2 /= dblen
    return bee, bee2
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

function predict(db0, db1, points)
    g0 = calc_g(points, db0, 20)
    g1 = calc_g(points, db1, 20)
    return vec((g0 .- g1) .< 0) * 1
end

function error_points(db0, db1, points, labels)
    return sum(abs(predict(db0, db1, points) - labels))/size(labels,1)
end

