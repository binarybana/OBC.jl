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

function sample!(cls::OBC.BinaryClassifier, iters=10000; verbose=false)
    t0 = time()
    OBC.sample!(cls.mcmc1, iters, verbose=verbose)
    t1 = time()
    OBC.sample!(cls.mcmc2, iters, verbose=verbose)
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

############################################################
# BEE Moment calculations
############################################################

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

        g1 = calc_g(points, db1, numlam, dmean=dmean1)
        g2 = calc_g(points, db2, numlam, dmean=dmean2)
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
        temp = dbpass/numpts/2 #FIXME c!=0.5 AND CHECK 4 constant
        beem1 += temp
        beem2 += temp^2
    end
    beem1 /= dblen
    beem2 /= dblen
    return beem1, beem2
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
    g0 = calc_g(points, cls.mcmc1.db, 20, dmean=dmean1)
    g1 = calc_g(points, cls.mcmc2.db, 20, dmean=dmean2)
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

    @assert length(cls.mcmc1.db) == length(cls.mcmc2.db)
    dbsize = length(cls.mcmc1.db)
    bee = 0.0
    bee2 = 0.0
    numpts = div(numpts,2)
    for i=1:dbsize

        pts1 = gen_points(numpts*2, dmean1, cls.mcmc1.db[i]) # x and z
        pts2 = gen_points(numpts*2, dmean2, cls.mcmc2.db[i]) # x and z

        allpts = hcat(pts1,pts2) 

        g1 = clamp(calc_g(allpts, [cls.mcmc1.db[i]], 20, dmean=dmean1), -1e6, Inf)
        g2 = clamp(calc_g(allpts, [cls.mcmc2.db[i]], 20, dmean=dmean2), -1e6, Inf)

        g1eff = clamp(calc_g(allpts, cls.mcmc1.db, 20, dmean=dmean1), -1e6, Inf)
        g2eff = clamp(calc_g(allpts, cls.mcmc2.db, 20, dmean=dmean2), -1e6, Inf)

        # Could inline and optimize calc_g here
        
        g1x = g1[[1:numpts,numpts*2+1:numpts*3]]
        g1z = g1[[numpts+1:numpts*2,numpts*3+1:end]]
        g2x = g2[[1:numpts,numpts*2+1:numpts*3]]
        g2z = g2[[numpts+1:numpts*2,numpts*3+1:end]]

        g1xeff = g1eff[[1:numpts,numpts*2+1:numpts*3]]
        g1zeff = g1eff[[numpts+1:numpts*2,numpts*3+1:end]]
        g2xeff = g2eff[[1:numpts,numpts*2+1:numpts*3]]
        g2zeff = g2eff[[numpts+1:numpts*2,numpts*3+1:end]]

        #agree = (g1x .< g2x) .== (g1z .< g2z)
        labelx = (g1xeff .< g2xeff)
        labelz = (g1zeff .< g2zeff)

        zx = Float64[logsum(g1x[j],g2x[j]) for j=1:numpts*2]
        zz = Float64[logsum(g1z[j],g2z[j]) for j=1:numpts*2]

        lpyx = min(g1x,g2x) .- zx
        lpyz = min(g1z,g2z) .- zz

        both1 = g1x .+ g1z .- zx .- zz
        both2 = g2x .+ g2z .- zx .- zz

        cross = min(both1,both2)

        ## BEE calculation
        temp = exp(logsum(vcat(lpyx,lpyz)))
        #@show exp(minimum(lpyx))
        #@show exp(minimum(lpyz))
        @show hist(exp(lpyx))
        @show extrema(exp(lpyx))
        @show temp/numpts/2
        #@show sum(exp(lpyx))/numpts/2
        @show size(lpyx)
        @show numpts
        #@show exp(lpyx[1:3])
        #@show exp(lpyz[1:3])
        bee += temp

        for j=1:numpts
            if labelx[j]==labelz[j]
                bee2 += exp(cross[j])
            end
        end
    end
    bee /= (numpts*4 * dbsize)
    bee2 /= (numpts*2 * dbsize) # FIXME not sure about this 2 constant

    return bee, bee2, bee2 - bee^2
end

function predict(db0, db1, points; dmean=10.0)
    g0 = calc_g(points, db0, 20, dmean=dmean)
    g1 = calc_g(points, db1, 20, dmean=dmean)
    return vec((g0 .- g1) .< 0) * 1.0
end

function error_points(cls::OBC.BinaryClassifier, points, labels; dmean=10.0)
    return sum(abs(predict(cls.mcmc1.db, cls.mcmc2.db, points, dmean=dmean) - labels))/size(labels,1)
end

