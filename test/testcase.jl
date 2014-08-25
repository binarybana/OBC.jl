reload("OBC")

data = rand(0:1000,48,18000)

normfac = vec(mapslices(x->quantile(x,0.75), data, 2))

cls = MPM.mpm_classifier(rand(0:100,4,3), rand(0:100,4,3); burn=50, thin=50, d1=rand(4), d2=rand(4), usepriors=true)
MPM.sample(cls,1)

###############################
##### Copy of the above #######
###############################

reload("OBC")

data = rand(0:1000,48,18000)

normfac = vec(mapslices(x->quantile(x,0.75), data, 2))

cls = MPM.mpm_classifier(rand(0:100,4,3), rand(0:100,4,3); burn=50, thin=50, d1=rand(4), d2=rand(4), usepriors=true)
MPM.sample(cls,1)

###############################
# Another copy. Probably unnecessary, 
# but I've sometimes had one or two reps work
# before before segfaulting
###############################

reload("OBC")

data = rand(0:1000,48,18000)

normfac = vec(mapslices(x->quantile(x,0.75), data, 2))

cls = MPM.mpm_classifier(rand(0:100,4,3), rand(0:100,4,3); burn=50, thin=50, d1=rand(4), d2=rand(4), usepriors=true)
MPM.sample(cls,1)
