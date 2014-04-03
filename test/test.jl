using Distributions
using OBC
import OBC: propose, energy, reject
using Base.Test

########### Test Sampler ########
type TestParams 
    mu :: Float64
    sigma :: Float64
end
TestParams() = TestParams(0.0,1.0)

function copy!(to::TestParams, from::TestParams)
    to.mu = from.mu
    to.sigma = from.sigma
end

type TestSampler <: Sampler
    curr :: TestParams
    old :: TestParams
    data :: Vector{Float64}
end
TestSampler() = TestSampler(TestParams(), TestParams(), [0.1, 0.1, 0.1]*3)

function propose(s::TestSampler)
    copy!(s.old, s.curr)
    s.curr.mu += randn()*0.2
    s.curr.sigma += randn()*0.2
    s.curr.sigma = clamp(s.curr.sigma, 0.1, Inf)
    0
end

function energy(s::TestSampler)
    -sum(logpdf(Normal(s.curr.mu, s.curr.sigma), s.data)) #Energies are -logp
end

function reject(s::TestSampler)
    copy!(s.curr, s.old)
end
########### /Test Sampler ########

s = TestSampler()
#x = set_energy_limits(s)
mymh = MHRecord(s, verbose=true)
sample(mymh, int(2e4))

@test length(mymh.db) == 20000
@test eltype(mymh.db) == Any
@test typeof(mymh.db[1]) == TestParams

