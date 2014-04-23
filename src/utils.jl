function logsum(x::AbstractArray)
    if length(x) == 1
        return x[1]
    else
        a = maximum(x)
        return a .+ log(sum_kbn(exp(x.-a)))
    end
end

function logsum(x::AbstractArray, w::AbstractArray)
    if length(x) == 1 && w[1] == 1.0
        return x[1]
    else
        a = maximum(x)
        return a .+ log(sum_kbn(w.*exp(x.-a)))
    end
end

function logsum(x::Float64, y::Float64)
    x > y ? x+log1p(exp(y-x)) : y+log1p(exp(x-y))
end

function get_bbox(data; factor=1.5, positive=true)
    D = size(data,2)

    #maxs = vec(maximum(data, 1))
    maxs = vec(mapslices(x->quantile(x,0.75), data, 1))
    maxs .+= 3*maximum(maxs) / 4

    if positive
        mins = zeros(D)
    else
        #mins = vec(minimum(data,1))
        maxs = vec(mapslices(x->quantile(x,0.25), data, 1))
    end

    spreads = maxs .- mins

    if !positive
        mins .-= spreads .* factor
    end
    maxs .+= spreads .* factor
    mins, maxs
end

function get_grid(data; factor=1.5, N=30, positive=true)
    mins, maxs = get_bbox(data, factor=factor, positive=positive)
    mins, maxs, gen_grid(mins, maxs, N)
end

gen_unit_grid(mins,maxs) = gen_grid(mins,maxs,typemax(Int))

function gen_grid(mins, maxs, N=30)
    D = length(mins)
    stepsizes = ceil(float(maxs .- mins) ./N)
    ranges = [mins[i]:stepsizes[i]:maxs[i] for i=1:D]
    grid = hcat(map(collect, collect(product(ranges...)))...)'
    return map(length,ranges), stepsizes, grid
end

