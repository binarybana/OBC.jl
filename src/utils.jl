
function logsum(x::AbstractArray)
    if length(x) == 1
        return x[1]
    else
        a = maximum(x)
        return a + log(sum_kbn(exp(x-a)))
    end
end

function logsum(x::AbstractArray, w::AbstractArray)
    if length(x) == 1 && w[1] == 1.0
        return x[1]
    else
        a = maximum(x)
        return a + log(sum_kbn(w.*exp(x-a)))
    end
end

function logsum(x::Float64, y::Float64)
    if x > y
        return x+log1p(exp(y-x))
    else
        return y+log1p(exp(x-y))
    end
end

