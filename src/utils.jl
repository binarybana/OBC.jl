
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

function get_grid(data; positive=true)
    @assert size(data,2) == 2
    lx,ly = minimum(data, 1)
    hx,hy = maximum(data, 1)
    xspread = max(hx-lx, 1.2)
    yspread = max(hy-ly, 1.2)

    if positive
        lx = clamp(lx-xspread*0.2, 0.0, Inf)
        ly = clamp(ly-yspread*0.2, 0.0, Inf)
    else
        lx = lx-xspread*0.2
        ly = ly-yspread*0.2
    end
    hx = hx+xspread*0.2
    hy = hy+yspread*0.2
    gext = [lx, hx, ly, hy]
    gext, gen_grid(gext)
end

function gen_grid(gext, N=30)
    lx,hx,ly,hy = gext
    xstride = float(hx-lx)/N
    ystride = float(hy-ly)/N
    grid = hcat([[i,j] for i in lx:xstride:hx, j in ly:ystride:hy]...)'
    return length(lx:xstride:hx), length(ly:ystride:hy), xstride, ystride, grid
end

