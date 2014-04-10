
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

function get_grid(data; positive=true, x=0.0, y=0.0)
    @assert size(data,2) == 2
    lx,ly = mapslices(minimum, data, 1)
    hx,hy = mapslices(maximum, data, 1)
    xspread = hx-lx
    yspread = hy-ly

    if positive
        lx = clamp(lx-xspread*0.2, 0.1, Inf)
        hx = hx+xspread*0.2
        ly = clamp(ly-yspread*0.2, 0.1, Inf)
        hy = hy+yspread*0.2
    else
        lx = lx-xspread*0.2
        hx = hx+xspread*0.2
        ly = ly-yspread*0.2
        hy = hy+yspread*0.2
    end
    if x > 0.0
        gext = [0.0, x, 0.0, y]
    else
        gext = [lx, hx, ly, hy]
    end
    gen_grid(map(iround,gext))
end

function gen_grid(gext, N=30)
    lx,hx,ly,hy = gext
    xstride = (hx-lx)/N
    ystride = (hy-ly)/N
    grid = hcat([[i,j] for i in lx:xstride:hx, j in ly:ystride:hy]...)'
    return length(lx:xstride:hx), length(ly:ystride:hy), xstride, ystride, grid
end

