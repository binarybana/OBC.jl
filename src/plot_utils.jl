using PyPlot

function plot_traces(db, names)
    if !issubtype(typeof(names), Array)
        names = [names]
    end
    for name in names
        figure()
        if issubtype(typeof(getfield(db[1],name)), Array)
            vals = [vec(getfield(x,name)) for x in db]
        else
            vals = [getfield(x,name) for x in db]
        end
        plot(vals)
        legend(map(string,1:length(vals[1])))
        ylabel(string(name))
        title(string(name))
        xlabel("samples")
    end
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

function gen_grid(gext)
    N = 30
    lx,hx,ly,hy = gext
    xstride = max(1, iround((hx-lx)/N))
    ystride = max(1, iround((hy-ly)/N))
    grid = hcat([[i,j] for i in lx:xstride:hx, j in ly:ystride:hy]...)'
    return length(lx:xstride:hx), length(ly:ystride:hy), gext, grid
end

