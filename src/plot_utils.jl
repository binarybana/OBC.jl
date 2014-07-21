using PyPlot

function plot_data(cls)
    figure()
    #mins, maxs, (lens, steps, grid) = MPM.get_grid(vcat(data,datab))

    #close("all")

    #ga = MPM.calc_g(grid, mymh_a.db, 20)
    #ga = exp(ga.-maximum(ga))
    #imshow(reshape(ga,lens...)', extent=[mins[1],maxs[1],mins[2],maxs[2]], aspect="equal", origin="lower")
    #colorbar()
    #contour(reshape(ga,n1,n2)', extent=gext, aspect="equal", origin="lower")
    #plot(data[:,1], data[:,2], "g.", alpha=0.8)

    #allnodes = collect(r)
    #for n in allnodes
        #if length(n.vals)!=0
            #plot(n.mins[1], n.mins[2], "g.")
        #end
    #end

    #figure()
    #gavg = ga - gb
    #imshow(reshape(gavg,n1,n2)', extent=gext, aspect="equal", origin="lower")
    #colorbar()
    #contour(reshape(gavg,n1,n2)', (0.0,), extent=gext, aspect="equal", origin="lower")
    data = cls.cls1.data
    datab = cls.cls2.data
    plot(data[:,1], data[:,2], "g.", alpha=0.8)
    plot(datab[:,1], datab[:,2], "r.", alpha=0.8)
end

function plot_traces(db, names; save=false)
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
        if save
            savefig("$(string(name)).png")
        end
    end
end

