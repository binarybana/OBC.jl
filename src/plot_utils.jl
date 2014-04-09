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

