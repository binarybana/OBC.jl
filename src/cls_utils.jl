
function plot_all(mha,mhb;hidden=false)
    colors=["g","r"]
    figure()
    for (i,mh) in enumerate((mha, mhb))
        bestparams = mh.mapvalue.curr
        bestpred = pred(bestparams.theta, bestparams.init1, bestparams.init2)
        for j in 1:size(mh.obj.data, 1)
            if j == 1
                plot(mh.obj.data[j,:]', colors[i], alpha=0.3, linewidth=2, label="Class $i")
            elseif j == 2
                plot(mh.obj.data[j,:]', colors[i], alpha=0.3, linewidth=2, label="Class $i")
            else
                plot(mh.obj.data[j,:]', colors[i], alpha=0.3, linewidth=2)
            end
        end
        if !hidden
            plot(bestpred, colors[i], linewidth=5)
        end
    end
end

function plot_trace(mha, mhb, name)
    colors=["g","r"]
    figure()
    for (i,mh) in enumerate((mha, mhb))
        accum = Float64[]
        for params in mh.db
            push!(accum, params.(symbol(name)))
        end
        plot(accum, colors[i])
        println("Mean: $(mean(accum)), stddev: $(std(accum))")
    end
end

function plot_posterior(mh)
    figure()
    for params in mh.db
        plot(pred(params.theta, params.init1, params.init2), "k")
    end
end
