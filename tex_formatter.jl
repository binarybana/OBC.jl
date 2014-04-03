
gnorms = mean(luad_norm, 2)
Tdesc = ["5th quantile", "Median", "95th quantile", "IQR", "Variance"]
order = sortperm(gnorms[genes])
genes = genes[order]
pvals = pvals[order,:]
true_ts = true_ts[order,:]
tvals = tvals[order,:,:]

ints = mapslices(x->quantile(x, [0.05,0.95]), tvals, 3)

for i=1:length(Ts)
   println("") 
   println("""\\begin{table*}[h!]
\\label{tab:pvals}
\\caption{Posterior predictive model diagnostic -- $(Tdesc[i]).}
      \\begin{tabular}{r||c|ccc}
        \\hline
        Gene ID & Mean expression (counts) & \$T(S_n)\$ & 95th int. \$T(x^{rep})\$ & P-value \\\\ 
        \\hline""")
   for j=1:rounds*D
       str = join([
        luad_norm_geneids[genes[j]], 
        (@sprintf "%.1f" gnorms[genes[j]]),
        (@sprintf "%.2f" true_ts[j,i]),
        (@sprintf "[%.2f, %.2f]" ints[j,i,1] ints[j,i,2]),
        (@sprintf "%.2f" pvals[j,i])], " & ") * " \\\\"
      println(str)
    end
    println("""
\\end{tabular}
\\end{table*}
    """)
end


        
       
       
    
