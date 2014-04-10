using HDF5, JLD

assert(length(ARGS) == 1)
assert(isfile(ARGS[1]))

basename = split(ARGS[1],".")[1]
println("Converting $(ARGS[1]) into $basename.jld")

isfile("$basename.jld") && rm("$basename.jld")
outfid = jldopen("$basename.jld", "w")

dat = readcsv(ARGS[1])
header = convert(Vector{UTF8String}, vec(dat[1,2:end]))
sampleids = convert(Vector{Int}, vec(dat[2:end, 1]))
dat = convert(Matrix{Int}, dat[2:end, 2:end])

write(outfid, "geneids", header)
write(outfid, "sampleids", sampleids)
write(outfid, "data", dat)

close(outfid)
