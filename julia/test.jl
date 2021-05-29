using Plots, LaTeXStrings
using DelimitedFiles
using Revise
using JLD2, FileIO

cd("julia/")
push!(LOAD_PATH, ".")
using VirusLoadCurve

tdata = [4.0/24.0, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.]
RAWDATA = readdlm("../../data/Smith2018/Virus_Best10.txt", ' ', Float64, '\n')

k = 6
ttdata = repeat(tdata, inner=1)
tend = ttdata[end]
vvdata = reshape(RAWDATA[k, :], 1, :)[:]

data = VirusLoadCurve.VirusLoadData(ttdata, vvdata)

par0 = [0.70, 2.88, 6.00, 7.60, 0.20, 5.0]
resultVLF = VirusLoadCurve.fitVLF(data; ϵ=0.1)
summary(resultVLF)
save("test.jld2", Dict("result" => resultVLF, "data" => data))

resultVLF = load("test.jld2", "result")

plot(data)
plot!(resultVLF, empirical=true, ylims=(-2, 9))
ylims!((-2, 9))

pt = VirusLoadCurve.Boxplots(resultVLF)
