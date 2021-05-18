using Plots, LaTeXStrings
using DelimitedFiles
using Revise

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

plot(data)
plot!(resultVLF, empirical=true)
