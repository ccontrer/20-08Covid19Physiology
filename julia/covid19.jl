using OrdinaryDiffEq, ParameterizedFunctions
using LsqFit, Optim
using Plots, LaTeXStrings
using DelimitedFiles
using Statistics

tdata = [range(0, 12 , step=1.0);]
tend = tdata[end]
data = readdlm("../data/Smith2018/Virus_Best10.txt", ' ', Float64, '\n')'

plt = plot(legend=:none, grid=:none, xlabel="Time (days)", ylabel=L"\log\,V(t)")
[plot!(plt, tdata, data[:, i], seriestype=:scatter) for i in 1:size(data, 2)]
plt

ttdata = reshape(repeat(tdata', 10), 1, :)[:]
vvdata = reshape(data', 1, :)[:]
plot(ttdata, vvdata, seriestype=:scatter,
    legend=:none, grid=:none, xlabel="Time (days)", ylabel=L"\log\,V(t)")