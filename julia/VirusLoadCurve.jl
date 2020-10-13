module VirusLoadCurve

    export fitVLF,
           LogVirusLoadFunction,
           fitVTM,
           LogViralTargetModel,
           VirusLoadData

    using OrdinaryDiffEq, ParameterizedFunctions
    using Optim, LsqFit
    using Plots, LaTeXStrings
    using DelimitedFiles, DataFrames
    using Statistics
    using Printf

    import Base.summary

    struct VirusLoadData{T<:Real}
        t::Vector{T}
        v::Vector{T}
    end

    @recipe function f(data::VirusLoadCurve.VirusLoadData)
        tmin, tmax = extrema(data.t)
        vmin, vmax = extrema(data.v)
        x := data.t
        y := data.v
        seriestype := :scatter
        xaxis --> ("Time (days)", (tmin, tmax))
        yaxis --> (L"\log\,V(t)", (vmin-0.5, vmax+0.5))
        legend -> :none
        grid --> :none
        ()
    end

    include("virus_load_function.jl")
    include("virus_taget_model.jl")

end #module