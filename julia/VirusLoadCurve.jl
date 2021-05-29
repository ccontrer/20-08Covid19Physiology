module VirusLoadCurve

    export fitVLF,
           LogVirusLoadFunction,
           fitVTM,
           LogViralTargetModel,
           VirusLoadData,
           Boxplots
    using OrdinaryDiffEq
    using Optim, LsqFit
    using Plots, LaTeXStrings
    using DelimitedFiles, DataFrames
    using Statistics
    using Turing
    using Printf
    using ProgressMeter

    import Base.summary

    struct VirusLoadData{T<:Real}
        t::Vector{T}
        v::Vector{T}
    end

    @recipe function f(data::VirusLoadData)
        tmin, tmax = extrema(data.t)
        vmin, vmax = extrema(data.v)
        x := data.t
        y := data.v
        seriestype := :scatter
        xaxis := ("Time (days)", (0.0, tmax), font(14))
        yaxis := (L"\log\,V(t)", font(14))
        label := "Data"
        grid := :none
        ()
    end

    include("virus_load_function.jl")
    include("virus_target_model.jl")

end #module