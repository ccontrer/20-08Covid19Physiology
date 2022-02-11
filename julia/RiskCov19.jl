module RiskCov19

using DataFrames: DataFrame, select!, dropmissing!, rename!, transform!, ByRow, Not
using OrdinaryDiffEq: ODEProblem, ODESolution, remake, Tsit5, solve
using CSV: read as CSVread
using Turing
using LinearAlgebra: I, diagm
using QuadGK: quadgk
using MLJ: partition
# using Reexport

# @reexport begin
#     using MLJ: partition
# end

include("construct_data.jl")
include("ode_models.jl")
include("turing_model.jl")

export RiskCov19Data, load_AHS_data
export model001
export train, turing_model
export data_partition

end #module