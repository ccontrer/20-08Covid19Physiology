using DelimitedFiles, DataFrames, CategoricalArrays
using StatsBase

function onehot!(df::DataFrame, col::Union{Symbol,String})
    cats = unique(df[!, col])
    col_names = string.(col, "_", cats)
    col_names = replace.(col_names, " " => "_")
    transform!(df, col => (x -> Float64.(x .== permutedims(cats))) => col_names)
    select!(df, Not(col))
end

function onehot!(df::DataFrame, cols::Vector)
    for col in cols
        onehot!(df, col)
    end
    df
end

onehot(df::DataFrame, col) = onehot!(copy(df), col)

struct RiskCov19Data
    df::DataFrame
    timecol::Symbol
    deadcol::Symbol
    covarcols::Vector{Symbol}
    odegrcols::Vector{Symbol}
end


# ## Data
# # load data
# data_filename = pwd()*"/../../data/AHS/Restricted/analysis.csv"
# data, header = readdlm(data_filename, ',', header=true)
# df = DataFrame(data, vec(header))
# # handeling data types
# df = identity.(df) # auto identify the data type
# df."Age group" = categorical(df."Age group", ordered=true)
# levels!(df."Age group", ["Under 1 year", "1-4 years", "5-9 years", "10-19 years", "20-29 years", 
#         "30-39 years", "40-49 years", "50-59 years", "60-69 years", "70-79 years", "80+ years"])
# transform!(df, :Disch_days => ByRow(x -> Float64(-x)) => :Time, 
#                :Dead => ByRow(x -> x == "True" ? 1.0 : 0.0) => :Dead)
# # Remove empty
# deleteat!(df, df.Sex.=="")
# # features
# features = ["Age group"]
# ode_groups = ["Sex"]
# select!(df, vcat(["Time", "Dead"], features, ode_groups))
# # one hot encoding
# onehot!(df, features)
# transform!(df, :Sex => ByRow(x -> x=="Male" ? 1 : 0) => :Male)
# new_features = names(df)[3:end]

# data = RiskCov19Data(df, :Time, :Dead, )

# writedlm("NTCM_data.csv", Iterators.flatten(([names(df2)], eachrow(df2))), ',')
# writedlm("NTCM_data_groups_info.csv", eachrow(groups_info), ',')