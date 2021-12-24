using DelimitedFiles, DataFrames, CategoricalArrays
using StatsBase

## Data
# load data
data_filename = pwd()*"/../../data/AHS/Restricted/analysis.csv"
data, header = readdlm(data_filename, ',', header=true)
df = DataFrame(data, vec(header))
# handeling data types
df = identity.(df) # auto identify the data type
df."Age group" = categorical(df."Age group", ordered=true)
levels!(df."Age group", ["Under 1 year", "1-4 years", "5-9 years", "10-19 years", "20-29 years", 
        "30-39 years", "40-49 years", "50-59 years", "60-69 years", "70-79 years", "80+ years"])
# data 
df2 = df[sample(1:size(df, 1), 3000, replace=false), ["Disch_days", "Dead", "Age group"]]
rename!(df2, :Disch_days => :time, :Dead => :dead, Symbol("Age group") => :group)#, Symbol("Age group") => :age_group)
# prepare data
function dataGroups(x)
    if x <= CategoricalValue("30-39 years", df2.group) 
        return 1
    elseif x <= CategoricalValue("60-69 years", df2.group)
        return 2
    else
        return 3
    end
end
groups_info = [
    "age < 40"
    "40 ≤ age < 70"
    "70 ≤ age"
]
transform!(df2, :time => ByRow(x -> Float64(-x)) => :time, 
                :dead => ByRow(x -> x == "True" ? 1.0 : 0.0) => :dead,
                :group => ByRow(dataGroups) => :group)
writedlm("NTCM_data.csv", Iterators.flatten(([names(df2)], eachrow(df2))), ',')
writedlm("NTCM_data_groups_info.csv", eachrow(groups_info), ',')