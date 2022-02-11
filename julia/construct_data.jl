struct RiskCov19Data
    df::DataFrame
    timecol::Symbol
    deadcol::Symbol
    lincovars::Vector{Symbol}
    odecovars::Vector{Symbol}
end

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

function load_AHS_data(lincovars, odecovars)
    data_filename = joinpath("..", "..", "data/AHS/Restricted/analysis.csv")
    df = CSVread(data_filename, DataFrame)
    # handeling data types
    df = identity.(df) # auto identify the data type
    features = ["Days to admission", "Disch_days", "Dead", 
            "Sex", "Age", "Age group", "Year Month",
            "Myocardial infarction",
            "Congestive Heart Failure", "Peripheral Vascular Disease", "Cerebrovascular Disease", 
            "Dementia", "Chronic Pulmonary Disease", "Rheumatic Disease", "Peptic Ulcer Disease", 
            "Liver disease – mild", "Diabetes without complications", "Diabetes with complications", 
            "Paraplegia and Hemiplegia", "Renal Disease", "Cancer", "Metastatic Carcinoma", 
            "Liver disease – moderate/severe",
            "Num. comorbidities", "Num. symptoms"]
    select!(df, features)
    # Remove missing
    dropmissing!(df)
    # Transform
    rename!(df, replace.(names(df), " " => "_"))
    rename!(df, replace.(names(df), "." => ""))
    transform!(df, :Disch_days => ByRow(x -> Float64(-x)) => :Time,
                :Dead => ByRow(x -> x == true ? 1.0 : 0.0) => :Dead,
                renamecols=true)
    onehot!(df, [:Sex])
    # data type
    return RiskCov19Data(df, :Time, :Dead, lincovars, odecovars)
end

function data_partition(data::RiskCov19Data, args...; kw...)
    dftrain, dftest = partition(data.df, args...; kw...)
    train = RiskCov19Data(dftrain, data.timecol, data.deadcol, data.lincovars, data.odecovars)
    test = RiskCov19Data(dftest, data.timecol, data.deadcol, data.lincovars, data.odecovars)
    return train, test
end