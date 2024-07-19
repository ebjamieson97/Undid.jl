module stage1

using Statistics
using LinearAlgebra
using DataFrames
using DelimitedFiles

export create_init_csv, create_diff_df

"""

Takes in vectors of equal lengths and creates the initializing .csv required for UNDID.

"""
function create_init_csv(names=[], start_times=[], end_times=[], treatment_times=[])
    
    # Check input vectors are the same length
    if length(names) == length(start_times) == length(end_times) == length(treatment_times)
        # Do nothing
    else 
        error("Input vectors for names, start_times, end_times and treatment_times must all be the same length.") 
    end 

    # Convert each element of the vectors to a string
    names = string.(names)
    start_times = string.(start_times)
    end_times = string.(end_times)
    treatment_times = string.(treatment_times)

    # Create the data based on input vectors
    header = ["Silo Name", "Start Time", "End Time", "Treatment Time"] 
    data = [header]
    for i in 1:length(names)
        push!(data, [names[i], start_times[i], end_times[i], treatment_times[i]])
    end

    # Define the file path
    filename = "init.csv"
    filepath = abspath(filename)

    # Save as init.csv
    open(filepath, "w") do file
        for i in 1:length(data)
            println(file, join(data[i][1]), ",", data[i][2], ",", data[i][3], ",", data[i][4])
        end
    end

    println("init.csv saved to")
    
    return filepath

end 


"""
Takes in a filepath to the .csv created with create_init_csv and uses that information to create a
dataframe thats shows all the differences that must be calculated for each group.
"""
function create_diff_df(csv)

    data = readdlm(csv, ',')
    header = data[1, :]
    rows = data[2:end, :]
    df = DataFrame(rows, Symbol.(header))

    
    
    header = ["Silo Name", "gvar", "treat", "diff_years"]
    column_types = [Any, Int, Int, Tuple]
    diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])


    
    for group in df[:, "Silo Name"]
    
        years = df[df[!, "Silo Name"] .== group, "Start Time"][1]:df[df[!, "Silo Name"] .== group, "End Time"][1]
        if df[df[!, "Silo Name"] .== group, "Treatment Time"][1] !== 0
            # These are the years that matter for treated groups
            combinations = [(y1, y2) for y1 in years, y2 in years if y1 > y2 && y1 - y2 >= 1 && df[df[!, "Silo Name"] .== group, "Treatment Time"][1] -1 == y2]
            treat = 1
        else 
            # These are the years that matter for untreated groups
            combinations = [(y1, y2) for y1 in years, y2 in years if y1 > y2 && y1 - y2 >= 1]
            treat = 0
        end 

        for i in 1:length(combinations)
            push!(diff_df, [group, combinations[i][2]+1, treat, combinations[i]])
        end 
   
    end

    df = diff_df

    # Remove any rows for which the difference is calcuable for a treatment group but does exist for a control group
    # and vice versa
    unpaired_years = union(setdiff(df[df[!, "treat"] .== 1, "diff_years"], df[df[!, "treat"] .== 0, "diff_years"]),setdiff(df[df[!, "treat"] .== 0, "diff_years"], df[df[!, "treat"] .== 1, "diff_years"]))
    filter!(row -> !(row.diff_years in unpaired_years), df)
    
    df.diff_estimate = Vector{Union{Missing, Float64}}(missing, nrow(df))
    df.diff_se = Vector{Union{Missing, Float64}}(missing, nrow(df))
    
    return df

end 

end 
