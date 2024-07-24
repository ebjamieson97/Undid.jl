module UNDID

using Statistics
using LinearAlgebra
using DataFrames
using DelimitedFiles

export create_init_csv, create_diff_df, fill_diff_df, read_diff_df, combine_silo_data, create_ATT_by_gt_by_silo, compute_ATT


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

function create_diff_df(csv)

    data = readdlm(csv, ',')
    header = data[1, :]
    rows = data[2:end, :]
    df = DataFrame(rows, Symbol.(header))

    
    
    header = ["Silo Name", "gvar", "treat", "diff_years", "(g,t)"]
    column_types = [Any, Int, Int, Tuple, Tuple]
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
            push!(diff_df, [group, combinations[i][2]+1, treat, combinations[i], (combinations[i][2]+1, combinations[i][1])])
        end 
   
    end

    df = diff_df

    # Remove any rows for which the difference is calcuable for a treatment group but does exist for a control group
    # and vice versa
    unpaired_years = union(setdiff(df[df[!, "treat"] .== 1, "diff_years"], df[df[!, "treat"] .== 0, "diff_years"]),setdiff(df[df[!, "treat"] .== 0, "diff_years"], df[df[!, "treat"] .== 1, "diff_years"]))
    filter!(row -> !(row.diff_years in unpaired_years), df)
    
    df.diff_estimate = Vector{Union{Missing, Float64}}(missing, nrow(df))
    df.diff_se = Vector{Union{Missing, Float64}}(missing, nrow(df))


    # Define the file path
    filename = "empty_diff_df.csv"
    filepath = abspath(filename)

    # Save as empty_diff_df.csv
    open(filepath, "w") do file
        println(file, join("Silo Name"), ";", "gvar", ";", "treat", ";", "diff_years", ";", "(g,t)", ";", "diff_estimate", ";",
                             "diff_se")
        for i in 1:nrow(df)
            println(file, join(df[!, "Silo Name"][i]), ";", df[!, "gvar"][i], ";", df[!, "treat"][i], ";", df[!, "diff_years"][i],
                                ";", df[!, "(g,t)"][i], ";", df[!, "diff_estimate"][i], ";", df[!, "diff_se"][i])
        end
    end

    println("empty_diff_df.csv saved to")
    formatted_path = replace(filepath, "\\" => "\\\\")
    println(formatted_path)

    
    return df
    

end 

function fill_diff_df(silo_name, empty_diff_df, local_data)

    empty_diff_df = empty_diff_df[empty_diff_df[!, "Silo Name"] .== silo_name,:]

    for diff_combo in empty_diff_df[empty_diff_df[!, "Silo Name"] .== silo_name,:].diff_years
        silo_subset = local_data[(local_data[!, "year"] .== diff_combo[1] .|| local_data[!, "year"] .== diff_combo[2]),:]
        silo_subset.year = replace(silo_subset.year, diff_combo[1] => 1, diff_combo[2] => 0)
        X = hcat(1 .- silo_subset.year, silo_subset.year)
        Y = silo_subset.coll
        beta_hat = X\Y
        gamma = beta_hat[2] - beta_hat[1]
        empty_diff_df[(empty_diff_df[!, "Silo Name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_years"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_years"])), "diff_estimate"] .= gamma
    end 

     # Define the file path
     filename = "filled_diff_df_$(silo_name).csv"
     filepath = abspath(filename)
 
     # Save as empty_diff_df.csv
     open(filepath, "w") do file
         println(file, join("Silo Name"), ";", "gvar", ";", "treat", ";", "diff_years", ";", "(g,t)", ";", "diff_estimate", ";",
                              "diff_se")
         for i in 1:nrow(empty_diff_df)
             println(file, join(empty_diff_df[!, "Silo Name"][i]), ";", empty_diff_df[!, "gvar"][i], ";", empty_diff_df[!, "treat"][i], ";", empty_diff_df[!, "diff_years"][i],
                                 ";", empty_diff_df[!, "(g,t)"][i], ";", empty_diff_df[!, "diff_estimate"][i], ";", empty_diff_df[!, "diff_se"][i])
         end
     end
 
     println("filled_diff_df_$(silo_name).csv saved to")
     formatted_path = replace(filepath, "\\" => "\\\\")
     println(formatted_path)

    return empty_diff_df

end 

function read_diff_df(filepath_to_diff_df)
    data = readdlm(filepath_to_diff_df, ';')
    header = data[1, :]
    rows = data[2:end, :]
    df_readin = DataFrame(rows, Symbol.(header))

    transform!(df_readin, :diff_years => ByRow(s -> tuple(parse.(Float64, split(strip(s, ['(', ')']), ','))...)) => :diff_years)
    transform!(df_readin, :"(g,t)" => ByRow(s -> tuple(parse.(Float64, split(strip(s, ['(', ')']), ','))...)) => :"(g,t)")

    return df_readin
end 

function combine_silo_data(dir_path)
    files = readdir(dir_path)
    matching_files = filter(file -> startswith(file, "filled_diff_df_") && endswith(file, ".csv"), files)

    data = read_diff_df("$dir_path\\$(matching_files[1])")
    for i in 2:length(matching_files)
        data = vcat(data, read_diff_df("$dir_path\\$(matching_files[i])"))
    end 

    return data
end

function create_ATT_by_gt_by_silo(combined_diff_data; save_as_csv=false)
    treated_silo = []
    control_silo = []
    diff_years = []
    gt_vec = []
    ATT = []

    
    for treated_group in unique(combined_diff_data[combined_diff_data[!,"treat"] .== 1,"Silo Name"])

        for time in unique(combined_diff_data[combined_diff_data[!, "Silo Name"] .== treated_group,"diff_years"])

            for control_group in combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_years"])) .&& (in.(time[2], combined_diff_data[!, "diff_years"])) .&& (combined_diff_data[!, "treat"] .== 0), "Silo Name"]
            
                control_diff = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_years"])) .&& (in.(time[2], combined_diff_data[!, "diff_years"])) .&& (combined_diff_data[!, "treat"] .== 0) .&& (combined_diff_data[!, "Silo Name"] .== control_group), "diff_estimate"]
                treated_diff = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_years"])) .&& (in.(time[2], combined_diff_data[!, "diff_years"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "Silo Name"] .== treated_group), "diff_estimate"]

                attgt = treated_diff - control_diff

                push!(treated_silo, treated_group)
                push!(control_silo, control_group)
                push!(diff_years, time)
                push!(gt_vec, combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_years"])) .&& (in.(time[2], combined_diff_data[!, "diff_years"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "Silo Name"] .== treated_group), "(g,t)"][1])
                push!(ATT, attgt[1])

            end 

        end 

    end 

    ATT_by_gt_by_silo_df = ATT_by_gt_by_silo_df = DataFrame("TreatedSilo" => treated_silo, "ControlSilo" => control_silo,
                                                            "diff_years" => diff_years, "ATT" => ATT, "(g,t)" => gt_vec)
    if save_as_csv == true
        
        filename = "ATT_by_gt_by_silo.csv"
        filepath = abspath(filename)

        # Save as empty_diff_df.csv
        open(filepath, "w") do file
            println(file, join("TreatedSilo"), ";", "ControlSilo", ";", "diff_years", ";", "ATT", ";", "(g,t)")
            for i in 1:nrow(ATT_by_gt_by_silo_df)
                println(file, join(ATT_by_gt_by_silo_df[!, "TreatedSilo"][i]), ";", ATT_by_gt_by_silo_df[!, "ControlSilo"][i], ";", ATT_by_gt_by_silo_df[!, "diff_years"][i], ";", ATT_by_gt_by_silo_df[!, "ATT"][i],
                                ";", ATT_by_gt_by_silo_df[!, "(g,t)"][i])
            end
        end

        println("ATT_by_gt_by_silo.csv saved to")
        formatted_path = replace(filepath, "\\" => "\\\\")
        println(formatted_path)
    end


    return ATT_by_gt_by_silo_df

end

function compute_ATT(ATT_by_gt_by_silo)

    ATTs_bysilo = []
    treated_group = []
    avg_ATT_across_silos = Vector{Any}(fill(nothing, length(unique(ATT_by_gt_by_silo[!,"TreatedSilo"]))))

    for group in unique(ATT_by_gt_by_silo[!,"TreatedSilo"])
        mean_att = mean(ATT_by_gt_by_silo[ATT_by_gt_by_silo[!, "TreatedSilo"] .== group, "ATT"])
    
        push!(ATTs_bysilo, mean_att)
        push!(treated_group, group)
    end 

    avg_ATT_across_silos[1] = mean(ATTs_bysilo)

    ATTs_bysilo_df = DataFrame(Group = treated_group, ATT_bysilo = ATTs_bysilo, avg_ATT_across_silos = avg_ATT_across_silos)

    return ATTs_bysilo_df
end 




end 
