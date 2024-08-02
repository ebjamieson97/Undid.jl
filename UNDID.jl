module UNDID

using Statistics
using LinearAlgebra
using DataFrames
using DelimitedFiles
using Dates

export create_init_csv, create_diff_df, fill_diff_df, read_csv_data, combine_diff_data, create_ATT_by_gt_by_silo, compute_ATT_staggered,
        create_trends_df, run_stage_two, combine_trends_data, compute_ATT_common, compute_covariance_matrix, save_df_as_csv, run_stage_one,
        run_stage_three

### Stage 1 Functions ###
function create_init_csv(names=[], start_times=[], end_times=[], treatment_times=[]; covariates = false)
    
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
    header = ["silo_name", "start_time", "end_time", "treatment_time"] 
    data = [header]
    for i in 1:length(names)
        push!(data, [names[i], start_times[i], end_times[i], treatment_times[i]])
    end

    
    # Handle covariate information
    if covariates != false
        if typeof(covariates) == Vector{String}
            cov_column = fill(covariates, length(end_times))    
            push!(data[1], "covariates")
            for i in 1:length(names)
                push!(data[i+1], join(cov_column[i], ";"))
            end 
        else 
            error("Covariates must be entered as a vector of strings or set to covariates = false to indicate no covariates.")
        end
    end  

    
    # Define the file path
    filename = "init.csv"
    filepath = abspath(filename)

    # Save as init.csv
    open(filepath, "w") do file
        for row in data
            println(file, join(row, ","))
        end
    end

    # Print the location where the .csv was saved to
    println("init.csv saved to")
    formatted_path = replace(filepath, "\\" => "\\\\")
    println(formatted_path)

    return filepath

end 

function create_diff_df(csv; covariates = false)

    # Read in the init.csv and make it a dataframe
    data = readdlm(csv, ',')
    header = data[1, :]
    rows = data[2:end, :]
    df = DataFrame(rows, Symbol.(header))

    # Setup for staggered adoption
    if length(unique(df.treatment_time)) > 2

        # Set up the skeleton of the full empty_diff_df
        header = ["silo_name", "gvar", "treat", "diff_times", "(g,t)"]
        column_types = [Any, Int, Int, Tuple, Tuple]
        diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])


        # Produces the appropriate year combinations that need to be calculated based on the init.csv
        for group in df[:, "silo_name"]
        
            years = df[df[!, "silo_name"] .== group, "start_time"][1]:df[df[!, "silo_name"] .== group, "end_time"][1]
            if df[df[!, "silo_name"] .== group, "treatment_time"][1] !== 0
                # These are the years that matter for treated groups
                combinations = [(y1, y2) for y1 in years, y2 in years if y1 > y2 && y1 - y2 >= 1 && df[df[!, "silo_name"] .== group, "treatment_time"][1] -1 == y2]
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



        # Remove any rows for which the difference is calcuable for a treatment group but does exist for a control group
        # and vice versa
        unpaired_years = union(setdiff(diff_df[diff_df[!, "treat"] .== 1, "diff_times"], diff_df[diff_df[!, "treat"] .== 0, "diff_times"]),setdiff(diff_df[diff_df[!, "treat"] .== 0, "diff_times"], diff_df[diff_df[!, "treat"] .== 1, "diff_times"]))
        filter!(row -> !(row.diff_times in unpaired_years), diff_df)

    # For common treatment time
    elseif length(unique(df.treatment_time)) <= 2

        # Set up the skeleton of the full empty_diff_df
        header = ["silo_name", "treat", "common_treatment_time"]
        column_types = [Any, Int, Int]
        diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])
        common_treatment_time = df[df[!, "treatment_time"] .!== 0, "treatment_time"][1]

        for silo_name in unique(df.silo_name)
            if df[df[!, "silo_name"] .== silo_name, "treatment_time"][1] !== 0
                treat = 1
            else 
                treat = 0
            end

            push!(diff_df, [silo_name, treat, common_treatment_time])

        end 
        
    end

    # Create the empty columns for the diff_esimates and diff_standard_errors
    diff_df.diff_estimate = Vector{Union{Missing, Float64}}(missing, nrow(diff_df))
    diff_df.diff_var = Vector{Union{Missing, Float64}}(missing, nrow(diff_df))
    diff_df.diff_estimate_covariates = Vector{Union{Missing, Float64}}(missing, nrow(diff_df))
    diff_df.diff_var_covariates = Vector{Union{Missing, Float64}}(missing, nrow(diff_df))
    
    # Add the covariates if they exist
    if covariates == false
        if "covariates" in DataFrames.names(df) 
            covariates = [String(s) for s in split(df[!,"covariates"][1], ";")] 
            diff_df.covariates = fill(covariates, nrow(diff_df))
        else
            covariates = "none"
            diff_df.covariates = fill(covariates, nrow(diff_df))
        end
    else
        if typeof(covariates) == Vector{String}
        diff_df.covariates = fill(covariates, nrow(diff_df))
        else
            error("Covariates must be entered as a vector of strings or set to covariates = false to get covariate information from the init.csv.")
        end
    end 

    # Define the file path
    filename = "empty_diff_df.csv"
    filepath = abspath(filename)

    

    # Save as empty_diff_df.csv
    open(filepath, "w") do file
        println(file, join(DataFrames.names(diff_df), ";"))
        for row in eachrow(diff_df)
            println(file, join(row, ";"))
        end
    end

    # Print location where empty_diff_df.csv was saved to
    println("empty_diff_df.csv saved to")
    formatted_path = replace(filepath, "\\" => "\\\\")
    println(formatted_path)

    # Returns the empty_diff_df
    return diff_df
    
end 

function run_stage_one(names, start_times, end_times, treatment_times; covariates = false)
    csv = create_init_csv(names, start_times, end_times, treatment_times)
    return create_diff_df(csv, covariates = covariates)
end

### Stage 2 Functions ### 
function read_csv_data(filepath_to_csv_data)
    
    # Reads in the empty_diff_df.csv
    data = readdlm(filepath_to_csv_data, ';')
    header = data[1, :]
    rows = data[2:end, :]
    df_readin = DataFrame(rows, Symbol.(header))

    # Converts the strings back into tuples for the (g,t) and diff_times columns
    # And converts the string of covariates back into a vector of strings
    if "diff_times" in DataFrames.names(df_readin)
        transform!(df_readin, :diff_times => ByRow(s -> tuple(parse.(Float64, split(strip(s, ['(', ')']), ','))...)) => :diff_times)
    end

    if "(g,t)" in DataFrames.names(df_readin)
        transform!(df_readin, :"(g,t)" => ByRow(s -> tuple(parse.(Float64, split(strip(s, ['(', ')']), ','))...)) => :"(g,t)")
    end

    if "covariates" in DataFrames.names(df_readin)
        transform!(df_readin, :covariates => ByRow(s -> [String(sub) for sub in split(replace(replace(strip(s, ['[' , ']', '"']), "\"" => ""), " " => ""), ",")]) => :covariates)
    end 

    return df_readin
end

function fill_diff_df(silo_name, empty_diff_df, silo_data; renaming_dictionary = false)

    # Select only the relevant silo from the empty_diff_df
    empty_diff_df = empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name,:]

    # Rename columns if necessary
    if renaming_dictionary == false
        # do nothing
    elseif typeof(renaming_dictionary) == Dict{String, String}
        # Check that all of the dictionary values (new column names) are in the empty_diff_df and all dictionary keys are in the local silo data column names and then rename
        if all(value -> value in empty_diff_df.covariates[1], values(renaming_dictionary)) && all(key -> key in DataFrames.names(silo_data), keys(renaming_dictionary))
            rename!(silo_data, renaming_dictionary)
        elseif all(value -> value in empty_diff_df.covariates[1], values(renaming_dictionary))
            error("Ensure that the keys in the renaming_dictionary are the covariate names in the local silo data.")
        elseif all(key -> key in DataFrames.names(silo_data), keys(renaming_dictionary))
            error("Ensure that the values in the renaming_dictionary are the covariate names in the empty_diff_df.")
        else
            error("Ensure that the keys in the renaming_dictionary are the covariate names in the local silo data and that the values of the renaming_dictionary are the corresponding covariate names in the empty_diff_df.")
        end  
    else
        error("If renaming columns ensure that the typeof(renaming_dictionary) == Dict{String, String}")
    end 

    # Compute diff (gamma) estimates and standard errors in the case of common treatment times across silos
    if "common_treatment_time" in DataFrames.names(empty_diff_df)
        
        silo_data_copy = copy(silo_data)

        # Use treatment time to construct dummy varaible indicating obs after treatment time
        treatment_time = empty_diff_df.common_treatment_time[1]
        silo_data_copy.year = ifelse.(silo_data_copy.year .>= treatment_time, 1, 0)

        # Construct X and Y for regression
        X = convert(Matrix{Float64},hcat(fill(1,nrow(silo_data_copy)), silo_data_copy.year))
        Y = convert(Vector{Float64},silo_data_copy.coll)

        # Tell the next for loop to calculate diff_estimate with covariates if cov_on == 2
        if empty_diff_df.covariates[1][1] == "none"
            cov_on = 0
        else 
            cov_on = 1
        end

        # Calculates diff_estimate and SE with and without covariates
        for i in 1:cov_on+1
            X = convert(Matrix{Float64}, X)
            beta_hat = X\Y 
            resid = Y - X*beta_hat
            sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
            cov_beta_hat = compute_covariance_matrix(X, sigma_sq)
            
            if i == 1 # For no covariates
                empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "diff_estimate"] .= beta_hat[2]
                empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "diff_var"] .= cov_beta_hat[2,2]
        
            elseif i == 2 # For covariates
                empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "diff_estimate_covariates"] .= beta_hat[2]
                empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "diff_var_covariates"] .= cov_beta_hat[2,2]
            end 

            if cov_on == 1
                try
                    X = hcat(X, Matrix(silo_data_copy[:, empty_diff_df.covariates[1]]))
                catch ex
                    if isa(ex, ArgumentError)
                        error("Please rename your covariates to align with the names used here: $(empty_diff_df.covariates[1])")
                    else
                        rethrow(ex)
                    end
                end
            end
            
        end


    else

        # Compute the local silo regression for each relevant year
        for diff_combo in empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name,:].diff_times

            # Filters the silo down to the relevant year combinations and replaces the values with 1 (post) and 0 (pre)
            silo_subset = silo_data[(silo_data[!, "year"] .== diff_combo[1] .|| silo_data[!, "year"] .== diff_combo[2]),:]
            silo_subset.year = replace(silo_subset.year, diff_combo[1] => 1, diff_combo[2] => 0)

            # Define X and Y matrix for regression and add covariates if applicable
            X = hcat(fill(1, length(silo_subset.year)), silo_subset.year)
            Y = silo_subset.coll

            # Perform regression and throw gamma and var(gamma) to the empty_diff_df
            beta_hat = X\Y
            resid = Y - X*beta_hat
            sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
            cov_beta_hat = compute_covariance_matrix(X, sigma_sq, diff_times = diff_combo)
            
            empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_times"])), "diff_estimate"] .= beta_hat[2]
            empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_times"])), "diff_var"] .= cov_beta_hat[2,2]


            if (empty_diff_df[!, "covariates"][1][1] == "none") 
                # do nothing
            else
                # Perform regression with covariates
                columns = Symbol.(empty_diff_df[!, "covariates"][1])
                
                try # Try appending columns and throw error message if column names are not correctly specified
                    X = hcat(X, Matrix(silo_subset[:, columns]))
                catch ex
                    if isa(ex, ArgumentError)
                        error("Please rename your covariates to align with the names used here: $(empty_diff_df.covariates[1])")
                    else
                        rethrow(ex)
                    end
                end

                beta_hat = X\Y  
                resid = Y - X*beta_hat
                sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
                cov_beta_hat = compute_covariance_matrix(X, sigma_sq, diff_times = diff_combo, covariates = empty_diff_df[!, "covariates"][1])              
                empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_times"])), "diff_estimate_covariates"] .= beta_hat[2]
                empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_times"])), "diff_var_covariates"] .= cov_beta_hat[2,2]
            end

        end 
    end

    # Save as csv
    silo_name = string(silo_name)
    save_df_as_csv("filled_diff_df_$silo_name.csv", empty_diff_df)
    

    return empty_diff_df

end 

function create_trends_df(silo_data, silo_name; covariates = ["none"], treatment_time = missing)

    silo_name = string(silo_name)

    # Define column headers
    header = ["silo_name", "treatment_time", "time", "mean_outcome", "mean_outcome_residualized", "covariates"]
    data = [header]

    # Push means and time to data
    if covariates == ["none"]
        for x in minimum(silo_data[!,"year"]):maximum(silo_data[!,"year"])
            push!(data, [silo_name, string(treatment_time), string(x), string(mean(silo_data[silo_data[!, "year"] .== x, "coll"])), "n/a", "none"])
        end
    else
        for x in minimum(silo_data[!,"year"]):maximum(silo_data[!,"year"])
            silo_subset = silo_data[silo_data[!, "year"] .== x,:]

            X = Matrix(silo_subset[:, covariates])
            Y = silo_subset.coll
            beta_hat = X\Y 
            Y_hat = X * beta_hat
            residuals = Y - Y_hat


            push!(data, [silo_name, string(treatment_time), string(x), string(mean(silo_data[silo_data[!, "year"] .== x, "coll"])), string(mean(residuals)), string(covariates)])
        end
    end

    # Define the file path
    filename = "trends_data_$(silo_name).csv"
    filepath = abspath(filename)

    # Save as trends_data_$(silo_name).csv
    open(filepath, "w") do file
        for row in data
            println(file, join(row, ";"))
        end
    end

    # Print the location where the .csv was saved to
    println("trends_data_$(silo_name).csv saved to")
    formatted_path = replace(filepath, "\\" => "\\\\")
    println(formatted_path)

    return filepath
end

function run_stage_two(filepath_to_empty_diff_df, silo_name, silo_data; renaming_dictionary = false)
    
    # Given a filepath to the empty_diff_df, the name of the local silo, and 
    # a dataframe of the local silo data, runs all the necessary Stage 2 functions

    empty_diff_df = read_csv_data(filepath_to_empty_diff_df)

    fill_diff_df(silo_name, empty_diff_df, silo_data, renaming_dictionary = renaming_dictionary)

    covariates = empty_diff_df[!, "covariates"][1]

    treated_or_not = empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "treat"][1]
    if treated_or_not == 1
        if "common_treatment_time" in DataFrames.names(empty_diff_df)
            treatment_time = empty_diff_df.common_treatment_time[1]
        else 
            treatment_time = empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "gvar"][1]
        end 
    elseif treated_or_not == 0
        treatment_time = "control"
    end 
    
    create_trends_df(silo_data, silo_name; covariates = covariates, treatment_time = treatment_time)

end 

### Stage 3 Functions ###
function combine_diff_data(dir_path; save_csv = false)
    
    # Collects all the filled_diff_df_... csv files
    files = readdir(dir_path)
    matching_files = filter(file -> startswith(file, "filled_diff_df_") && endswith(file, ".csv"), files)

    # Uses the read_csv_data function to read in the csv's and appends them all together
    data = read_csv_data("$dir_path\\$(matching_files[1])")
    for i in 2:length(matching_files)
        data = vcat(data, read_csv_data("$dir_path\\$(matching_files[i])"))
    end 

    # Save as csv if save_csv == true
    if save_csv == true
        save_df_as_csv("combined_diff_data.csv", data)        
    end

    # Returns df
    return data
end

function combine_trends_data(dir_path; save_csv = false)
    
    # Collects all the trends_data_... csv files
    files = readdir(dir_path)
    matching_files = filter(file -> startswith(file, "trends_data_") && endswith(file, ".csv"), files)

    # Uses the read_csv_data function to read in the csv's and appends them all together
    data = read_csv_data("$dir_path\\$(matching_files[1])")
    for i in 2:length(matching_files)
        data = vcat(data, read_csv_data("$dir_path\\$(matching_files[i])"))
    end 

    # Save as csv if save_csv == true
    if save_csv == true
        save_df_as_csv("combined_trends_data.csv", data)        
    end

    # Returns df
    return data

end 

function create_ATT_by_gt_by_silo(combined_diff_data; save_csv = false)
    
    # Creating some preallocation vectors
    treated_silo = []
    control_silo = []
    diff_times = []
    gt_vec = []
    ATT = []
    SE_ATT = []
    ATT_covariates = []
    SE_ATT_covariates = []
    covariates = []
    covariate = combined_diff_data[!, "covariates"][1]

    # For every treated silo
    for treated_group in unique(combined_diff_data[combined_diff_data[!,"treat"] .== 1,"silo_name"])

        # For every relevant time combination for said treated silo
        for time in unique(combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_group,"diff_times"])

            # For every control group that has a corresponding time combination
            for control_group in combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 0), "silo_name"]
            
                # Compute the ATT(g,t) for that treated-control combination
                control_diff = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 0) .&& (combined_diff_data[!, "silo_name"] .== control_group), "diff_estimate"]
                treated_diff = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "silo_name"] .== treated_group), "diff_estimate"]
                attgt = treated_diff - control_diff
                
                se_control_diff = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 0) .&& (combined_diff_data[!, "silo_name"] .== control_group), "diff_var"]
                se_treated_diff = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "silo_name"] .== treated_group), "diff_var"]
                
                # Compute the ATT(g,t) for that treated-control combination (with covariates)
                if covariate[1] != "none"
                    control_diff_with_covariates = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 0) .&& (combined_diff_data[!, "silo_name"] .== control_group), "diff_estimate_covariates"]
                    treated_diff_with_covariates = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "silo_name"] .== treated_group), "diff_estimate_covariates"]
                    attgt_covariates = treated_diff_with_covariates - control_diff_with_covariates
                    push!(ATT_covariates, attgt_covariates[1])

                    se_control_diff_with_covariates = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 0) .&& (combined_diff_data[!, "silo_name"] .== control_group), "diff_var_covariates"]
                    se_treated_diff_with_covariates = combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "silo_name"] .== treated_group), "diff_var_covariates"]
                    push!(SE_ATT_covariates, sqrt((se_control_diff_with_covariates + se_treated_diff_with_covariates)[1]))
                else
                    push!(ATT_covariates, "n/a")
                    push!(SE_ATT_covariates, "n/a")
                end 
                
                # Push all results to the preallocation vectors
                push!(treated_silo, treated_group)
                push!(control_silo, control_group)
                push!(diff_times, time)
                push!(gt_vec, combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "silo_name"] .== treated_group), "(g,t)"][1])
                push!(ATT, attgt[1])
                push!(SE_ATT,  sqrt((se_control_diff + se_treated_diff)[1]))              
                push!(covariates, covariate)

            end 

        end 

    end 

    # Turn the filled preallocation vectors into a dataframe
    ATT_by_gt_by_silo_df = ATT_by_gt_by_silo_df = DataFrame("treated_silo" => treated_silo, "control_silo" => control_silo,
    "diff_times" => diff_times, "ATT" => ATT, "SE_ATT" => SE_ATT, "ATT_covariates" => ATT_covariates, "SE_ATT_covariates" => SE_ATT_covariates,
    "(g,t)" => gt_vec, "covariates" => covariates)

    # If save_as_csv == true then save as a csv as well
    if save_csv == true
        save_df_as_csv("ATT_by_gt_by_silo.csv", ATT_by_gt_by_silo_df)        
    end


    return ATT_by_gt_by_silo_df

end

function compute_ATT_staggered(ATT_by_gt_by_silo)

    # Preallocation vectors
    ATTs_bysilo = []
    ATT_bysilo_covariates = []
    treated_group = []
    avg_ATT_across_silos = Vector{Any}(fill(nothing, length(unique(ATT_by_gt_by_silo[!,"treated_silo"]))))
    avg_ATT_across_silos_covariates = Vector{Any}(fill(nothing, length(unique(ATT_by_gt_by_silo[!,"treated_silo"]))))
    covariates = Vector{Any}(fill(nothing, length(unique(ATT_by_gt_by_silo[!,"treated_silo"]))))

    # Compute ATT by silo
    for group in unique(ATT_by_gt_by_silo[!,"treated_silo"])
        
        # Compute ATT by silo without covariates
        mean_att = mean(ATT_by_gt_by_silo[ATT_by_gt_by_silo[!, "treated_silo"] .== group, "ATT"])
        push!(ATTs_bysilo, mean_att)
        push!(treated_group, group)

        # Compute ATT by silo with covariates (if covariates exist)
        if ATT_by_gt_by_silo[!, "covariates"][1][1] != "none"
            mean_att_covariates = mean(ATT_by_gt_by_silo[ATT_by_gt_by_silo[!, "treated_silo"] .== group, "ATT_covariates"])
            push!(ATT_bysilo_covariates, mean_att_covariates)
        elseif ATT_by_gt_by_silo[!, "covariates"][1][1] == "none"
            push!(ATT_bysilo_covariates, "n/a")
        end 
        
    end 

    # Take average ATT across silos and display covariates
    avg_ATT_across_silos[1] = mean(ATTs_bysilo)
    if ATT_bysilo_covariates[1] != "n/a"
        avg_ATT_across_silos_covariates[1] = mean(ATT_bysilo_covariates)
        covariates[1] = ATT_by_gt_by_silo[!, "covariates"][1]
    end

    # Construct dataframe
    ATTs_bysilo_df = DataFrame(treated_silo = treated_group, ATT_bysilo = ATTs_bysilo, ATT_bysilo_covariates = ATT_bysilo_covariates,
     avg_ATT_across_silos = avg_ATT_across_silos, avg_ATT_across_silos_covariates = avg_ATT_across_silos_covariates, covariates = covariates)

     # Save as csv
     save_df_as_csv("UNDID_results.csv", ATTs_bysilo_df)

    # Return dataframe
    return ATTs_bysilo_df
end 

function compute_ATT_common(combined_diff_data)

    
    ATTs_bysilos = []
    SE_ATTs_bysilos = []    
    treated_group = []
    control_group = []
    jacknives = [] 
    jackknives_covariates = []

    # Telling the next for loop whether or not to include covariate related calculations
    if combined_diff_data[!,"covariates"][1] == ["none"]
        cov_on = false        
    else
        cov_on = true   
        ATTs_bysilos_covariates = []
        SE_ATTs_bysilos_covariates = []        
    end

    
    # This for loop calculates the ATTs by silo pairings
    if cov_on == true
        for treated_silo in combined_diff_data[combined_diff_data[!, "treat"] .== 1, "silo_name"]
            for control_silo in combined_diff_data[combined_diff_data[!, "treat"] .== 0, "silo_name"]
                push!(treated_group, treated_silo)
                push!(control_group, control_silo)
                push!(ATTs_bysilos, (combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_silo, "diff_estimate"] - combined_diff_data[combined_diff_data[!, "silo_name"] .== control_silo, "diff_estimate"])[1])            
                push!(SE_ATTs_bysilos, sqrt((combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_silo, "diff_var"] + combined_diff_data[combined_diff_data[!, "silo_name"] .== control_silo, "diff_var"])[1]))                
                push!(ATTs_bysilos_covariates, (combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_silo, "diff_estimate_covariates"] - combined_diff_data[combined_diff_data[!, "silo_name"] .== control_silo, "diff_estimate_covariates"])[1])
                push!(SE_ATTs_bysilos_covariates, sqrt((combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_silo, "diff_var_covariates"] + combined_diff_data[combined_diff_data[!, "silo_name"] .== control_silo, "diff_var_covariates"])[1]))
            end 
        end
    elseif cov_on == false
        for treated_silo in combined_diff_data[combined_diff_data[!, "treat"] .== 1, "silo_name"]
            for control_silo in combined_diff_data[combined_diff_data[!, "treat"] .== 0, "silo_name"]
                push!(treated_group, treated_silo)
                push!(control_group, control_silo)
                push!(ATTs_bysilos, (combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_silo, "diff_estimate"] - combined_diff_data[combined_diff_data[!, "silo_name"] .== control_silo, "diff_estimate"])[1])            
                push!(SE_ATTs_bysilos, sqrt((combined_diff_data[combined_diff_data[!, "silo_name"] .== treated_silo, "diff_var"] + combined_diff_data[combined_diff_data[!, "silo_name"] .== control_silo, "diff_var"])[1]))
            end 
        end
        ATTs_bysilos_covariates = Vector{Any}(fill(nothing, length(ATTs_bysilos)))
        SE_ATTs_bysilos_covariates = Vector{Any}(fill(nothing, length(ATTs_bysilos)))  
    end


    
    # This portion calculates the aggregate ATTs
    ATTs = Vector{Any}(fill(nothing, length(ATTs_bysilos)))
    ATTs_covariates = Vector{Any}(fill(nothing, length(ATTs_bysilos)))
    jackknife_SE_aggregate_ATT = Vector{Any}(fill(nothing, length(ATTs_bysilos)))
    jackknife_SE_aggregate_ATT_covariates = Vector{Any}(fill(nothing, length(ATTs_bysilos)))

    X = convert(Matrix{Float64}, hcat(fill(1, nrow(combined_diff_data)), combined_diff_data.treat))
    Y = convert(Vector{Float64}, combined_diff_data.diff_estimate)
    beta_hat = X\Y
    ATTs[1] = beta_hat[2]

    cov_check = combined_diff_data.covariates[1]
    if cov_check != ["none"]
        Y = convert(Vector{Float64}, combined_diff_data.diff_estimate_covariates)
        beta_hat = X\Y
        ATTs_covariates[1] =  beta_hat[2]    
    end

    # Calculate the jackknife SE(ATT) if there are >= 2 control and >=2 treated silos
    if nrow(combined_diff_data[combined_diff_data[!, "treat"] .== 0,:]) >= 2 && nrow(combined_diff_data[combined_diff_data[!, "treat"] .== 1,:]) >= 2
        # Compute first for the case of no covariates
        for silo in combined_diff_data[!, "silo_name"]
            X = convert(Matrix{Float64}, hcat(fill(1, nrow(combined_diff_data[combined_diff_data[!, "silo_name"] .!== silo,:])), combined_diff_data[combined_diff_data[!, "silo_name"] .!== silo, "treat"]))
            Y = convert(Vector{Float64}, combined_diff_data[combined_diff_data[!, "silo_name"] .!== silo, "diff_estimate"])
            push!(jacknives, (X\Y)[2])
        end
        jackknife_SE_aggregate_ATT[1] = 0
        for i in 1:length(jacknives)
            jackknife_SE_aggregate_ATT[1] += (ATTs[1] - jacknives[i])^2
        end
        jackknife_SE_aggregate_ATT[1] = sqrt(jackknife_SE_aggregate_ATT[1])
        
        # Compute as well for the case of covariates
        if cov_check != ["none"]
            for silo in combined_diff_data[!, "silo_name"]
                X = convert(Matrix{Float64}, hcat(fill(1, nrow(combined_diff_data[combined_diff_data[!, "silo_name"] .!== silo,:])), combined_diff_data[combined_diff_data[!, "silo_name"] .!== silo, "treat"]))
                Y = convert(Vector{Float64}, combined_diff_data[combined_diff_data[!, "silo_name"] .!== silo, "diff_estimate_covariates"])
                push!(jackknives_covariates, (X\Y)[2])
            end
            jackknife_SE_aggregate_ATT_covariates[1] = 0
            for i in 1:length(jackknives_covariates)
                jackknife_SE_aggregate_ATT_covariates[1] += (ATTs[1] - jackknives_covariates[i])^2
            end
            jackknife_SE_aggregate_ATT_covariates[1] = sqrt(jackknife_SE_aggregate_ATT_covariates[1])
        end 
    end 


    # Turn all those vectors into a Dataframe
    undid_results = DataFrame(treated_group = treated_group, control_group = control_group, ATT_bysilo = ATTs_bysilos, SE_ATT_bysilo = SE_ATTs_bysilos,
    ATT_bysilo_covariates = ATTs_bysilos_covariates, SE_ATT_bysilos_covariates = SE_ATTs_bysilos_covariates, aggregate_ATT = ATTs,
    jackknife_SE_aggregate_ATT = jackknife_SE_aggregate_ATT, aggregate_ATT_covariates = ATTs_covariates, jackknife_SE_aggregate_ATT_covariates = jackknife_SE_aggregate_ATT_covariates)

    # Save as csv
    save_df_as_csv("UNDID_results.csv", undid_results)


    # Return the dataframe
    return undid_results    
end

function run_stage_three(dir_path; save_all_csvs = false)

    combined_diff_data = combine_diff_data(dir_path, save_csv = save_all_csvs)

    combined_trends_data = combine_trends_data(dir_path, save_csv = save_all_csvs)

    if "common_treatment_time" in DataFrames.names(combined_diff_data)
        results = compute_ATT_common(combined_diff_data)
    else
        ATT_by_gt_by_silo = create_ATT_by_gt_by_silo(combined_diff_data, save_csv = save_all_csvs)
        results = compute_ATT_staggered(ATT_by_gt_by_silo)
    end 

    return results

end

### Misc Functions ###
function compute_covariance_matrix(X, sigma_sq; diff_times = false, covariates = false)
    X = convert(Matrix{Float64}, X)
    cov_beta_hat = zeros(size(X, 2), size(X, 2))
    try
        cov_beta_hat = inv(X' * X) * X' * sigma_sq * X * inv(X' * X)        
    catch ex 
        if isa(ex, SingularException)
            det_Gram = det(X' * X)
            if diff_times !== false
                if covariates != false
                    println("Warning!! Gram matrix (X' * X) for diff_times $diff_times and covariates $covariates is singular (det = $det_Gram), using pseudoinverse instead.")
                else
                    println("Warning!! Gram matrix (X' * X) for diff_times $diff_times and no covariates is singular (det = $det_Gram), using pseudoinverse instead.")
                end
            else 
                println("Warning!! Gram matrix (X' * X) is singular (det = $det_Gram), using pseudoinverse instead.")
            end 
            cov_beta_hat = pinv(X' * X) * X' * sigma_sq * X * pinv(X' * X)
        else
            println("Unexpected error occurred:", ex)
            rethrow(ex)
        end
    end

    return cov_beta_hat

end 

function save_df_as_csv(filename, df)

    filepath = abspath(filename)

    # Save as .csv and print filepath
    open(filepath, "w") do file
        println(file, join(DataFrames.names(df), ";"))
        for row in eachrow(df)
            println(file, join(row, ";"))
        end
    end

    println("$filename.csv saved to")
    formatted_path = replace(filepath, "\\" => "\\\\")
    println(formatted_path)
end 


end 
