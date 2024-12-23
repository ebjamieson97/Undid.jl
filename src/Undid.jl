module Undid

# All of these are part of base Julia with the exception of DataFrames
using Statistics
using LinearAlgebra
using DataFrames # Required for julia.ado (i.e. the Stata-Julia interface)
using DelimitedFiles
using Dates
using Random

export create_init_csv, create_diff_df, run_stage_one, # Stage 1 Functions
fill_diff_df, create_trends_df,  run_stage_two, # Stage 2 Functions
combine_diff_data, combine_trends_data, calculate_agg_att_df, run_stage_three, # Stage 3 Functions
read_csv_data, compute_covariance_matrix, save_as_csv, parse_string_to_date, parse_date_to_string, parse_freq, match_dates # Misc Functions

### Stage 1 Functions ###
function create_init_csv(names=[], start_times=[], end_times=[], treatment_times=[]; covariates = false, printmessage::Bool = false)
    
    # Check input vectors are the same length
    if length(names) == length(start_times) == length(end_times) == length(treatment_times)
        # Do nothing
    else 
        error("Input vectors for names, start_times, end_times and treatment_times must all be the same length.") 
    end 

    # Convert each element of the vectors to a string
    names, start_times, end_times, treatment_times = string.(names), string.(start_times), string.(end_times), string.(treatment_times)

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
    # This saves and prints out the file location of the created csv
    filepath = save_as_csv("init.csv", data, "array_of_arrays", printmessage)   
    return filepath
end 

function create_diff_df(csv::AbstractString; covariates = false, date_format = false, freq = false, freq_multiplier = false, weights = "standard")

    # This is for flexible date handling for the parse_string_to_date function
    possible_formats_UNDID = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm", "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy", "mm-dd-yyyy", "mmddyyyy",
    "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy", "ddmonyyyy", "yyyym00"]
    month_map_UNDID = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04", "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08", "sep" => "09", "oct" => "10",
    "nov" => "11", "dec" => "12")

    # Ensure that date_format is manually selected
    if !(date_format in possible_formats_UNDID)
        error("Please ensure the start_times, end_times and treatment_times are formatted identically and set the argument date_format equal to one of the following: [\"yyyy/mm/dd\", \"yyyy-mm-dd\", \"yyyymmdd\", \"yyyy/dd/mm\", \"yyyy-dd-mm\", \"yyyyddmm\", \"dd/mm/yyyy\", \"dd-mm-yyyy\", \"ddmmyyyy\", \"mm/dd/yyyy\", \"mm-dd-yyyy\", \"mmddyyyy\",
        \"mm/yyyy\", \"mm-yyyy\", \"mmyyyy\", \"yyyy\", \"ddmonyyyy\", \"yyyym00\"]")
    end 

    # Define the frequency of data
    if freq == false
        error("Please enter a frequency of data: 'daily', 'weekly', 'monthly', or 'yearly'")
    elseif freq == "daily"
        freq = Day(1)
    elseif freq == "weekly"
        freq = Week(1)
    elseif freq == "monthly"
        freq = Month(1)
    elseif freq == "yearly"
        freq = Year(1)
    else
        error("Please enter a frequency of data as a string: 'daily', 'monthly', 'weekly' or 'yearly'")
    end 

    # Allow for weird frequencies like every 4 days or every 7 months etc
    if freq_multiplier == false
        # do nothing
    elseif typeof(freq_multiplier) <: Integer && freq_multiplier !== 0
        freq = freq*freq_multiplier
    else
        error("freq_multiplier must be a non-zero integer or set to false.")
    end

    # Read in the init.csv and make it a dataframe
    data = readdlm(csv, ',')
    header, rows  = data[1, :], data[2:end, :]
    df = DataFrame(rows, Symbol.(header))
    
    # Converting the dates into date type objects
    try 
        df.start_time = parse_string_to_date.(string.(df.start_time), date_format, Ref(possible_formats_UNDID), Ref(month_map_UNDID))
    catch ex
        if isa(ex, ArgumentError)
            println("Ensure date_format (e.g. 'yyyy/mm/dd') reflects the actual format the dates are stored as in the init.csv.")
        else 
            throw(ex)
        end 
    end 
    df.end_time = parse_string_to_date.(string.(df.end_time), date_format, Ref(possible_formats_UNDID), Ref(month_map_UNDID))
    df.treatment_time = parse_string_to_date.(string.(df.treatment_time), date_format, Ref(possible_formats_UNDID), Ref(month_map_UNDID))

    # Setup for staggered adoption
    if length(unique(df.treatment_time)) > 2
        
        # Set up the skeleton of the full empty_diff_df
        header = ["silo_name", "gvar", "treat", "diff_times", "gt"]
        column_types = [Any, String, Int, Tuple, Tuple]
        diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])
        
        # Produces the appropriate year combinations that need to be calculated based on the init.csv
        for group in df[:, "silo_name"]
            
            times = df[df[!, "silo_name"] .== group, "start_time"][1]:freq:df[df[!, "silo_name"] .== group, "end_time"][1]       
            if df[df[!, "silo_name"] .== group, "treatment_time"][1] !== "control"
                # These are the years that matter for treated groups                
                combinations = [(t1, t2) for t1 in times, t2 in times if t1 > t2 && Dates.value(t1 - t2) >= Dates.value(freq) && Dates.value(df[df[!, "silo_name"] .== group, "treatment_time"][1] - freq) == Dates.value(t2)]
                treat = 1                
            elseif df[df[!, "silo_name"] .== group, "treatment_time"][1] == "control"
                # These are the years that matter for untreated groups
                combinations = [(t1, t2) for t1 in times, t2 in times if t1 > t2 && Dates.value(t1 - t2) >= Dates.value(freq)]
                treat = 0
            end 

            for i in 1:length(combinations)
                push!(diff_df, [group, parse_date_to_string(combinations[i][2]+freq, date_format), treat, parse_date_to_string.(combinations[i], date_format), parse_date_to_string.((combinations[i][2]+freq, combinations[i][1]), date_format)])
            end            
        end
      
        # Remove any rows for which the difference is calcuable for a treatment group but does exist for a control group
        # and vice versa
        diff_set_1 = setdiff(diff_df[diff_df[!, "treat"] .== 1, "diff_times"], diff_df[diff_df[!, "treat"] .== 0, "diff_times"])
        diff_set_2 = setdiff(diff_df[diff_df[!, "treat"] .== 0, "diff_times"], diff_df[diff_df[!, "treat"] .== 1, "diff_times"])
        unpaired_times = union(diff_set_1, diff_set_2)        
        filter!(row -> !(row.diff_times in unpaired_times), diff_df)

        # Change the , in the date tuples to a ;
        diff_df.diff_times = map(x -> join(x, ";"), diff_df.diff_times)
        diff_df[!,"gt"] = map(x -> join(x, ";"), diff_df[!,"gt"])

    # For common treatment time
    elseif length(unique(df.treatment_time)) <= 2

        # Set up the skeleton of the full empty_diff_df
        header = ["silo_name", "treat", "common_treatment_time", "start_time", "end_time"]
        column_types = [Any, Int, String, Date, Date]
        diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])
        common_treatment_time = parse_date_to_string(df[df[!, "treatment_time"] .!== "control", "treatment_time"][1], date_format)

        for silo_name in unique(df.silo_name)
            if df[df[!, "silo_name"] .== silo_name, "treatment_time"][1] !== "control"
                treat = 1
                start_time = df[df[!, "silo_name"] .== silo_name, "start_time"][1]
                end_time = df[df[!, "silo_name"] .== silo_name, "end_time"][1]
            elseif df[df[!, "silo_name"] .== silo_name, "treatment_time"][1] == "control"
                treat = 0
                start_time = df[df[!, "silo_name"] .== silo_name, "start_time"][1]
                end_time = df[df[!, "silo_name"] .== silo_name, "end_time"][1]
            end
            push!(diff_df, [silo_name, treat, common_treatment_time, start_time, end_time])
        end
        if weights == "standard"
            diff_df.weights = fill("standard", nrow(diff_df))
        else
            error("Please select a valid weighting method. Options include: \"standard\"")
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
            covariates = df[!,"covariates"][1]
            diff_df.covariates = fill(covariates, nrow(diff_df))
        else
            covariates = "none"
            diff_df.covariates = fill(covariates, nrow(diff_df))
        end
    else
        if typeof(covariates) == Vector{String}
        diff_df.covariates = fill(join(covariates, ";"), nrow(diff_df))
        else
            error("Covariates must be entered as a vector of strings or set to covariates = false to get covariate information from the init.csv.")
        end
    end 

    # Note the formatting of the date related information
    diff_df.date_format = Vector{Any}(fill(date_format, nrow(diff_df)))    
    diff_df.freq = Vector{Any}(fill(freq, nrow(diff_df))) 

    # Finally, this if loop adds the necessary rows needed to carry out the RI inference later for staggered adoption
    if "gvar" in DataFrames.names(diff_df)
        diff_df.RI = fill(0, nrow(diff_df))
        diff_df.treat = convert(Vector{Int}, diff_df.treat)
        for silo in unique(diff_df[diff_df.treat .== 1, "silo_name"])
            diff_times_to_add = setdiff(unique(diff_df.diff_times), diff_df[diff_df.silo_name .== silo ,:].diff_times)
            for time in diff_times_to_add
                filtered_df = unique(filter(row -> row.diff_times == time, diff_df)[:, Not([:silo_name, :RI, :treat])])
                gvar = filtered_df.gvar[1]
                filtered_df = filtered_df[:, Not(:gvar)]
                insertcols!(filtered_df, 1, :silo_name => [silo])
                insertcols!(filtered_df, 2, :gvar => [gvar])
                insertcols!(filtered_df, 3, :treat => [-1])
                filtered_df[!, :RI] = [1]
                append!(diff_df, filtered_df)
            end
        diff_df.start_time = Vector{Date}(undef, nrow(diff_df)) 
        diff_df.end_time = Vector{Date}(undef, nrow(diff_df)) 
        for silo in unique(diff_df.silo_name)
            diff_df[diff_df.silo_name .== silo, "start_time"] .= df[df.silo_name .== silo, "start_time"][1]
            diff_df[diff_df.silo_name .== silo, "end_time"] .= df[df.silo_name .== silo, "end_time"][1]
        end 

        end 
    end

    filepath = save_as_csv("empty_diff_df.csv", diff_df, "df", false)
    return filepath, diff_df 
    
end 

function run_stage_one(names, start_times, end_times, treatment_times; covariates = false, date_format = false, freq = false, freq_multiplier = false, weights = "standard")
    csv = create_init_csv(names, start_times, end_times, treatment_times, printmessage = true)
    return create_diff_df(csv, covariates = covariates, date_format = date_format, freq = freq, freq_multiplier = freq_multiplier, weights = weights)
end

### Stage 2 Functions ### 
function fill_diff_df(silo_name::AbstractString, empty_diff_df::DataFrame, silo_data::DataFrame; treatment_time = false, consider_covariates::Bool = true)
       
    empty_diff_df_date_format = empty_diff_df.date_format[1]

    # Allows for the option of ignoring covariates, even if specified initially in stage one.
    if consider_covariates == true
        covariates = empty_diff_df.covariates[1]
    elseif consider_covariates == false
        covariates = ["none"]
    end     

    # Compute diff (gamma) estimates and standard errors in the case of common treatment times across silos
    if "common_treatment_time" in DataFrames.names(empty_diff_df)
        if treatment_time == false
            error("Please specify a common treatment time")
        end

        start_date = Date(empty_diff_df[empty_diff_df.silo_name .== silo_name, "start_time"][1], "yyyy-mm-dd")
        end_date = Date(empty_diff_df[empty_diff_df.silo_name .== silo_name, "end_time"][1], "yyyy-mm-dd")
        
        # Use treatment time to construct dummy varaible indicating obs after treatment time
        silo_data = silo_data[(silo_data.time .>= start_date .&& silo_data.time .<= end_date),:]
        silo_data_copy = copy(silo_data)
        silo_data_copy.time = ifelse.(silo_data_copy.time .>= treatment_time, 1, 0)
        
        # Calculate weight for silo and throw to empty_diff_df
        if empty_diff_df.weights[1] == "standard"
            weight = sum(silo_data_copy.time) / nrow(silo_data_copy)
            empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "weights"] .= weight
        end 

        # Construct X and Y for regression
        X = convert(Matrix{Float64},hcat(fill(1,nrow(silo_data_copy)), silo_data_copy.time))
        Y = convert(Vector{Float64},silo_data_copy.outcome)


        
        # Tell the next for loop to calculate diff_estimate with covariates if cov_on == 2
        if covariates == ["none"]
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

            if cov_on == 1 && i == 1
                try
                    X = hcat(X, Matrix(silo_data[:, covariates]))
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
            
            # Check if either diff_combo[1] or diff_combo[2] is not found in the "time" column of silo_data
            if !(diff_combo[1] in silo_data[!, "time"]) || !(diff_combo[2] in silo_data[!, "time"])
                continue  # Skip to the next diff_combo if either date is not found
            end


            # Filters the silo down to the relevant year combinations and replaces the values with 1 (post) and 0 (pre)
            silo_subset = silo_data[(silo_data[!, "time"] .== diff_combo[1] .|| silo_data[!, "time"] .== diff_combo[2]),:]
            silo_subset.time = replace(silo_subset.time, diff_combo[1] => 1, diff_combo[2] => 0)
            
            # Define X and Y matrix for regression and add covariates if applicable
            X = convert(Matrix{Float64},hcat(fill(1, length(silo_subset.time)), silo_subset.time))
            Y = convert(Vector{Float64}, silo_subset.outcome)
            
            # Perform regression and throw gamma and var(gamma) to the empty_diff_df
            beta_hat = X\Y
            resid = Y - X*beta_hat
            sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))            
            cov_beta_hat = compute_covariance_matrix(X, sigma_sq, diff_times = diff_combo)
            
            empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_times"])), "diff_estimate"] .= beta_hat[2]
            empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& (in.(diff_combo[2], empty_diff_df[!, "diff_times"])), "diff_var"] .= cov_beta_hat[2,2]


            if (covariates == ["none"]) 
                # do nothing
            else
                # Perform regression with covariates
                columns = Symbol.(covariates)
                
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

        # Return date objects to strings
        empty_diff_df.gvar = parse_date_to_string.(empty_diff_df.gvar, empty_diff_df_date_format)
        empty_diff_df[!, "gt"] = [join((parse_date_to_string(date1, empty_diff_df_date_format), parse_date_to_string(date2, empty_diff_df_date_format)), ";") for (date1, date2) in empty_diff_df[!, "gt"]]
        empty_diff_df.diff_times = [join((parse_date_to_string(date1, empty_diff_df_date_format), parse_date_to_string(date2, empty_diff_df_date_format)), ";") for (date1, date2) in empty_diff_df.diff_times]

    end  

    # Save as csv and restructure covariates so that they are delimited by ; makes it easier to read back in from csv
    silo_name = string(silo_name)
    empty_diff_df.covariates = fill(join(empty_diff_df.covariates[1], ";"), nrow(empty_diff_df))
      
    
    filepath = save_as_csv("filled_diff_df_$silo_name.csv", empty_diff_df, "df", false)
    return filepath, empty_diff_df

       
end 

function create_trends_df(empty_diff_df::DataFrame, silo_name::AbstractString, silo_data::DataFrame, freq; covariates::Vector{String} = ["none"], treatment_time = missing, date_format::AbstractString, consider_covariates::Bool = true)    

    # Define column headers
    header = ["silo_name", "treatment_time", "time", "mean_outcome", "mean_outcome_residualized", "covariates", "date_format", "freq"]
    trends_df = DataFrame(Symbol.(header) .=> [[] for column in header])
    
    # Allows for the option of ignoring covariates, even if specified initially in stage one.
    if consider_covariates == false
        covariates = ["none"]
    end

    # Grab start and end times
    start_time = Date(empty_diff_df.start_time[1])
    end_time = Date(empty_diff_df.end_time[1])

    # Push means and time to data
    if covariates == ["none"]
        for x in start_time:freq:end_time
            push!(trends_df, [silo_name, string(treatment_time), parse_date_to_string(x, date_format), string(mean(silo_data[silo_data[!, "time"] .== x, "outcome"])), "missing", ["none"], string(date_format), string(freq)])
        end
    else        
        for x in start_time:freq:end_time       
            
            silo_subset = silo_data[silo_data[!, "time"] .== x,:]
            

            X = Matrix(silo_subset[:, covariates])
            X = convert(Matrix{Float64}, X) # Ensure matrix is of type Float64 as oppose to Any. Important for wrappers
            Y = convert(Vector{Float64},silo_subset.outcome)
           
            beta_hat = X\Y 
            Y_hat = X * beta_hat
            residuals = Y - Y_hat

            push!(trends_df, [silo_name, string(treatment_time), parse_date_to_string(x, date_format), string(mean(silo_data[silo_data[!, "time"] .== x, "outcome"])), string(mean(residuals)), covariates, string(date_format), string(freq)])
        end
    end
    
    trends_df.covariates = fill(join(trends_df.covariates[1], ";"), nrow(trends_df))

    filepath = save_as_csv("trends_data_$(silo_name).csv", trends_df, "df", false)
    return filepath, trends_df

end

function run_stage_two(filepath_to_empty_diff_df::AbstractString, silo_name::AbstractString, silo_data::DataFrame, time_column::AbstractString, outcome_column::AbstractString, date_format_local::AbstractString; renaming_dictionary = false, consider_covariates::Bool = true)
    
    # Given a filepath to the empty_diff_df, the name of the local silo, and 
    # a dataframe of the local silo data, runs all the necessary Stage 2 functions
    empty_diff_df = read_csv_data(filepath_to_empty_diff_df)
    empty_diff_df.silo_name = string.(empty_diff_df.silo_name)
    covariates = empty_diff_df[!, "covariates"][1]
    
    # First doing some pre-processing and error checking

    # Check that the silo_name is in the empty_diff_df
    if !(silo_name in unique(empty_diff_df.silo_name))
        error("Please ensure that the silo_name exists in the silo_name column of the empty_diff_df.")
    end 

    # Check that some of the arguments are strings
    if isa(time_column, String) && isa(outcome_column, String) && isa(date_format_local, String)
        # do nothing
    else 
        error("Please ensure that time_column, outcome_column, and date_format are entered as strings.")
    end

    # Rename columns if necessary
    if renaming_dictionary == false
        # do nothing
    elseif typeof(renaming_dictionary) == Dict{String, String}
        # Check that all of the dictionary values (new column names) are in the empty_diff_df and all dictionary keys are in the local silo data column names and then rename
        if all(value -> value in covariates, values(renaming_dictionary)) && all(key -> key in DataFrames.names(silo_data), keys(renaming_dictionary))
            rename!(silo_data, renaming_dictionary)
        elseif all(value -> value in covariates, values(renaming_dictionary))
            error("Ensure that the keys in the renaming_dictionary are the covariate names in the local silo data.")
        elseif all(key -> key in DataFrames.names(silo_data), keys(renaming_dictionary))
            error("Ensure that the values in the renaming_dictionary are the covariate names in the empty_diff_df.")
        else
            error("Ensure that the keys in the renaming_dictionary are the covariate names in the local silo data and that the values of the renaming_dictionary are the corresponding covariate names in the empty_diff_df.")
        end  
    else
        error("If renaming columns ensure that the typeof(renaming_dictionary) == Dict{String, String}")
    end 

    # Select only the relevant silo from the empty_diff_df also save the date format
    empty_diff_df = empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name,:]
    empty_diff_df_date_format = empty_diff_df.date_format[1]

    # Rename the time and treatment columns appropriately
    rename!(silo_data, time_column => "time")
    rename!(silo_data, outcome_column => "outcome")
    

    # This is for flexible date handling for the parse_string_to_date function
    possible_formats_UNDID = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm", "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy", "mm-dd-yyyy", "mmddyyyy",
    "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy", "ddmonyyyy", "yyyym00"]
    month_map_UNDID = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04", "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08", "sep" => "09", "oct" => "10",
    "nov" => "11", "dec" => "12")

    # Convert some string date information into Date type objects
    if "diff_times" in DataFrames.names(empty_diff_df)      
        treatment_time = false
    end     
    if "common_treatment_time" in DataFrames.names(empty_diff_df)
        treatment_time = parse_string_to_date(string(empty_diff_df.common_treatment_time[1]), empty_diff_df.date_format[1], possible_formats_UNDID, month_map_UNDID)
    end
    
    # Convert the time column to a date object
    if date_format_local == "yyyy" && typeof(silo_data.time[1]) <: AbstractFloat
        silo_data.time = Int.(silo_data.time)
    elseif date_format_local == "yyyy" && typeof(silo_data.time[1]) <: Integer
        # do nothing
    elseif !(typeof(silo_data.time[1]) <: String)
        error("Ensure your time column contains string values.")
    end
    silo_data.time = parse_string_to_date.(lowercase.(string.(silo_data.time)), date_format_local, Ref(possible_formats_UNDID), Ref(month_map_UNDID))
    
    # Grab frequency of dates
    freq = parse_freq(string(empty_diff_df.freq[1]))

    # Go through dates matching procedure if necessary
    if "diff_times" in DataFrames.names(empty_diff_df)
        empty_diff_dates = unique([x for x in empty_diff_df.diff_times for x in x])
        silo_diff_times = unique(silo_data.time)
        silo_unmatched_dates = [d for d in silo_diff_times if d ∉ empty_diff_dates]
        if length(silo_unmatched_dates) > 0 
            date_map = match_dates(empty_diff_dates, silo_unmatched_dates, freq)
            silo_data.time = get.(Ref(date_map), silo_data.time, silo_data.time)
        end 
    end

    filepath_and_diff_df = fill_diff_df(silo_name, empty_diff_df, silo_data, treatment_time = treatment_time, consider_covariates = consider_covariates)

    treated_or_not = empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (empty_diff_df[!, "treat"] .!= -1), "treat"][1]
    if treated_or_not == 1
        if "common_treatment_time" in DataFrames.names(empty_diff_df)
            treatment_time = empty_diff_df.common_treatment_time[1]
        else 
            treatment_time = empty_diff_df[(empty_diff_df[!, "silo_name"] .== silo_name) .&& (empty_diff_df[!, "treat"] .!= -1), "gvar"][1]
        end 
    elseif treated_or_not == 0
        treatment_time = "control"
    end 
    
    filepath_and_trends_df = create_trends_df(empty_diff_df, silo_name, silo_data, freq, covariates = covariates, treatment_time = treatment_time, date_format = empty_diff_df_date_format, consider_covariates = consider_covariates)

    return filepath_and_diff_df, filepath_and_trends_df

end 

### Stage 3 Functions ###
function combine_diff_data(dir_path::AbstractString; save_csv::Bool = false, interpolation = false)
    
    # Collects all the filled_diff_df_... csv files
    files = readdir(dir_path)
    matching_files = filter(file -> startswith(file, "filled_diff_df_") && endswith(file, ".csv"), files)

    # Uses the read_csv_data function to read in the csv's and appends them all together
    data = read_csv_data("$dir_path\\$(matching_files[1])")
    for i in 2:length(matching_files)
        data = vcat(data, read_csv_data("$dir_path\\$(matching_files[i])"))
    end 

    # This block performs linear interpolation/extrapolation for diff_estimate if interpolation is set to true
    if any(==("missing"), data.diff_estimate)
        indices = findall(==( "missing"), data.diff_estimate)
        println("Missing diff_estimate for rows: $indices")
        if interpolation == false
            println("Consider setting interpolation = \"linear_function\".")
        elseif interpolation == "linear_function"
            data.local_period = Vector{Any}(fill("missing", nrow(data)))
            for index in indices
                silo = data[index,:].silo_name
                gvar = data[index,:].gvar
                gt = data[index,"gt"]
                println("Using a linear function to fill missing values of diff_estimate for silo: $silo (g,t): $gt")
                periods = sort(unique([x for x in data[(data.silo_name .== silo) .&& (data.gvar .== gvar),"gt"] for x in x]))
                for i in 1:length(periods)
                    data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& ((getindex.(data."gt", 2)) .== periods[i]), "local_period"] .= i
                end 
                temp_df = data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate .!= "missing"),:]
                temp_df_missing = data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate .== "missing"),:]
                if nrow(temp_df) == 0
                    error("Silo: $silo for gvar: $gvar has no non-missing entries, unable to perform linear interpolation or extrapolation.")
                elseif nrow(temp_df) == 1
                    value = temp_df.diff_estimate[1]
                    indices = findall((data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate .== "missing"))
                    data.diff_estimate[indices] .= value
                else                
                    Y = convert(Vector{Float64}, temp_df.diff_estimate)
                    X = convert(Matrix{Float64}, hcat(fill(1, nrow(temp_df)), temp_df.local_period))
                    beta_hat = X \ Y
                    interpolations = convert(Matrix{Float64}, hcat(fill(1, length(temp_df_missing.local_period)), temp_df_missing.local_period)) * beta_hat
                    indices = findall((data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate .== "missing"))
                    data.diff_estimate[indices] .= interpolations
                end
            end 
        else
            error("Set interpolation to false to \"linear_function\".")
        end 
    end

    # This block performs linear interpolation/extrapolation for diff_estimate_covariates if interpolation is on
    if nrow(data[data.diff_estimate_covariates .!= "missing",:]) > 0
        indices = findall(==( "missing"), data.diff_estimate_covariates)
        if length(indices) > 0 
            println("Missing diff_estimate_covariates for rows: $indices")
            if interpolation == false
                println("Consider setting interpolation = \"linear_function\".")
            elseif interpolation == "linear_function"
                data.local_period = Vector{Any}(fill("missing", nrow(data)))
                for index in indices
                    silo = data[index,:].silo_name
                    gvar = data[index,:].gvar
                    gt = data[index,"gt"]
                    println("Using a linear function to fill missing values of diff_estimate_covariates for silo: $silo (g,t): $gt")
                    periods = sort(unique([x for x in data[(data.silo_name .== silo) .&& (data.gvar .== gvar),"gt"] for x in x]))
                    for i in 1:length(periods)
                        data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& ((getindex.(data."gt", 2)) .== periods[i]), "local_period"] .= i
                    end 
                    temp_df = data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate_covariates .!= "missing"),:]
                    temp_df_missing = data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate_covariates .== "missing"),:]

                    if nrow(temp_df) == 0
                        error("Silo: $silo for gvar: $gvar has no non-missing entries, unable to perform linear interpolation or extrapolation.")
                    elseif nrow(temp_df) == 1
                        value = temp_df.diff_estimate_covariates[1]
                        indices = findall((data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate_covariates .== "missing"))
                        data.diff_estimate_covariates[indices] .= value
                    else                
                        Y = convert(Vector{Float64}, temp_df.diff_estimate_covariates)
                        X = convert(Matrix{Float64}, hcat(fill(1, nrow(temp_df)), temp_df.local_period))
                        beta_hat = X \ Y
                        interpolations = convert(Matrix{Float64}, hcat(fill(1, length(temp_df_missing.local_period)), temp_df_missing.local_period)) * beta_hat
                        indices = findall((data.silo_name .== silo) .&& (data.gvar .== gvar) .&& (data.diff_estimate_covariates .== "missing"))
                        data.diff_estimate_covariates[indices] .= interpolations
                    end
                end 
             
            else
                error("Set interpolation to false or to \"linear_function\".")
            end 
        end 
    
    end

    # Save as csv if save_csv == true
    if save_csv == true
        copy_df = copy(data)
        copy_df.covariates = fill(join(copy_df.covariates[1], ";"), nrow(copy_df))           
        # Return date objects to strings
        copy_df_date_format = data.date_format[1]

        if "diff_times" in DataFrames.names(copy_df)
            copy_df[!, "gt"] = [join((parse_date_to_string(date1, copy_df_date_format), parse_date_to_string(date2, copy_df_date_format)), ";") for (date1, date2) in copy_df[!, "gt"]]
            copy_df.diff_times = [join((parse_date_to_string(date1, copy_df_date_format), parse_date_to_string(date2, copy_df_date_format)), ";") for (date1, date2) in copy_df.diff_times]
            copy_df.gvar = parse_date_to_string.(copy_df.gvar, copy_df_date_format)
        end 


        save_as_csv("combined_diff_data.csv", copy_df, "df")        
    end

    # Returns df
    return data
end

function combine_trends_data(dir_path::AbstractString; save_csv::Bool = false)
    
    # Collects all the trends_data_... csv files
    files = readdir(dir_path)
    matching_files = filter(file -> startswith(file, "trends_data_") && endswith(file, ".csv"), files)

    # Uses the read_csv_data function to read in the csv's and appends them all together
    data = read_csv_data("$dir_path\\$(matching_files[1])")
    for i in 2:length(matching_files)
        data = vcat(data, read_csv_data("$dir_path\\$(matching_files[i])"))
    end 

    if any(isnan, data.mean_outcome)
        data= filter(row -> !isnan(row.mean_outcome), data)
    end

    # Save as csv if save_csv == true
    if save_csv == true
        copy_df = copy(data)
        copy_df.covariates = fill(join(copy_df.covariates[1], ";"), nrow(copy_df))           
        # Return date objects to strings
        copy_df_date_format = data.date_format[1]
        copy_df.time = parse_date_to_string.(copy_df.time, copy_df_date_format)        
        copy_df[copy_df.treatment_time .!= "control", "treatment_time"] = parse_date_to_string.(copy_df[copy_df.treatment_time .!= "control",:].treatment_time, copy_df_date_format)
        save_as_csv("combined_trends_data.csv", copy_df, "df")        
    end

    # Returns df
    return data

end 

function run_stage_three(dir_path::AbstractString; agg::AbstractString = "silo", covariates::Bool = false, save_all_csvs::Bool = false, interpolation = false, weights::Bool = true, nperm = 1000)

    # Set nperm to be a positive integer
    nperm = Int(round(nperm))
    if nperm < 1
        error("`nperm` must be an integer > 0")
    end 
    
    combined_diff_data = combine_diff_data(dir_path, save_csv = save_all_csvs, interpolation = interpolation) 
    combined_diff_data.treat = convert(Vector{Float64}, combined_diff_data.treat)
    if covariates == false
        combined_diff_data.y = convert(Vector{Float64}, combined_diff_data.diff_estimate)
    elseif covariates == true
        try
            combined_diff_data.y = convert(Vector{Float64}, combined_diff_data.diff_estimate_covariates)
        catch ex
            if isa(ex, MethodError)
                error("No covariates. Set covariates = false.")
            else
                rethrow(ex)
            end 
        end 
    end

    # Generate all the necessary matrices to do the randomization inference
    if "RI" in DataFrames.names(combined_diff_data)
        # Filter for just the data that is reflective of the silos true treatment status 
        diff_matrix_no_RI = combined_diff_data[combined_diff_data.RI .== 0,:]
        
        # Calculate att by specified aggregation method 
        results = calculate_agg_att_df(diff_matrix_no_RI; agg = agg, covariates = covariates, save_all_csvs = save_all_csvs, printinfo = true, weights = weights, nperm = nperm)  

        # Compute p-value based on randomization inference and add to results as a column
        original_ATT = results.agg_att[1]
        gts = unique(combined_diff_data.gt)
        RI_ATT = []
        treated_silos = unique(combined_diff_data[combined_diff_data.treat .== 1, ["silo_name", "gvar", "treat"]])
        treated_silos.gvar = convert(Vector{Any}, treated_silos.gvar)
        control_silos = unique(combined_diff_data[combined_diff_data.treat .== 0, ["silo_name", "treat"]])
        control_silos.gvar = fill("not_treated", nrow(control_silos))
        init = vcat(treated_silos, control_silos)
        unique_permutations = binomial(nrow(init), nrow(treated_silos))
        if nperm > unique_permutations
            @warn "nperm was set to $nperm, but the number of unique permutations is only $unique_permutations. Overriding nperm to $unique_permutations."
            nperm = unique_permutations
        end 
        # Create hash table, log initial gvar 
        # Enforce each permutation to be unique
        key_type = typeof(hash("example"))
        seen = Dict{key_type, Bool}()
        seen[hash(init.gvar)] = true
        i = 1
        while i < nperm
            new_perm = shuffle(init.gvar)
            key = hash(new_perm)
            if !haskey(seen, key)
                init[:, string("gvar_randomized_", i)] = new_perm
                seen[key] = true
                i += 1
            end                
        end

        # Create preallocation vector to store agg_ATT from each permutation
        RI_ATT = Vector{Float64}(undef, nperm)

        # Loop through gvar permutations nperm times
        for i in 1:nperm-1
            # Selects control silos
            mask_gvar = isa.(init[:, string("gvar_randomized_", i)], String)

            # Grabs gvars in the same order as the silos they are assigned to
            new_gvars = init[.!mask_gvar, string("gvar_randomized_", i)]

            new_treated_silos = init.silo_name[.!mask_gvar]
            new_control_silos = init.silo_name[mask_gvar]

            ri_df = combined_diff_data[in.(combined_diff_data.silo_name, Ref(new_control_silos)), :]
            ri_df.treat .= 0

            # Construct the rest of this iteration's treated silos
            for j in 1:length(new_gvars)
                silo = new_treated_silos[j]
                gvar = new_gvars[j]

                mask = (combined_diff_data.silo_name .== silo) .& (combined_diff_data.gvar .== gvar)
                ri_df = vcat(ri_df, combined_diff_data[mask, :])
                ri_mask = (ri_df.silo_name .== silo) .& (ri_df.gvar .== gvar)
                ri_df.treat[ri_mask] .= 1
            end
            
            # Compute agg_ATT for agg == "silo"
            if agg == "silo"
                att_vector = Vector{Float64}(undef, length(new_treated_silos))
                for k in 1:length(new_treated_silos)
                    silo = new_treated_silos[k]
                    treated_mask = (ri_df.silo_name .== silo) .& (ri_df.treat .== 1)
                    gvar_treated = ri_df[treated_mask, "gvar"][1]
                    control_mask = (ri_df.treat .== 0) .& (ri_df.gvar .== gvar_treated)
                    subset_df = vcat(ri_df[treated_mask, :], ri_df[control_mask, :])
                    X = convert(Matrix{Float64}, hcat(ones(nrow(subset_df)), subset_df.treat))
                    Y = convert(Vector{Float64}, subset_df.y)
                    att_vector[k] = (X \ Y)[2]
                end 
                RI_ATT[i] = mean(att_vector)
            end 
            
            # Compute agg_ATT if agg == "g"
            if agg == "g"
                att_vector = Vector{Float64}(undef, length(new_gvars))
                for k in 1:length(new_gvars)
                    gvar = new_gvars[k]
                    mask = ri_df.gvar .== gvar
                    X = ri_df[mask, "treat"]
                    X = convert(Matrix{Float64}, hcat(ones(length(X)), X))
                    Y = convert(Vector{Float64}, ri_df[mask, "y"])
                    att_vector[k] = (X \ Y)[2]
                end 
                RI_ATT[i] = mean(att_vector)
            end 

            # Compute agg_ATT if agg == "gt"
            if agg == "gt"
                att_vector = Vector{Float64}(undef, length(gts))
                for k in 1:length(gts)
                    gt = gts[k]
                    mask = map(x -> all(x .== gt), ri_df.gt)
                    X = ri_df[mask, "treat"]
                    X = convert(Matrix{Float64}, hcat(ones(length(X)), X))
                    Y = convert(Vector{Float64}, ri_df[mask, "y"])
                    att_vector[k] = (X \ Y)[2]
                end 
                RI_ATT[i] = mean(att_vector)
            end 
        end 

        p_value_RI = sum(abs.(RI_ATT) .> abs(results.agg_att[1])) / length(RI_ATT)
        results.p_value_RI = vcat(p_value_RI, fill(NaN, nrow(results)-1))
        save_as_csv("UNDID_results.csv", results, "df", true)
        return results
    end 

    results = calculate_agg_att_df(combined_diff_data; agg = agg, covariates = covariates, save_all_csvs = save_all_csvs, weights = weights, nperm = nperm) 
    save_as_csv("UNDID_results.csv", results, "df", true)
    return results 
end

function calculate_agg_att_df(combined_diff_data::DataFrame; agg::AbstractString = "silo", covariates::Bool = false, save_all_csvs::Bool = false, printinfo::Bool = true, weights::Bool = true, nperm = 1000)

    if "common_treatment_time" in DataFrames.names(combined_diff_data)
        println("Calcualting aggregate ATT with covariates set to $covariates.")
        nrow_combined_diff_data = nrow(combined_diff_data)
        if weights == true && nrow_combined_diff_data > 2
            W = convert(Vector{Float64}, combined_diff_data.weights)
            sum_W = sum(W)
            W = W ./ sum_W
            W = Diagonal(W)
        elseif weights == false || nrow_combined_diff_data == 2 
            W = Diagonal(fill(1, nrow_combined_diff_data))
        end 
        X = hcat(fill(1.0, nrow_combined_diff_data), combined_diff_data.treat)
        Y = combined_diff_data.y
        beta_hat = (inv(X' * W * X) * X' * W * Y)
        resid = Y - X*beta_hat
        sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
        ATT_se = sqrt(compute_covariance_matrix(X, sigma_sq)[2,2])
        ATT = beta_hat[2]
        treatment_time = combined_diff_data.common_treatment_time[1]
        results = DataFrame(treatment_time = treatment_time, agg_att = ATT, agg_att_se = ATT_se)
        # Compute jackknife SE when there are at least 2 controls and 2 treatment silos
        if (sum(combined_diff_data.treat .== 1) >= 2 && sum(combined_diff_data.treat .== 0) >= 2)
            jackknives_common = []
            if sum(combined_diff_data.treat .== 1) == 1
                silo_names = combined_diff_data[combined_diff_data.treat .== 0, "silo_name"]
            elseif sum(combined_diff_data.treat .== 0) == 1
                silo_names = combined_diff_data[combined_diff_data.treat .== 1, "silo_name"]
            else 
                silo_names = combined_diff_data.silo_name
            end
            for silo in silo_names
                subset = combined_diff_data[combined_diff_data.silo_name .!= silo,:]
                X_sub = hcat(fill(1.0, length(subset.treat)), subset.treat)
                Y_sub = subset.y
                W_sub = convert(Vector{Float64}, subset.weights)
                if weights == true 
                    sum_W_sub = sum(W_sub)
                    W_sub = W_sub ./ sum_W_sub
                    W_sub = Diagonal(W_sub)
                elseif weights == false
                    W_sub = Diagonal(fill(1, nrow_combined_diff_data-1))
                end

                push!(jackknives_common, (inv(X_sub' * W_sub * X_sub) * X_sub' * W_sub * Y_sub)[2])
            end 
            jackknife_SE = sqrt(sum((jackknives_common .- ATT).^2) * ((length(jackknives_common) - 1)/length(jackknives_common)))
            results.jackknife_se = [jackknife_SE]
            results.dof = [length(jackknives_common)-1]
            results.dof = Float64.(results.dof)

        elseif sum(combined_diff_data.treat .== 1) == 1 && sum(combined_diff_data.treat .== 0) == 1
            if covariates == false 
                SE = sqrt(sum(convert(Vector{Float64}, combined_diff_data.diff_var)))
            elseif covariates == true 
                SE = sqrt(sum(convert(Vector{Float64}, combined_diff_data.diff_var_covariates)))
            end
            select!(results, Not(:agg_att_se))
            results.stderr = [SE]
        end


        # Do randomization inference procedure
        unique_permutations = binomial(nrow(combined_diff_data), sum(combined_diff_data.treat .== 1))
        if nperm > unique_permutations
            @warn "nperm was set to $nperm, but the number of unique permutations is only $unique_permutations. Overriding nperm to $unique_permutations."
            nperm = unique_permutations
        end 
        ri_atts = Vector{Float64}(undef, nperm)
        intercept_column = ones(nrow(combined_diff_data))
        
        # Enforce unique randomizations
        perm_matrix = Matrix{Float64}(undef, nrow(combined_diff_data), nperm)
        key_type = typeof(hash(X[:,2]))
        seen = Dict{key_type, Bool}()
        seen[hash(X[:,2])] = true
        i = 1
        while i < nperm
            new_perm = shuffle(X[:, 2])
            key = hash(new_perm)
            if !haskey(seen, key)
                perm_matrix[:, i] = new_perm
                seen[key] = true
                i += 1
            end 
        end
        for i in 1:nperm - 1
            X_ri = hcat(intercept_column, perm_matrix[:, i])
            beta_hat = (inv(X_ri' * W * X_ri) * X_ri' * W * Y)
            ri_atts[i] = beta_hat[2]
        end 

        ri_pval = sum(abs.(ri_atts) .> abs(results.agg_att[1])) / length(ri_atts)
        results.p_value_RI = [ri_pval]
                
        return results
    else
        if !(agg == "silo" || agg == "g" || agg == "gt")
            error("Please specify a weighting scheme for ATT: \"silo\", \"gt\", or \"g\".")
        end
        if printinfo == true
            println("Calcualting aggregate ATT with covariates set to $covariates and aggregation method set to $agg.")
        end 

        if agg == "silo"            
            treated_silos = unique(combined_diff_data[combined_diff_data.treat .== 1, "silo_name"])
            silos_att = []
            silos_n = []
            jack_n = []
            silos_att_se = []
            silos_att_se_jack = []
            for silo in treated_silos
                treated = combined_diff_data[combined_diff_data.silo_name .== silo,:]
                control = combined_diff_data[combined_diff_data.treat .== 0,:]            
                gvar_treated = treated.gvar[1]
                control = control[control.gvar .== gvar_treated,:]          
                local_df = vcat(control, treated)
                push!(silos_n, nrow(local_df))
                X = hcat(fill(1.0, length(local_df.treat)), local_df.treat)
                Y = local_df.y
                beta_hat = X \ Y
                resid = Y - X*beta_hat
                sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
                push!(silos_att_se, sqrt(compute_covariance_matrix(X, sigma_sq)[2,2]))
                push!(silos_att, beta_hat[2])

                jackknife_sub = []           
                for i in 1:length(Y)
                    try
                        X_sub = vcat(X[1:i-1, :], X[i+1:end, :])
                        Y_sub = vcat(Y[1:i-1], Y[i+1:end])
                        beta_hat_sub = X_sub \ Y_sub
                        push!(jackknife_sub, beta_hat_sub[2])
                    catch e
                        continue 
                    end 
                end
                push!(jack_n, length(jackknife_sub))
                push!(silos_att_se_jack ,sqrt(sum((jackknife_sub .- beta_hat[2]).^2) * ((length(jackknife_sub) - 1) / length(jackknife_sub))))
            end
            jackknives_silo = []
            for i in 1:length(silos_att)
                # Exclude the i-th element and calculate the mean of the remaining elements
                mean_excluding_i = mean([silos_att[1:i-1]..., silos_att[i+1:end]...])
                # Push the result to the jackknives vector
                push!(jackknives_silo, mean_excluding_i)
            end
            ATT_silo = mean(silos_att)
            ATT_silo_se = sqrt(var(silos_att)/length(silos_att))
            jackknife_SE = sqrt(sum((jackknives_silo .-ATT_silo).^2) * ((length(jackknives_silo) - 1)/length(jackknives_silo)))
            results = DataFrame(silos = unique(combined_diff_data[combined_diff_data.treat .== 1, "silo_name"]), silo_n = Float64.(silos_n), jack_n = Float64.(jack_n), att_s = Float64.(silos_att), att_s_se = Float64.(silos_att_se), att_s_se_jackknife = Float64.(silos_att_se_jack),
                                agg_att = vcat([ATT_silo], fill(NaN, length(silos_att) - 1)), agg_att_se = vcat([ATT_silo_se], fill(NaN, length(silos_att) - 1)), jackknife_se = vcat([jackknife_SE], fill(NaN, length(silos_att) - 1)))          
        elseif agg == "gt" || agg == "g"            
            ATT_vec = []
            ATT_n = []
            jack_n = []
            ATT_vec_se = []
            ATT_vec_se_jack = []
            gt_vec = []
            for gt in unique(combined_diff_data[!, "gt"])                
                subset = filter(row -> row[Symbol("gt")] == gt, combined_diff_data)
                push!(ATT_n, nrow(subset))
                X = hcat(fill(1.0, length(subset.treat)),subset.treat)
                Y = subset.y
                beta_hat =  X \ Y
                resid = Y - X*beta_hat
                sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
                push!(ATT_vec_se, sqrt(compute_covariance_matrix(X, sigma_sq)[2,2]))
                push!(ATT_vec, beta_hat[2])                
                push!(gt_vec, gt)
                jackknife_sub = []
                for i in 1:length(Y)
                    try
                        X_sub = vcat(X[1:i-1, :], X[i+1:end, :])
                        Y_sub = vcat(Y[1:i-1], Y[i+1:end])
                        beta_hat_sub = X_sub \ Y_sub
                        push!(jackknife_sub, beta_hat_sub[2])
                    catch e
                        continue 
                    end 
                end
                push!(jack_n, length(jackknife_sub))
                push!(ATT_vec_se_jack ,sqrt(sum((jackknife_sub .- beta_hat[2]).^2) * ((length(jackknife_sub) - 1) / length(jackknife_sub))))

            end
            if agg == "gt"
                # Compute the jackknife means
                jackknives_gt = []
                for i in 1:length(ATT_vec)
                    # Exclude the i-th element and calculate the mean of the remaining elements
                    mean_excluding_i = mean([ATT_vec[1:i-1]..., ATT_vec[i+1:end]...])
                    # Push the result to the jackknives vector
                    push!(jackknives_gt, mean_excluding_i)
                end
                ATT_gt = mean(ATT_vec)
                ATT_gt_se = sqrt(var(ATT_vec)/length(ATT_vec))
                jackknife_SE = sqrt(sum((jackknives_gt .-ATT_gt).^2) * ((length(jackknives_gt) - 1)/length(jackknives_gt)))
                results = DataFrame(gt = gt_vec, gt_n = Float64.(ATT_n), jack_n = Float64.(jack_n), att_gt = Float64.(ATT_vec), att_gt_se = Float64.(ATT_vec_se), att_gt_se_jackknife = Float64.(ATT_vec_se_jack), agg_att = vcat([ATT_gt], fill(NaN, length(ATT_vec) - 1)), agg_ATT_se = vcat([ATT_gt_se], fill(NaN, length(ATT_vec) - 1)), jackknife_se = vcat([jackknife_SE], fill(NaN, length(ATT_vec) - 1)))
                results.gt = [join((parse_date_to_string(date1, combined_diff_data.date_format[1]), parse_date_to_string(date2, combined_diff_data.date_format[1])), ";") for (date1, date2) in results[!, "gt"]]
            end 
            if agg == "g"                
                gvars = sort(unique(combined_diff_data.gvar))
                ATT_g = []
                ATT_g_n = []
                jack_n = []
                ATT_g_se = []
                ATT_g_se_jack = []
                for g in gvars
                    subset = filter(row -> row.gvar == g, combined_diff_data)
                    push!(ATT_g_n, nrow(subset))
                    X = convert(Matrix{Float64}, hcat(fill(1, nrow(subset)), subset.treat))
                    Y = convert(Vector{Float64}, subset.diff_estimate)
                    beta_hat = X \ Y 
                    push!(ATT_g, beta_hat[2])
                    resid = Y - X*beta_hat
                    sigma_sq = sum(resid.^2) / (length(resid) - length(beta_hat))
                    push!(ATT_g_se, sqrt(compute_covariance_matrix(X, sigma_sq)[2,2]))
                    jackknife_sub = []
                    for i in 1:length(Y)
                        try
                            X_sub = vcat(X[1:i-1, :], X[i+1:end, :])
                            Y_sub = vcat(Y[1:i-1], Y[i+1:end])
                            beta_hat_sub = X_sub \ Y_sub
                            push!(jackknife_sub, beta_hat_sub[2])
                        catch e
                            continue 
                        end 
                    end
                    push!(jack_n, length(jackknife_sub))
                    push!(ATT_g_se_jack ,sqrt(sum((jackknife_sub .- beta_hat[2]).^2) * ((length(jackknife_sub) - 1) / length(jackknife_sub))))
                end 
                # Compute the jackknife means
                jackknives_g = []
                for i in 1:length(ATT_g)
                    # Exclude the i-th element and calculate the mean of the remaining elements
                    mean_excluding_i = mean([ATT_g[1:i-1]..., ATT_g[i+1:end]...])
                    # Push the result to the jackknives vector
                    push!(jackknives_g, mean_excluding_i)
                end
                agg_ATT = mean(ATT_g)
                agg_ATT_se = sqrt(var(ATT_g)/length(ATT_g))
                jackknife_SE = sqrt(sum((jackknives_g .- agg_ATT).^2) * ((length(jackknives_g) - 1)/length(jackknives_g)))
                results = DataFrame(g = parse_date_to_string.(gvars, combined_diff_data.date_format[1]), g_n = Float64.(ATT_g_n), jack_n = Float64.(jack_n), att_g = Float64.(ATT_g), att_g_se = Float64.(ATT_g_se), att_g_se_jackknife = Float64.(ATT_g_se_jack),
                                    agg_att = vcat([agg_ATT], fill(NaN, length(ATT_g) - 1)), agg_ATT_se = vcat([agg_ATT_se], fill(NaN, length(ATT_g) - 1)), jackknife_se = vcat([jackknife_SE], fill(NaN, length(ATT_g) - 1)))
            end
        end 
        return results
    end   

end   

### Misc Functions ###
function read_csv_data(filepath_to_csv_data::AbstractString)

    # This is for flexible date handling for the parse_string_to_date function
    possible_formats_UNDID = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm", "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy", "mm-dd-yyyy", "mmddyyyy",
    "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy", "ddmonyyyy", "yyyym00"]
    month_map_UNDID = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04", "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08", "sep" => "09", "oct" => "10",
    "nov" => "11", "dec" => "12")

    # Reads in the empty_diff_df.csv
    data = readdlm(filepath_to_csv_data, ',')
    header = data[1, :]    
    rows = data[2:end, :]
    df_readin = DataFrame(rows, Symbol.(header))    
    if "silo_name" in DataFrames.names(df_readin)
        df_readin.silo_name = string.(df_readin.silo_name)
    end
    if "covariates" in DataFrames.names(df_readin)
        transform!(df_readin, :covariates => ByRow(s -> [String(sub) for sub in split(replace(replace(strip(s, ['[' , ']', '"']), "\"" => ""), " " => ""), ";")]) => :covariates)
    end
    if "diff_times" in DataFrames.names(df_readin)        
        df_readin_date_format = df_readin.date_format[1]

        df_readin.diff_times = split.(df_readin.diff_times, ";")        
        transform!(df_readin, :diff_times => ByRow(time -> (parse_string_to_date(time[1], df_readin_date_format, possible_formats_UNDID, month_map_UNDID), parse_string_to_date(time[2], df_readin_date_format, possible_formats_UNDID, month_map_UNDID))) => :diff_times)

        if "gvar" in DataFrames.names(df_readin)
            df_readin.gvar = parse_string_to_date.(string.(df_readin.gvar), df_readin_date_format, Ref(possible_formats_UNDID), Ref(month_map_UNDID))
        end 
        df_readin[!, "gt"] = split.(df_readin[!, "gt"], ";")
        transform!(df_readin, :"gt" => ByRow(time -> (parse_string_to_date(time[1], df_readin_date_format, possible_formats_UNDID, month_map_UNDID), parse_string_to_date(time[2], df_readin_date_format, possible_formats_UNDID, month_map_UNDID))) => :"gt")
    end 
    if "time" in DataFrames.names(df_readin)
        
        df_readin.time = parse_string_to_date.(string.(df_readin.time), df_readin.date_format[1], Ref(possible_formats_UNDID), Ref(month_map_UNDID))

        # Now based on Date object we can created a column "period" that will translate each date to an associated integer period
        dates_vector  = sort(unique(df_readin.time))
        date_dict = Dict(zip(dates_vector, 1:length(dates_vector)))
        df_readin.period = [date_dict[date] for date in df_readin.time]

        # Likewise we define the treatment_time as a period as well
        df_readin[df_readin.treatment_time .!= "control", "treatment_time"] = parse_string_to_date.(string.(df_readin[df_readin.treatment_time .!= "control",:].treatment_time), df_readin.date_format[1], Ref(possible_formats_UNDID), Ref(month_map_UNDID))
        df_readin.treatment_period = [x isa Date ? date_dict[x] : missing for x in df_readin.treatment_time]
    end 
    return df_readin
end

function compute_covariance_matrix(X, sigma_sq; diff_times = false, covariates = false)
    X = convert(Matrix{Float64}, X)
    cov_beta_hat = zeros(size(X, 2), size(X, 2))

    try      
        cov_beta_hat = inv(X' * X) * X' * sigma_sq * X * inv(X' * X)        
    catch ex 
        if (isa(ex, SingularException)) || (isa(ex, LAPACKException)) 
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

function save_as_csv(filename, data, datatype::AbstractString, printmessage::Bool = true)

    # Get the absolute filepath i.e. the save location for the csv
    filepath = abspath(filename)
    
    # Save the data as a , delimited csv (slightly different loop depending on the data structure)
    if datatype == "array_of_arrays"        
        open(filepath, "w") do file
            for row in data
                println(file, join(row, ","))
            end
        end        
    elseif datatype == "df"        
        open(filepath, "w") do file
            println(file, join(DataFrames.names(data), ","))
            for row in eachrow(data)
                println(file, join(row, ","))
            end
        end
    else 
        error("Specify if you are saving an 'array_of_arrays' or a 'df'.")
    end 

    # if printmessage == true then print filepath
    println("$filename saved to")
    if printmessage == true
        formatted_path = replace(filepath, "\\" => "\\\\")
        println(formatted_path)
    end 
    

    return filepath
end 

function parse_string_to_date(date, date_format::AbstractString, possible_formats = false, month_map = false)
    
    # Ensure that a dictionary relating month abbrevations to integers from 1 to 12.
    # Likewise ensure that an array of possible date formats is passed. 
    # This function is always called from create_diff_df, run_stage_two, or run_stage_three so this shouldn't be an issue as the 
    # dictionary and array will be specified within those functions.
    if possible_formats == false || month_map == false
        error("Please pass a dictionary such as Dict(\"jan\" => \"01\", \"feb\" => \"02\", \"mar\" => \"03\", \"apr\" => \"04\", \"may\" => \"05\", \"jun\" => \"06\", \"jul\" => \"07\", \"aug\" => \"08\", \"sep\" => \"09\", \"oct\" => \"10\", \"nov\" => \"11\", \"dec\" => \"12\") as well as the possible formats [\"yyyy/mm/dd\", \"yyyy-mm-dd\", \"yyyymmdd\", \"yyyy/dd/mm\", \"yyyy-dd-mm\", \"yyyyddmm\", \"dd/mm/yyyy\", \"dd-mm-yyyy\", \"ddmmyyyy\", \"mm/dd/yyyy\", \"mm-dd-yyyy\", \"mmddyyyy\",
        \"mm/yyyy\", \"mm-yyyy\", \"mmyyyy\", \"yyyy\", \"ddmonyyyy\", \"yyyym00\"]")
    end
    
    # First check if the date is a lack-of-a-date and if so simply return "control"
    if date in ["control", "n/a", "never", "not treated", "not_treated", "", "Control", "CONTROL", "N/A", "NEVER", "NOT TREATED",
        "not_treated", "NOT_TREATED", "missing", "nothing"]
        output = "control"
        return output
    end 

    # Convert strings to date objects depending on the specified date_format
    if date_format == "ddmonyyyy"
        # The first case is dealing with formats such as 25dec1991 
        # which is what dates sometimes look like when converting to strings in Stata
        date = lowercase(date)
        day = date[1:2]
        month_str = date[3:5]
        year = date[6:end]
        month = month_map[month_str]
        output = Date("$day/$month/$year", "dd/mm/yyyy")
    elseif date_format == "yyyym00"
        # This is for another common Stata formatting of dates when converted to strings
        date = lowercase(date)
        info = split(date, "m")        
        output = Date("$(info[2])/$(info[1])", "mm/yyyy")
    elseif date_format in possible_formats  
        # Other formats are handled natively by Dates      
        output = Date(date, date_format)
    else
        error("Please specify a date_format listed here: $possible_formats. Format 'ddmonyyyy' should look like '25dec2020' and format yyyym00 should look like '2020m12'.")
    end 

    return output 
end

function parse_date_to_string(date, date_format::AbstractString)
    
    # This function is basically a wrapper for parse_date_to_string()
    # except this adds a bit more functionality so that it can 
    # take date strings in Stata formats (e.g. 25dec2020 or 2020m12) and return those as strings
    if date_format == "ddmonyyyy"
        month_dict = Dict("01" => "jan", "02" => "feb", "03" => "mar", "04" => "apr", "05" => "may", 
        "06" => "jun", "07" => "jul", "08" => "aug", "09" => "sep", "10" => "oct", 
        "11" => "nov", "12" => "dec")
        vectorized_date_object = split(string(date), "-")
        return "$(vectorized_date_object[3])$(month_dict[vectorized_date_object[2]])$(vectorized_date_object[1])"        
    elseif date_format == "yyyym00"
        month_dict = Dict("01" => "m1", "02" => "m2", "03" => "m3", "04" => "m4", "05" => "m5", 
        "06" => "m6", "07" => "m7", "08" => "m8", "09" => "m9", "10" => "m10", 
        "11" => "m11", "12" => "m12")
        vectorized_date_object = split(string(date), "-")
        return "$(vectorized_date_object[1])$(month_dict[vectorized_date_object[2]])"
    else
        return Dates.format(date, date_format)
    end 
end

function parse_freq(period_str::AbstractString)
    parts = split(period_str)
    value = parse(Int, parts[1])
    period_type = parts[2]
    
    if period_type in ["week", "weeks"]
        return Week(value)
    elseif period_type in ["day", "days"]
        return Day(value)
    elseif period_type in ["month", "months"]
        return Month(value)
    elseif period_type in ["year", "years"]
        return Year(value)
    else
        throw(ArgumentError("Unsupported period type: $period_type, try day(s), week(s), month(s), or year(s)."))
    end
end

function match_dates(empty_diff_dates::Vector{Date}, silo_unmatched_dates::Vector{Date}, freq)
    
    # Initialize empty dictionary
    date_dict = Dict{Date, Date}()
    
    # First ensure that there are enough dates to match to the silo_unmatched_dates
    min_diff_date = minimum(empty_diff_dates)
    max_diff_date = maximum(empty_diff_dates)
    backward_dates = Date[]
    date = min_diff_date - freq
    while date >= minimum(silo_unmatched_dates) - freq
        push!(backward_dates, date)
        date -= freq
    end
    forward_dates = Date[]
    date = max_diff_date + freq
    while date <= maximum(silo_unmatched_dates) + freq
        push!(forward_dates, date)
        date += freq
    end
    empty_diff_dates = sort(vcat(backward_dates, empty_diff_dates, forward_dates))

    for target in silo_unmatched_dates
        if target >= minimum(empty_diff_dates)
            date_dict[target] = maximum([d for d in empty_diff_dates if d <= target])
        end
    end 
    return date_dict

end

end 
