module Undid

# All of these are part of base Julia with the exception of DataFrames
using Statistics
using LinearAlgebra
using DataFrames
using DelimitedFiles
using Dates


export create_init_csv, create_diff_df, run_stage_one, # Stage 1 Functions
fill_diff_df, create_trends_df,  run_stage_two, # Stage 2 Functions
combine_diff_data, combine_trends_data, calculate_agg_att_df, run_stage_three, # Stage 3 Functions
read_csv_data, compute_covariance_matrix, save_as_csv, parse_string_to_date, parse_date_to_string, parse_freq, match_dates, combine_fuzzy_data # Misc Functions

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

function create_diff_df(csv::AbstractString; covariates = false, date_format = false, freq = false, freq_multiplier = false, confine_matching::Bool = true, return_filepath::Bool = false)

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
        header = ["silo_name", "gvar", "treat", "diff_times", "(g;t)"]
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
        diff_df[!,"(g;t)"] = map(x -> join(x, ";"), diff_df[!,"(g;t)"])

    # For common treatment time
    elseif length(unique(df.treatment_time)) <= 2

        # Set up the skeleton of the full empty_diff_df
        header = ["silo_name", "treat", "common_treatment_time", "start_time", "end_time"]
        column_types = [Any, Int, String, Date, Date]
        diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])
        common_treatment_time = parse_date_to_string(df[df[!, "treatment_time"] .!== 0, "treatment_time"][1], date_format)

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

    # And if confine_matching == false then enable fuzzy matching of dates
    if confine_matching == false 
        diff_df.fuzzy_matching = Vector{String}(fill("fuzzy", nrow(diff_df))) 
    end

    # Finally, this if loop adds the necessary rows needed to carry out the RI inference later
    if "gvar" in DataFrames.names(diff_df) && length(unique(diff_df.gvar)) == 2
        diff_df.RI_inference = fill(0, nrow(diff_df))
        diff_df.treat = convert(Vector{Any}, diff_df.treat)
        for silo in unique(diff_df[diff_df.treat .== 1, "silo_name"])
            diff_times_to_add = setdiff(unique(diff_df.diff_times),diff_df[diff_df.silo_name .== silo ,:].diff_times)
            for time in diff_times_to_add
                filtered_df = unique(filter(row -> row.diff_times == time, diff_df)[:, Not([:silo_name, :RI_inference, :treat])])
                gvar = filtered_df.gvar[1]
                filtered_df = filtered_df[:, Not(:gvar)]
                insertcols!(filtered_df, 1, :silo_name => [silo])
                insertcols!(filtered_df, 2, :gvar => [gvar])
                insertcols!(filtered_df, 3, :treat => ["RI"])
                filtered_df[!, :RI_inference] = [1]
                append!(diff_df, filtered_df)
            end 

        end 
    end

    # Save as empty_diff_df.csv
    if return_filepath == false
        save_as_csv("empty_diff_df.csv", diff_df, "df")
        return diff_df 
    elseif return_filepath == true
        filepath = save_as_csv("empty_diff_df.csv", diff_df, "df", false)
        return filepath
    end
    
end 

function run_stage_one(names, start_times, end_times, treatment_times; covariates = false, date_format = false, freq = false, freq_multiplier = false, confine_matching::Bool = true)
    csv = create_init_csv(names, start_times, end_times, treatment_times, printmessage = true)
    return create_diff_df(csv, covariates = covariates, date_format = date_format, freq = freq, freq_multiplier = freq_multiplier, confine_matching = confine_matching)
end

### Stage 2 Functions ### 
function fill_diff_df(silo_name::AbstractString, empty_diff_df::DataFrame, silo_data::DataFrame; treatment_time = false, return_filepath::Bool = false, consider_covariates::Bool = true)
       
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
        empty_diff_df[!, "(g;t)"] = [join((parse_date_to_string(date1, empty_diff_df_date_format), parse_date_to_string(date2, empty_diff_df_date_format)), ";") for (date1, date2) in empty_diff_df[!, "(g;t)"]]
        empty_diff_df.diff_times = [join((parse_date_to_string(date1, empty_diff_df_date_format), parse_date_to_string(date2, empty_diff_df_date_format)), ";") for (date1, date2) in empty_diff_df.diff_times]

    end  

    # Save as csv and restructure covariates so that they are delimited by ; makes it easier to read back in from csv
    silo_name = string(silo_name)
    empty_diff_df.covariates = fill(join(empty_diff_df.covariates[1], ";"), nrow(empty_diff_df))
      
    
    # Save as empty_diff_df.csv
    if return_filepath == false
        save_as_csv("filled_diff_df_$silo_name.csv", empty_diff_df, "df")
        return empty_diff_df   
    elseif return_filepath == true
        filepath = save_as_csv("filled_diff_df_$silo_name.csv", empty_diff_df, "df", false)
        return filepath
    end
       
end 

function create_trends_df(silo_name::AbstractString, silo_data::DataFrame, freq; covariates::Vector{String} = ["none"], treatment_time = missing, date_format::AbstractString, return_filepath::Bool = false, consider_covariates::Bool = true)    

    # Define column headers
    header = ["silo_name", "treatment_time", "time", "mean_outcome", "mean_outcome_residualized", "covariates", "date_format", "freq"]
    trends_df = DataFrame(Symbol.(header) .=> [[] for column in header])
    
    # Allows for the option of ignoring covariates, even if specified initially in stage one.
    if consider_covariates == false
        covariates = ["none"]
    end

    # Push means and time to data
    if covariates == ["none"]
        for x in minimum(silo_data[!,"time"]):freq:maximum(silo_data[!,"time"])
            push!(trends_df, [silo_name, string(treatment_time), parse_date_to_string(x, date_format), string(mean(silo_data[silo_data[!, "time"] .== x, "outcome"])), "n/a", ["none"], string(date_format), string(freq)])
        end
    else        
        for x in minimum(silo_data[!,"time"]):freq:maximum(silo_data[!,"time"])            
            
            silo_subset = silo_data[silo_data[!, "time"] .== x,:]
            

            X = Matrix(silo_subset[:, covariates])
            Y = convert(Vector{Float64},silo_subset.outcome)
           
            beta_hat = X\Y 
            Y_hat = X * beta_hat
            residuals = Y - Y_hat

            push!(trends_df, [silo_name, string(treatment_time), parse_date_to_string(x, date_format), string(mean(silo_data[silo_data[!, "time"] .== x, "outcome"])), string(mean(residuals)), covariates, string(date_format), string(freq)])
        end
    end
    
    trends_df.covariates = fill(join(trends_df.covariates[1], ";"), nrow(trends_df))
    
    # Save as empty_diff_df.csv
    if return_filepath == false
        save_as_csv("trends_data_$(silo_name).csv", trends_df, "df")  
        return trends_df  
    elseif return_filepath == true
        filepath = save_as_csv("trends_data_$(silo_name).csv", trends_df, "df", false)
        return filepath
    end

end

function run_stage_two(filepath_to_empty_diff_df::AbstractString, silo_name::AbstractString, silo_data::DataFrame, time_column::AbstractString, outcome_column::AbstractString, date_format_local::AbstractString; renaming_dictionary = false, return_filepath::Bool = false, consider_covariates::Bool = true)
    
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
        gvar = empty_diff_df[empty_diff_df.treat .!== "RI", "gvar"][1]
        if empty_diff_df[empty_diff_df.treat .!== "RI", "treat"][1] == 1
            ever_treated = true
        else
            ever_treated = false
        end
        if "fuzzy_matching" in DataFrames.names(empty_diff_df)
            confine_matching = false
        else
            confine_matching = true
        end 
        if length(silo_unmatched_dates) > 0 
            date_map = match_dates(empty_diff_dates, silo_unmatched_dates, gvar, silo_name, ever_treated, freq, confine_matching = confine_matching)
            silo_data.time = get.(Ref(date_map), silo_data.time, silo_data.time)
        end 
    end

    if return_filepath == false
        fill_diff_df(silo_name, empty_diff_df, silo_data, treatment_time = treatment_time, consider_covariates = consider_covariates)
    elseif return_filepath == true
        filepath_diff_df = fill_diff_df(silo_name, empty_diff_df, silo_data, treatment_time = treatment_time, return_filepath = true, consider_covariates = consider_covariates)
    end 

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
    
    if return_filepath == false
        create_trends_df(silo_name, silo_data, freq, covariates = covariates, treatment_time = treatment_time, date_format = empty_diff_df_date_format, consider_covariates = consider_covariates)
    elseif return_filepath == true
        filepath_trends_df = create_trends_df(silo_name, silo_data, freq, covariates = covariates, treatment_time = treatment_time, date_format = empty_diff_df_date_format, return_filepath = true, consider_covariates = consider_covariates)
    end 

    if return_filepath == true
        return filepath_diff_df, filepath_trends_df
    end 
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
                gt = data[index,"(g;t)"]
                println("Using a linear function to fill missing values of diff_estimate for silo: $silo (g,t): $gt")
                periods = sort(unique([x for x in data[(data.silo_name .== silo) .&& (data.gvar .== gvar),"(g;t)"] for x in x]))
                for i in 1:length(periods)
                    data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& ((getindex.(data."(g;t)", 2)) .== periods[i]), "local_period"] .= i
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

    # This block performs linear interpolation/extrapolation for diff_estimate_covariates if interpolation is set to true
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
                    gt = data[index,"(g;t)"]
                    println("Using a linear function to fill missing values of diff_estimate_covariates for silo: $silo (g,t): $gt")
                    periods = sort(unique([x for x in data[(data.silo_name .== silo) .&& (data.gvar .== gvar),"(g;t)"] for x in x]))
                    for i in 1:length(periods)
                        data[(data.silo_name .== silo) .&& (data.gvar .== gvar) .&& ((getindex.(data."(g;t)", 2)) .== periods[i]), "local_period"] .= i
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
                error("Set interpolation to false to \"linear_function\".")
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
            copy_df[!, "(g;t)"] = [join((parse_date_to_string(date1, copy_df_date_format), parse_date_to_string(date2, copy_df_date_format)), ";") for (date1, date2) in copy_df[!, "(g;t)"]]
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

function run_stage_three(dir_path::AbstractString; agg::AbstractString = "silo", covariates::Bool = false, save_all_csvs::Bool = false, interpolation = false)

    combined_diff_data = combine_diff_data(dir_path, save_csv = save_all_csvs, interpolation = interpolation) 
    

    # Generate all the necessary matrices to do the randomization inference
    if "RI_inference" in DataFrames.names(combined_diff_data)
        # Filter for just the data that is reflective of the silos true treatment status 
        diff_matrix_no_RI = combined_diff_data[combined_diff_data.RI_inference .== 0,:]

        # Grab the silo names and their associated gvar(s) and treatment status 
        silos = unique(diff_matrix_no_RI.silo_name)
        gvar_treat = []
        for silo in silos
            push!(gvar_treat, (unique(diff_matrix_no_RI[diff_matrix_no_RI.silo_name .== silo, "gvar"]), unique(diff_matrix_no_RI[diff_matrix_no_RI.silo_name .== silo, "treat"])))
        end
        
        # Create a df indicating the original gvar(s) and treatment status of silos
        # And use circshift to shuffle the gvar(s) and treatment status until a full loop is completed
        df = DataFrame(silo_name = silos, D = gvar_treat)
        for i in 1:nrow(df)
            df[!, Symbol("D_$(i)")] = circshift(df.D, i)
        end

        # Create a dictionary of empty dataframes for each rotation in the circshift
        RI_matrices = Dict{Int, DataFrame}()
        for i in 1:nrow(df)
            empty_matrix = combined_diff_data[combined_diff_data.gvar .== false, :]  
            RI_matrices[i] = empty_matrix
        end

        # Fill each one of those dataframes with the gvar(s) and treatment status based on the circshift assignment
        for i in 1:length(RI_matrices)
            D_star = Symbol("D_$i")
            for silo in silos
                
                if df[df.silo_name .== silo, D_star][1][2][1] == 0
                    rows_to_add = combined_diff_data[combined_diff_data.silo_name .== silo .&& in.(combined_diff_data.gvar, Ref(df[df.silo_name .== silo, D_star][1][1])),:]
                    RI_matrices[i] = vcat(RI_matrices[i], rows_to_add)
                    
                    RI_matrices[i][RI_matrices[i].silo_name .== silo, "treat"] .= 0
                else
                    rows_to_add = combined_diff_data[combined_diff_data.silo_name .== silo .&& in.(combined_diff_data.gvar, Ref(df[df.silo_name .== silo, D_star][1][1])),:]
                    RI_matrices[i] = vcat(RI_matrices[i], rows_to_add)
                    
                    RI_matrices[i][RI_matrices[i].silo_name .== silo, "treat"] .= 1
                end 
            end 
        end 
        
        # Calculate att by specified aggregation method 
        results = calculate_agg_att_df(diff_matrix_no_RI; agg = agg, covariates = covariates, save_all_csvs = save_all_csvs, printinfo = true)  

        # Compute p-value based on randomization inference and add to results as a column
        original_ATT = results.agg_ATT[1]
        RI_ATT = []
        for i in 1:length(RI_matrices)
            RI_ATT_single = calculate_agg_att_df(RI_matrices[i]; agg = agg, covariates = covariates, save_all_csvs = false, printinfo = false).agg_ATT[1]
            push!(RI_ATT, RI_ATT_single)
        end 

        p_value_RI = 0
        for ATT in RI_ATT
            p_value_RI += Int(abs(ATT) >= abs(original_ATT))
        end 
        p_value_RI = p_value_RI/ length(RI_ATT)
        results.p_value_RI = fill(p_value_RI, nrow(results))
        save_as_csv("UNDID_results.csv", results, "df", true)
        return results
    end 

    results = calculate_agg_att_df(combined_diff_data; agg = agg, covariates = covariates, save_all_csvs = save_all_csvs) 
    save_as_csv("UNDID_results.csv", results, "df", true)
    return results  
end



function calculate_agg_att_df(combined_diff_data::DataFrame; agg::AbstractString = "silo", covariates::Bool = false, save_all_csvs::Bool = false, printinfo::Bool = true)

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


    if "common_treatment_time" in DataFrames.names(combined_diff_data)
        println("Calcualting aggregate ATT with covariates set to $covariates.")
        ATT = (hcat(fill(1.0, length(combined_diff_data.treat)), combined_diff_data.treat) \ combined_diff_data.y)[2]
        treatment_time = combined_diff_data.common_treatment_time[1]
        results = DataFrame(treatment_time = treatment_time, ATT = ATT)
        # Compute jackknife SE if there are at least 2 controls and 2 treatment silos
        if sum(combined_diff_data.treat .== 1) >= 2 && sum(combined_diff_data.treat .== 0) >= 2
            jackknives_common = []
            for silo in combined_diff_data.silo_name
                subset = combined_diff_data[combined_diff_data.silo_name .!= silo,:]
                push!(jackknives_common, (hcat(fill(1.0, length(subset.treat)), subset.treat) \ subset.y)[2])
            end 
            jackknife_SE = sqrt(sum((jackknives_common .- ATT).^2) / length(jackknives_common))
            results.jackknife_SE = [jackknife_SE]

        elseif sum(combined_diff_data.treat .== 1) == 1 && sum(combined_diff_data.treat .== 0) == 1
            if covariates == false 
                SE = sqrt(sum(convert(Vector{Float64}, combined_diff_data.diff_var)))
            elseif covariates == true 
                SE = sqrt(sum(convert(Vector{Float64}, combined_diff_data.diff_var_covariates)))
            end 
            results.SE = [SE]
         
        # Do randomization inference in two cases: 1 treated and >1 controls, >1 treated and 1 control 
        elseif (sum(combined_diff_data.treat .== 0) == 1 && sum(combined_diff_data.treat .== 1) >= 2) || (sum(combined_diff_data.treat .== 0) >= 2 && sum(combined_diff_data.treat .== 1) == 1)
            if sum(combined_diff_data.treat .== 0) == 1
                unique_value = 0
                non_unique = 1
            elseif sum(combined_diff_data.treat .== 1) == 1
                unique_value = 1
                non_unique = 0
            end
            beta_r = []
            for i in 1:nrow(combined_diff_data)
                combined_diff_data.treat .= non_unique
                combined_diff_data.treat[i] = unique_value
            
                beta = (hcat(fill(1.0, length(combined_diff_data.treat)), combined_diff_data.treat) \ combined_diff_data.y)[2]
                push!(beta_r, beta)
            end
        
            p_value_RI = 0
            for beta in beta_r 
                p_value_RI += Int(abs(beta) >= abs(ATT))
            end 
            p_value_RI = p_value_RI/ length(beta_r)
            results.p_value_RI = [p_value_RI]
        end

                
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
            for silo in treated_silos
                treated = combined_diff_data[combined_diff_data.silo_name .== silo,:]
                control = combined_diff_data[combined_diff_data.treat .== 0,:]            
                diff_times_treated = treated.diff_times
                control = filter(row -> row.diff_times in diff_times_treated, control)            
                local_df = vcat(control, treated)
                beta_hat = hcat(fill(1.0, length(local_df.treat)), local_df.treat) \ local_df.y
                push!(silos_att, beta_hat[2])                
            end
            jackknives_silo = []
            for i in 1:length(silos_att)
                # Exclude the i-th element and calculate the mean of the remaining elements
                mean_excluding_i = mean([silos_att[1:i-1]..., silos_att[i+1:end]...])
                # Push the result to the jackknives vector
                push!(jackknives_silo, mean_excluding_i)
            end
            ATT_silo = mean(silos_att)
            jackknife_SE = sqrt(sum((jackknives_silo .-ATT_silo).^2) / length(jackknives_silo))
            results = DataFrame(silos = unique(combined_diff_data[combined_diff_data.treat .== 1, "silo_name"]), ATT_s = silos_att, agg_ATT = fill(ATT_silo, length(silos_att)), jackknife_SE = fill(jackknife_SE, length(silos_att)))          
        elseif agg == "gt" || agg == "g"            
            ATT_vec = []
            gt_vec = []
            for gt in unique(combined_diff_data[!, "(g;t)"])                
                subset = combined_diff_data[(in.(gt[1], combined_diff_data[!, "(g;t)"])) .&& (in.(gt[2], combined_diff_data[!, "(g;t)"])),:]
                beta_hat = hcat(fill(1.0, length(subset.treat)),subset.treat) \ subset.y
                push!(ATT_vec, beta_hat[2])                
                push!(gt_vec, gt)
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
                jackknife_SE = sqrt(sum((jackknives_gt .-ATT_gt).^2) / length(jackknives_gt))
                results = DataFrame(gt = gt_vec, ATT_gt = ATT_vec, agg_ATT = fill(ATT_gt, length(ATT_vec)), jackknife_SE = fill(jackknife_SE, length(ATT_vec)))
                results.gt = [join((parse_date_to_string(date1, combined_diff_data.date_format[1]), parse_date_to_string(date2, combined_diff_data.date_format[1])), ";") for (date1, date2) in results[!, "gt"]]
            end 
            if agg == "g"                
                gt_g = unique(combined_diff_data.gvar)
                counters = Dict{Date, Int}()
                for g in gt_g
                    # Initialize g counter if it doesn't exist yet
                    if !haskey(counters, g)
                        counters[g] = 0
                    end
                    # Add counts of g
                    for gt in gt_vec
                        if g == gt[1]
                            counters[g] += 1
                        end 
                    end
                end                 
                ATT_gt_df = DataFrame(ATT_gt = ATT_vec, gt = gt_vec)
                g_vec = []
                for g in ATT_gt_df.gt
                    push!(g_vec, g[1])
                end
                ATT_gt_df.g = g_vec
                g_weight = []
                for g in ATT_gt_df.g
                    push!(g_weight, 1/counters[g])
                end
                ATT_gt_df.g_weight = g_weight
                ATT_gt_df.weighted_ATT = ATT_gt_df.g_weight .* ATT_gt_df.ATT_gt
                ATT_by_gvar_weighted = []
                for g in unique(ATT_gt_df.g)
                    push!(ATT_by_gvar_weighted, sum(ATT_gt_df[ATT_gt_df.g .== g,"weighted_ATT"]))
                end
                # Compute the jackknife means
                jackknives_g = []
                for i in 1:length(ATT_by_gvar_weighted)
                    # Exclude the i-th element and calculate the mean of the remaining elements
                    mean_excluding_i = mean([ATT_by_gvar_weighted[1:i-1]..., ATT_by_gvar_weighted[i+1:end]...])
                    # Push the result to the jackknives vector
                    push!(jackknives_g, mean_excluding_i)
                end
                ATT_g = mean(ATT_by_gvar_weighted)
                jackknife_SE = sqrt(sum((jackknives_g .-ATT_g).^2) / length(jackknives_g))
                results = DataFrame(g = parse_date_to_string.(unique(ATT_gt_df.g), combined_diff_data.date_format[1]), ATT_g = ATT_by_gvar_weighted, agg_ATT = fill(ATT_g, length(ATT_by_gvar_weighted)), jackknife_SE = fill(jackknife_SE, length(ATT_by_gvar_weighted)))
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
        df_readin[!, "(g;t)"] = split.(df_readin[!, "(g;t)"], ";")
        transform!(df_readin, :"(g;t)" => ByRow(time -> (parse_string_to_date(time[1], df_readin_date_format, possible_formats_UNDID, month_map_UNDID), parse_string_to_date(time[2], df_readin_date_format, possible_formats_UNDID, month_map_UNDID))) => :"(g;t)")
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
    if date in ["control", "n/a", "never", "not treated", "not_treated", ""]
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

function match_dates(empty_diff_dates::Vector{Date}, silo_unmatched_dates::Vector{Date}, gvar::Date, silo_name, ever_treated::Bool, freq; confine_matching::Bool = true)
    
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

    # If confine matching is true then any date in silo_date gets matched to the most recently passed date in empty_diff_dates
    if confine_matching == true
        println("Matching local silo dates to the most recently passed date in empty_diff_df.")
        for target in silo_unmatched_dates
            if target >= minimum(empty_diff_dates)
                date_dict[target] = maximum([d for d in empty_diff_dates if d <= target])
            end
        end 
        return date_dict
    end 

    println("Fuzzy matching dates!")
    # The first two if statements ensure that we are only matching one date prior the minimum date 
    # in the empty_diff_dates and only one date post the maximum date in the empty_diff_dates.
    # Otherwise, every date in silo_unmatched_dates prior to the minimum date in empty_diff_dates would be matched to that minimum date.
    if maximum(silo_unmatched_dates) > maximum(empty_diff_dates)
        upper_bound = minimum([d for d in silo_unmatched_dates if d > maximum(empty_diff_dates)])
        filter!(d -> d <= upper_bound, silo_unmatched_dates)
    end
    
    if minimum(silo_unmatched_dates) < minimum(empty_diff_dates)
        lower_bound = maximum([d for d in silo_unmatched_dates if d < minimum(empty_diff_dates)])
        filter!(d -> d >= lower_bound, silo_unmatched_dates)
    end 

    for target in silo_unmatched_dates
        # Determine the valid dates based on gvar relative to the target date
        # Of course, this only matters if the silo is ever treated. Otherwise we don't care about it crossing the untreated/treated
        # time barrier, as it would be always be untreated
        if ever_treated == true
            if gvar > target
                valid_dates = filter(d -> d < gvar, empty_diff_dates) # Only allow dates up to 1 day behind gvar and not after gvar
            else
                valid_dates = filter(d -> d >= gvar, empty_diff_dates)  # Only allow dates that are no earlier than gvar
            end
        else
            valid_dates = empty_diff_dates
        end
        # Calculate absolute differences between each date and target date
        differences = abs.(valid_dates .- target)

        # Grab index of the smallest difference
        closest_index = argmin(differences)

        # Insert into dictionary
        date_dict[target] = valid_dates[closest_index]
    end 

    # This block is to keep track of which dates got fuzzy matched in which silos
    original_dates = parse_date_to_string.(keys(date_dict), "dd/mm/yyyy")
    matched_to_dates = parse_date_to_string.(values(date_dict), "dd/mm/yyyy")
    df = DataFrame(silo_name = fill(silo_name, length(date_dict)), original_date =  original_dates, matched_to = matched_to_dates)
    save_as_csv("fuzzy_dates_$(silo_name).csv", df ,"df", true)
    
    return date_dict
end

function combine_fuzzy_data(dir_path::String; save_csv::Bool = false)
    # Collects all the filled_diff_df_... csv files
    files = readdir(dir_path)
    matching_files = filter(file -> startswith(file, "fuzzy_dates_") && endswith(file, ".csv"), files)

    # Uses the read_csv_data function to read in the csv's and appends them all together
    data = read_csv_data("$dir_path\\$(matching_files[1])")
    for i in 2:length(matching_files)
        data = vcat(data, read_csv_data("$dir_path\\$(matching_files[i])"))
    end 

    if save_csv == true
        save_as_csv("combined_fuzzy_dates.csv", data, "df")
    end

    # This is for flexible date handling for the parse_string_to_date function
    possible_formats_UNDID = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm", "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy", "mm-dd-yyyy", "mmddyyyy",
    "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy", "ddmonyyyy", "yyyym00"]
    month_map_UNDID = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04", "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08", "sep" => "09", "oct" => "10",
    "nov" => "11", "dec" => "12")

    data.original_date = parse_string_to_date.(data.original_date, "dd/mm/yyyy", Ref(possible_formats_UNDID), Ref(month_map_UNDID))
    data.matched_to = parse_string_to_date.(data.matched_to, "dd/mm/yyyy", Ref(possible_formats_UNDID), Ref(month_map_UNDID))

    return data
end 

end 
