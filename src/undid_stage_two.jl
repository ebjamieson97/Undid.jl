function undid_stage_two(filepath_to_empty_diff_df::AbstractString,
                         silo_name::AbstractString,
                         silo_data::DataFrame,
                         time_column::AbstractString,
                         outcome_column::AbstractString, 
                         date_format_local::AbstractString; 
                         renaming_dictionary = false, 
                         consider_covariates::Bool = true,
                         anonymize_weights::Bool = false,
                         anonymize_size::Number = 5)
    
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

    filepath_and_diff_df = fill_diff_df(silo_name, empty_diff_df, silo_data,
                                       treatment_time = treatment_time,
                                       consider_covariates = consider_covariates,
                                       anonymize_weights = anonymize_weights,
                                       anonymize_size = anonymize_size)

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

function fill_diff_df(silo_name::AbstractString,
                      empty_diff_df::DataFrame,
                      silo_data::DataFrame;
                      treatment_time = false,
                      consider_covariates::Bool = true,
                      anonymize_weights::Bool = false,
                      anonymize_size::Number = 5)
       
    # Grab date format from empty_diff_df
    empty_diff_df_date_format = empty_diff_df.date_format[1]

    # Assign weights and throw error commmon to both staggered adoption and common adoption
    weights = empty_diff_df.weights[1]
    if in(weights, ["diff", "att", "both"]) && !anonymize_weights
        @warn "Setting weights to \"diff\", \"both\", or \"att\" will output a CSV file with counts of obs. for each difference calculation.\n 
        Consider setting `anonymize_weights = true` to round counts to the nearest $anonymize_size."
    end 

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
        if in(weights, ["none", "diff", "att", "both", "standard"])
            if weights == "standard"
                @warn "\"standard\" is a deprecated weighting option!"
                n = sum(silo_data_copy.time) / nrow(silo_data_copy)
                empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "n"] .= n
            elseif in(weights, ["diff", "both", "att"])
                n = nrow(silo_data_copy)
                if anonymize_weights
                    n = max(anonymize_size, (anonymize_size * round(Int, n / anonymize_size)))
                end
                empty_diff_df[empty_diff_df[!, "silo_name"] .== silo_name, "n"] .= n
            end 
        else
            error("Please select a valid weighting method. Options include:\"none\", \"diff\", \"att\", \"both\".")
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
            cov_beta_hat = compute_covariance_matrix(X, resid)
            
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
            cov_beta_hat = compute_covariance_matrix(X, resid, diff_times = diff_combo)

            mask = (empty_diff_df[!, "silo_name"] .== silo_name) .&&
                   (in.(diff_combo[1], empty_diff_df[!, "diff_times"])) .&& 
                   (in.(diff_combo[2], empty_diff_df[!, "diff_times"]))
            
            empty_diff_df[mask, "diff_estimate"] .= beta_hat[2]
            empty_diff_df[mask, "diff_var"] .= cov_beta_hat[2,2]
            
            # Assign weights
            if in(weights, ["none", "diff", "att", "both"])
                if in(weights, ["diff", "both"])
                    n = length(Y)
                    if anonymize_weights
                        n = max(anonymize_size, (anonymize_size * round(Int, n / anonymize_size)))
                    end
                    empty_diff_df[mask, "n"] .= n
                end 
                if in(weights, ["att", "both"])
                    n_t = sum((silo_subset.time .== 1))
                    if anonymize_weights
                        n_t = max(anonymize_size, (anonymize_size * round(Int, n_t / anonymize_size)))
                    end 
                    empty_diff_df[mask, "n_t"] .= n_t
                end 
            else
                error("Please select a valid weighting method. Options include:\"none\", \"diff\", \"att\", \"both\".")
            end

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
                cov_beta_hat = compute_covariance_matrix(X, resid, diff_times = diff_combo, covariates = empty_diff_df[!, "covariates"][1])              
                empty_diff_df[mask, "diff_estimate_covariates"] .= beta_hat[2]
                empty_diff_df[mask, "diff_var_covariates"] .= cov_beta_hat[2,2]
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
