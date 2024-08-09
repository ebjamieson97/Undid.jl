module UNDID

using Statistics
using LinearAlgebra
using DataFrames
using DelimitedFiles
using Dates

export create_init_csv, create_diff_df, fill_diff_df, read_csv_data, combine_diff_data, create_ATT_by_gt_by_silo, compute_ATT_staggered,
        create_trends_df, run_stage_two, combine_trends_data, compute_ATT_common, compute_covariance_matrix, save_as_csv, run_stage_one,
        run_stage_three, parse_string_to_date, parse_freq, parse_date_to_string

### Stage 1 Functions ###
function create_init_csv(names=[], start_times=[], end_times=[], treatment_times=[]; covariates = false)
    
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
    save_as_csv("init.csv", data, "array_of_arrays")   
end 

function create_diff_df(csv; covariates = false, date_format = false, freq = false, freq_multiplier = false)

    # Ensure that date_format is manually selected
    if date_format == false
        error("Please ensure the start_times, end_times and treatment_times are formatted identically and set date_format equal to one of the following:[\"yyyy/mm/dd\", \"yyyy-mm-dd\", \"yyyymmdd\", \"yyyy/dd/mm\", \"yyyy-dd-mm\", \"yyyyddmm\", \"dd/mm/yyyy\", \"dd-mm-yyyy\", \"ddmmyyyy\", \"mm/dd/yyyy\", \"mm-dd-yyyy\", \"mmddyyyy\",
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
        error("Please enter a frequency of data: 'daily', 'monthly', 'weekly' or 'yearly'")
    end 

    # Allow for weird frequencies like every 4 days or every 7 months etc
    if freq_multiplier == false
        # do nothing
    elseif typeof(freq_multiplier) <: Integer && freq_multiplier !== 0
        freq = freq*freq_multiplier
    else
        error("freq_multiplier must be an integer or set to false.")
    end

    # Read in the init.csv and make it a dataframe
    data = readdlm(csv, ',')
    header, rows  = data[1, :], data[2:end, :]
    df = DataFrame(rows, Symbol.(header))

    # This is for flexible date handling for the parse_string_to_date function
    possible_formats_UNDID = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm", "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy", "mm-dd-yyyy", "mmddyyyy",
    "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy", "ddmonyyyy", "yyyym00"]
    month_map_UNDID = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04", "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08", "sep" => "09", "oct" => "10",
    "nov" => "11", "dec" => "12")
    
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
            if df[df[!, "silo_name"] .== group, "treatment_time"][1] !== Date(0)
                # These are the years that matter for treated groups                
                combinations = [(t1, t2) for t1 in times, t2 in times if t1 > t2 && Dates.value(t1 - t2) >= Dates.value(freq) && Dates.value(df[df[!, "silo_name"] .== group, "treatment_time"][1] - freq) == Dates.value(t2)]
                treat = 1                
            elseif df[df[!, "silo_name"] .== group, "treatment_time"][1] == Date(0)
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
        header = ["silo_name", "treat", "common_treatment_time"]
        column_types = [Any, Int, String]
        diff_df = DataFrame([Symbol(header[i]) => Vector{column_types[i]}() for i in 1:length(header)])
        common_treatment_time = Dates.format(df[df[!, "treatment_time"] .!== 0, "treatment_time"][1], date_format)

        for silo_name in unique(df.silo_name)
            if df[df[!, "silo_name"] .== silo_name, "treatment_time"][1] !== Date(0)
                treat = 1
            elseif df[df[!, "silo_name"] .== silo_name, "treatment_time"][1] == Date(0)
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

    # Save as empty_diff_df.csv
    save_as_csv("empty_diff_df.csv", diff_df, "df")
    
    return diff_df
    
end 

function run_stage_one(names, start_times, end_times, treatment_times; covariates = false, date_format = false, freq = false, freq_multiplier = false)
    csv = create_init_csv(names, start_times, end_times, treatment_times)
    return create_diff_df(csv, covariates = covariates, date_format = date_format, freq = freq, freq_multiplier = freq_multiplier)
end

### Stage 2 Functions ### 
function read_csv_data(filepath_to_csv_data)
    # Reads in the empty_diff_df.csv
    data = readdlm(filepath_to_csv_data, ',')
    header = data[1, :]    
    rows = data[2:end, :]
    df_readin = DataFrame(rows, Symbol.(header))    
    df_readin.silo_name = string.(df_readin.silo_name)
    if "covariates" in DataFrames.names(df_readin)
        transform!(df_readin, :covariates => ByRow(s -> [String(sub) for sub in split(replace(replace(strip(s, ['[' , ']', '"']), "\"" => ""), " " => ""), ";")]) => :covariates)
    end
    if "diff_times" in DataFrames.names(df_readin)        
        df_readin_date_format = df_readin.date_format[1]

        df_readin.diff_times = split.(df_readin.diff_times, ";")        
        transform!(df_readin, :diff_times => ByRow(time -> (Date(time[1], df_readin_date_format), Date(time[2], df_readin_date_format))) => :diff_times)

        df_readin.gvar = Date.(string.(df_readin.gvar), df_readin_date_format)

        df_readin[!, "(g;t)"] = split.(df_readin[!, "(g;t)"], ";")
        transform!(df_readin, :"(g;t)" => ByRow(time -> (Date(time[1], df_readin_date_format), Date(time[2], df_readin_date_format))) => :"(g;t)")
    end 
    if "time" in DataFrames.names(df_readin)
        # This is for flexible date handling for the parse_string_to_date function
        possible_formats_UNDID = ["yyyy/mm/dd", "yyyy-mm-dd", "yyyymmdd", "yyyy/dd/mm", "yyyy-dd-mm", "yyyyddmm", "dd/mm/yyyy", "dd-mm-yyyy", "ddmmyyyy", "mm/dd/yyyy", "mm-dd-yyyy", "mmddyyyy",
        "mm/yyyy", "mm-yyyy", "mmyyyy", "yyyy", "ddmonyyyy", "yyyym00"]
        month_map_UNDID = Dict("jan" => "01", "feb" => "02", "mar" => "03", "apr" => "04", "may" => "05", "jun" => "06", "jul" => "07", "aug" => "08", "sep" => "09", "oct" => "10",
        "nov" => "11", "dec" => "12")
        df_readin.time = parse_string_to_date.(string.(df_readin.time), df_readin.date_format[1], Ref(possible_formats_UNDID), Ref(month_map_UNDID))

        # Now based on Date object we can created a column "period" that will translate each date to an associated integer period
        dates_vector  = sort(unique(df_readin.time))
        date_dict = Dict(zip(dates_vector, 1:length(dates_vector)))
        df_readin.period = [date_dict[date] for date in df_readin.time]

        # Likewise we define the treatment_time as a period as well
        df_readin[df_readin.treatment_time .!= "control", "treatment_time"] = Date.(string.(df_readin[df_readin.treatment_time .!= "control",:].treatment_time), "yyyy")
        df_readin.treatment_period = [x isa Date ? date_dict[x] : missing for x in df_readin.treatment_time]
    end 
    return df_readin
end

function fill_diff_df(silo_name, empty_diff_df, silo_data; treatment_time = false)
       
    empty_diff_df_date_format = empty_diff_df.date_format[1]
    # Compute diff (gamma) estimates and standard errors in the case of common treatment times across silos
    if "common_treatment_time" in DataFrames.names(empty_diff_df)
        if treatment_time == false
            error("Please specify a common treatment time")
        end

        # Use treatment time to construct dummy varaible indicating obs after treatment time
        silo_data_copy = copy(silo_data)
        silo_data_copy.time = ifelse.(silo_data_copy.time .>= treatment_time, 1, 0)

        # Construct X and Y for regression
        X = convert(Matrix{Float64},hcat(fill(1,nrow(silo_data_copy)), silo_data_copy.time))
        Y = convert(Vector{Float64},silo_data_copy.outcome)

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

            if cov_on == 1 && i == 1
                try
                    X = hcat(X, Matrix(silo_data[:, empty_diff_df.covariates[1]]))
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

        # Return date objects to strings
        empty_diff_df.gvar = Dates.format.(empty_diff_df.gvar, empty_diff_df_date_format)
        empty_diff_df[!, "(g;t)"] = [join((Dates.format(date1, empty_diff_df_date_format), Dates.format(date2, empty_diff_df_date_format)), ";") for (date1, date2) in empty_diff_df[!, "(g;t)"]]
        empty_diff_df.diff_times = [join((Dates.format(date1, empty_diff_df_date_format), Dates.format(date2, empty_diff_df_date_format)), ";") for (date1, date2) in empty_diff_df.diff_times]

    end  

    # Save as csv and restructure covariates so that they are delimited by ; makes it easier to read back in from csv
    silo_name = string(silo_name)
    empty_diff_df.covariates = fill(join(empty_diff_df.covariates[1], ";"), nrow(empty_diff_df))
    save_as_csv("filled_diff_df_$silo_name.csv", empty_diff_df, "df")     
    
    return empty_diff_df        
end 

function create_trends_df(silo_name, silo_data, freq; covariates = ["none"], treatment_time = missing, date_format)    

    # Define column headers
    header = ["silo_name", "treatment_time", "time", "mean_outcome", "mean_outcome_residualized", "covariates", "date_format", "freq"]
    trends_df = DataFrame(Symbol.(header) .=> [[] for column in header])
    
    # Push means and time to data
    if covariates == ["none"]
        for x in minimum(silo_data[!,"time"]):freq:maximum(silo_data[!,"time"])
            push!(trends_df, [silo_name, string(treatment_time), Dates.format(x, date_format), string(mean(silo_data[silo_data[!, "time"] .== x, "outcome"])), "n/a", ["none"], string(date_format), string(freq)])
        end
    else        
        for x in minimum(silo_data[!,"time"]):freq:maximum(silo_data[!,"time"])            
            
            silo_subset = silo_data[silo_data[!, "time"] .== x,:]

            X = Matrix(silo_subset[:, covariates])
            Y = convert(Vector{Float64},silo_subset.outcome)
           
            beta_hat = X\Y 
            Y_hat = X * beta_hat
            residuals = Y - Y_hat

            push!(trends_df, [silo_name, string(treatment_time), Dates.format(x, date_format), string(mean(silo_data[silo_data[!, "time"] .== x, "outcome"])), string(mean(residuals)), covariates, string(date_format), string(freq)])
        end
    end
    
    trends_df.covariates = fill(join(trends_df.covariates[1], ";"), nrow(trends_df))
    save_as_csv("trends_data_$(silo_name).csv", trends_df, "df")       
    return trends_df
end

function run_stage_two(filepath_to_empty_diff_df, silo_name, silo_data, time_column, outcome_column, date_format_local; renaming_dictionary = false)
    
    # Given a filepath to the empty_diff_df, the name of the local silo, and 
    # a dataframe of the local silo data, runs all the necessary Stage 2 functions
    empty_diff_df = read_csv_data(filepath_to_empty_diff_df)
    empty_diff_df.silo_name = string.(empty_diff_df.silo_name)
    covariates = empty_diff_df[!, "covariates"][1]
    
    # First doing some pre-processing and error checking

    # Check that the silo_name is in the empty_diff_df
    silo_name = string(silo_name)
    if silo_name in unique(empty_diff_df.silo_name)
        # do nothing
    else
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

    fill_diff_df(silo_name, empty_diff_df, silo_data, treatment_time = treatment_time)    

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
    
    freq = parse_freq(string(empty_diff_df.freq[1]))

    create_trends_df(silo_name, silo_data, freq, covariates = covariates, treatment_time = treatment_time, date_format = empty_diff_df_date_format)

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
        save_as_csv("combined_diff_data.csv", data, "df")        
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
        save_as_csv("combined_trends_data.csv", data, "df")        
    end

    # Returns df
    return data

end 

function create_ATT_by_gt_by_silo(combined_diff_data; save_csv = false, running_from_stage_3_wrapper = false)
    
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
                push!(gt_vec, combined_diff_data[(in.(time[1], combined_diff_data[!, "diff_times"])) .&& (in.(time[2], combined_diff_data[!, "diff_times"])) .&& (combined_diff_data[!, "treat"] .== 1) .&& (combined_diff_data[!, "silo_name"] .== treated_group), "(g;t)"][1])
                push!(ATT, attgt[1])
                push!(SE_ATT,  sqrt((se_control_diff + se_treated_diff)[1]))              
                push!(covariates, covariate)

            end 

        end 

    end 

    # Turn the filled preallocation vectors into a dataframe
    ATT_by_gt_by_silo_df = ATT_by_gt_by_silo_df = DataFrame("treated_silo" => treated_silo, "control_silo" => control_silo,
    "diff_times" => diff_times, "ATT" => ATT, "SE_ATT" => SE_ATT, "ATT_covariates" => ATT_covariates, "SE_ATT_covariates" => SE_ATT_covariates,
    "(g;t)" => gt_vec, "covariates" => covariates)

    # If save_as_csv == true then save as a csv as well
    if save_csv == true
        save_as_csv("ATT_by_gt_by_silo.csv", ATT_by_gt_by_silo_df, "df")        
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
    save_as_csv("UNDID_results.csv", ATTs_bysilo_df, "df")

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
        jackknife_SE_aggregate_ATT[1] = sqrt((1/length(jacknives))*jackknife_SE_aggregate_ATT[1])
        
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
    save_as_csv("UNDID_results.csv", undid_results, "df")


    # Return the dataframe
    return undid_results    
end

function run_stage_three(dir_path; save_all_csvs = false)

    combined_diff_data = combine_diff_data(dir_path, save_csv = save_all_csvs)    
    combined_trends_data = combine_trends_data(dir_path, save_csv = save_all_csvs)

    if "common_treatment_time" in DataFrames.names(combined_diff_data)
        results = compute_ATT_common(combined_diff_data)

        return combined_trends_data, combined_diff_data, results
    else
        ATT_by_gt_by_silo = create_ATT_by_gt_by_silo(combined_diff_data, save_csv = save_all_csvs, running_from_stage_3_wrapper = true)
        results = compute_ATT_staggered(ATT_by_gt_by_silo)

        return combined_trends_data, combined_diff_data, ATT_by_gt_by_silo, results
    end 

    

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

function save_as_csv(filename, data, datatype)

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

    println("$filename.csv saved to")
    formatted_path = replace(filepath, "\\" => "\\\\")
    println(formatted_path)

    return filepath
end 

function parse_string_to_date(date, date_format, possible_formats = false, month_map = false)
    
    if possible_formats == false || month_map == false
        error("Please pass a dictionary such as Dict(\"jan\" => \"01\", \"feb\" => \"02\", \"mar\" => \"03\", \"apr\" => \"04\", \"may\" => \"05\", \"jun\" => \"06\", \"jul\" => \"07\", \"aug\" => \"08\", \"sep\" => \"09\", \"oct\" => \"10\", \"nov\" => \"11\", \"dec\" => \"12\") as well as the possible formats [\"yyyy/mm/dd\", \"yyyy-mm-dd\", \"yyyymmdd\", \"yyyy/dd/mm\", \"yyyy-dd-mm\", \"yyyyddmm\", \"dd/mm/yyyy\", \"dd-mm-yyyy\", \"ddmmyyyy\", \"mm/dd/yyyy\", \"mm-dd-yyyy\", \"mmddyyyy\",
        \"mm/yyyy\", \"mm-yyyy\", \"mmyyyy\", \"yyyy\", \"ddmonyyyy\", \"yyyym00\"]")
    end
    
    if date_format == "ddmonyyyy"
        date = lowercase(date)
        day = date[1:2]
        month_str = date[3:5]
        year = date[6:end]
        month = month_map_UNDID[month_str]
        output = Date("$day/$month/$year", "dd/mm/yyyy")
    elseif date_format == "yyyym00"
        date = lowercase(date)
        info = split(date, "m")        
        output = Dates.format(Date("$(info[2])/$(info[1])", "mm/yyyy"), "mm/yyyy")
    elseif date_format in possible_formats        
        output = Date(date, date_format)        
    else
        error("Please specify a date_format listed here: $possible_formats_UNDID. Format 'ddmonyyyy' should look like '25dec2020' and format yyyym00 should look like '2020m12'.")
    end 

    return output 
end

function parse_date_to_string(date, date_format)
    
    # This function is basically a wrapper for Dates.format()
    # except this adds a bit more functionality so that it can 
    # take date strings in Stata formats (e.g. 25dec2020 or 2020m12) and return those as strings

    if date_format == "ddmonyyyy"
        month_dict = Dict("01" => "jan", "02" => "feb", "03" => "mar", "04" => "apr", "05" => "may", 
        "06" => "jun", "07" => "jul", "08" => "aug", "09" => "sep", "10" => "oct", 
        "11" => "nov", "12" => "dec")
        vectorized_date_object = split(string(date_object), "-")
        output = "$(vectorized_date_object[3])$(month_dict[vectorized_date_object[2]])$(vectorized_date_object[1])"        
    elseif date_format == "yyyym00"
        month_dict = Dict("01" => "m1", "02" => "m2", "03" => "m3", "04" => "m4", "05" => "m5", 
        "06" => "m6", "07" => "m7", "08" => "m8", "09" => "m9", "10" => "m10", 
        "11" => "m11", "12" => "m12")
        vectorized_date_object = split(string(date_object), "-")
        output = "$(vectorized_date_object[1])$(month_dict[vectorized_date_object[2]])"
    else
        output = Dates.format(date, date_format)
    end 

    return output
end

function parse_freq(period_str::String)
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

end 
