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
    
    # This function is basically a wrapper Dates.format()
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
