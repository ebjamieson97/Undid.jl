"""
    create_diff_df(csv::AbstractString; covariates = false, date_format = false,
                   freq = false, freq_multiplier = false, weights = "standard")

Creates a `empty_diff_df.csv` which lists all of the differences that
need to calculated at each silo in order to compute the aggregate ATT.
The `empty_diff_df.csv` is then to be sent out to each silo to be filled out.

# Details
Ensure that dates in the initializing CSV file, by default named `init.csv`,
are entered consistently in the same date format. Call `undid_date_formats()`
to see a list of valid date formats.

Covariates specified when calling `create_diff_df()` will override any covariates
specified in the `init.csv`.

# Parameters
- `csv::AbstractString` 
    Filepath to the `init.csv`.
- `date_format::AbstractString`
    A string specifying the date format used in the `init.csv`. Call
    `undid_date_formats()` to see a list of valid date formats.
- `freq::AbstractString`
    A string indicating the length of the time periods to be used
    when computing the differences in mean outcomes between periods at each
    silo. Options are: `"yearly"`, `"monthly"`, `"weekly"`, or `"daily"`.
- `covariates::Union{Vector{AbstractString}, AbstractString, Bool}` 
    A vector of strings or a single string specifying covariate(s) to be
    considered at each silo. If `false` (default) uses covariates from the `init.csv`.
- `freq_multiplier::Union{Real, Bool}`
    A Real value or `false` (default) which specifies if the frequency
    should be multiplied by a non-zero integer.
- `weights::AbstractString`
    A string indicating the weighting to use in the case of
    common adoption. The `"standard"` (default) weight is calculated as:

    ``w_s = \\frac{N_s^{\\text{post}}}{N_s^{\\text{post}} + N_s^{\\text{pre}}}``
    
    Options are: `"standard"`.
- `filename::AbstractString`
    A string filename for the created CSV file. Defaults to `"empty_diff_df.csv"`.
- `filepath::AbstractString`
    Filepath to save the CSV file. Defaults to `tempdir()`.

# Returns
A DataFrame detailing the silo and time combinations for which
differences must be calculated in order to compute the aggregate ATT. A
CSV copy is saved to the specified directory which is then to be sent out
to each silo.
"""
function create_diff_df(csv::AbstractString, date_format::AbstractString, freq::AbstractString;
                        covariates = false,   freq_multiplier = false, weights = "standard",
                        filename = "empty_diff_df.csv", filepath = tempdir())

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
