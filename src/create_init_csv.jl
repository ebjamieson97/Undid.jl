"""
    create_init_csv(names = String[], start_times = String[],
                    end_times = String[], treatment_times = String[];
                    covariates = false, filename = "init.csv",
                    filepath = tempdir())

The `create_init_csv()` function generates a CSV file with information
on each silo's start times, end times, and treatment times.
If parameters are not specified, the function generates a blank CSV with only
the headers.

# Details
Ensure dates are entered consistently in the same date format.
Call `undid_date_formats()` to view supported date formats. Control silos
should be marked as `"control"` in the `treatment_times` vector. If
`covariates` is `false`, no covariate column will be included in the CSV.

# Parameters
- `names::Vector{AbstractString}` 
    Vector of silo names.
- `start_times::Union{Vector{AbstractString}, AbstractString}`
    A vector of start times or a single start time.
- `end_times::Union{Vector{AbstractString}, AbstractString}`
    A vector of end times or a single end time.
- `treatment_times::Vector{AbstractString}`
    A vector of treatment times.
- `covariates::Union{Vector{AbstractString}, Bool}` 
    A vector of covariates, or, `false` (default).
- `filename::AbstractString` 
    A filename for the created initializing CSV file. Defaults to `"init.csv"`.
- `filepath::AbstractString`
    Filepath to save the CSV file. Defaults to `tempdir()`.

# Returns
A DataFrame containing the contents written to the CSV file.
The CSV file is saved in the specified directory (or in a temporary
directory by default) with the default filename `init.csv`.

"""
function create_init_csv(names = String[], start_times = String[], end_times = String[],
                         treatment_times = String[]; covariates = false, filename = "init.csv",
                         filepath = tempdir())
    
    # Check input vectors are the same length
    if !(length(names) == length(treatment_times))
        error("Input vectors for names, and treatment_times must be the same length.")       
    end 
    if length(end_times) == 1 || !(end_times isa Vector)
        if !(end_times isa Vector)
            end_times = [end_times]
        end
        end_times = end_times[1]
        end_times = fill(end_times, length(names))
    end 
    if length(start_times) == 1 || !(start_times isa Vector)
        if !(start_times isa Vector)
            start_times = [start_times]
        end
        start_times = start_times[1]
        start_times = fill(start_times, length(names))
    end 
    if (length(names) != length(start_times))
        error("If start_times is a different length than names, it must be of length 1.")
    end 
    if (length(names) != length(end_times))
        error("If start_times is a different length than names, it must be of length 1.")
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
    filepath = save_as_csv("init.csv", data, "array_of_arrays")   
    return filepath
end 