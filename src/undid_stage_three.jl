function undid_stage_three(dir_path::AbstractString;
                           agg::AbstractString = "g",
                           covariates::Bool = false,
                           interpolation::Union{Bool, AbstractString} = false,
                           nperm::Number = 1000,
                           verbose::Bool = true,
                           seed::Number = rand(1:1000000),
                           weighting::Union{Nothing, AbstractString} = nothing,
                           save_diff_data::Bool = false)

    # Set nperm to be a positive integer
    nperm = Int(round(nperm))
    if nperm < 1
        error("`nperm` must be an integer > 0")
    end 
    
    # Load in the diff dfs and convert treat column to an Int
    combined_diff_data = combine_diff_data(dir_path, save_csv = save_diff_data, interpolation = interpolation) 
    combined_diff_data.treat = convert(Vector{Int64}, combined_diff_data.treat)

    # Set diff_estimate or diff_estimate_covariates to be the Y column depending on user selection
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

    # Collect weighting option and perform checks
    badstr = Set(["missing", "na", "nan"])
    if isnothing(weighting)
        weighting = combined_diff_data.weights[1]
    end 
    if !in(weighting, ["none", "diff", "att", "both"])
            error("weighting must be one of: none, diff, att, both")
    end 
    if in(weighting, ["diff", "both"]) && any(x -> ismissing(x) || (x isa AbstractString && lowercase(x) in badstr), combined_diff_data.n)
        error("""
        Found missing values of `n` which are needed with weighting options "diff" and "both"!
        Try setting `weighting = "none"` or `weighting = "att"`
        """)
    end 
    if in(weighting, ["att", "both"]) && any(x -> ismissing(x) || (x isa AbstractString && lowercase(x) in badstr), combined_diff_data.n_t)
        error("""
        Found missing values of `n_t` which are needed with weighting options "att" and "both"!
        Try setting `weighting = "none"` or `weighting = "diff"`
        """)
    end 

    # Determine if diff_dfs represent a common or staggered adoption scenario
    if "common_treatment_time" in DataFrames.names(combined_diff_data)
        common_adoption = true
        staggered_adoption = false
    elseif "diff_times" in DataFrames.names(combined_diff_data)
        common_adoption = false
        staggered_adoption = true
    end

    # Compute ATTs and run randomization inference
    if common_adoption
        results = calculate_agg_att(combined_diff_data, agg, weighting,
                                    "common_adoption", covariates)
        results = randomization_inference_v2_undid(combined_diff_data, nperm, results,
                                                   agg, verbose, seed, weighting)
    elseif staggered_adoption
        combined_diff_data.t = [gt[2] for gt in combined_diff_data.gt]
        diff_matrix_no_RI = combined_diff_data[combined_diff_data.RI .== 0,:]
        results = calculate_agg_att(diff_matrix_no_RI, agg, weighting,
                                    "staggered_adoption")
        results =  randomization_inference_v2_undid(combined_diff_data, nperm, results,
                                                    agg, verbose, seed, weighting)
    end 

    return results 
end

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

function calculate_agg_att(combined_diff_data::DataFrame,
                           agg::AbstractString,
                           weighting::AbstractString,
                           scenario::AbstractString,
                           covariates::Bool = false)

    # Detect scenario and throw to lower level functions
    if scenario == "staggered_adoption"
        if !(in(agg, ["g", "gt", "silo", "sgt", "none"]))
            error("Please specify an aggregation scheme for ATT: \"silo\", \"gt\", \"g\", \"sgt\", or \"none\" .")
        end
        if agg == "g"
            results = calculate_agg_att_g(combined_diff_data, weighting)
        elseif agg == "gt"
            results = calculate_agg_att_gt(combined_diff_data, weighting)
        elseif agg == "silo"
            results = calculate_agg_att_silo(combined_diff_data, weighting)
        elseif agg == "sgt"
            results = calculate_agg_att_sgt(combined_diff_data, weighting)
        elseif agg == "none"
            results = calculate_agg_att_none(combined_diff_data, weighting, scenario)
        end 
        return results
    elseif scenario == "common_adoption"
        results = calculate_agg_att_none(combined_diff_data, weighting, scenario, covariates)
    end 

    return results
end   

function randomization_inference_v2_undid(diff_df::DataFrame,
                                          nperm::Int, 
                                          results::DataFrame,
                                          agg::AbstractString, 
                                          verbose::Bool, 
                                          seed::Number,
                                          weighting::AbstractString)

    # PART ONE: CREATE RANDOMIZED TREATMENT COLUMNS
    # Catch common_treatment_time
    if "common_treatment_time" in DataFrames.names(diff_df)
        diff_df.gvar = diff_df.common_treatment_time
    end
    original_treated = unique(diff_df[diff_df.treat .== 1, [:silo_name, :gvar]])
    treatment_times = original_treated.gvar
    k = nrow(original_treated)  
    treatment_states = original_treated.silo_name
    all_states = unique(diff_df.silo_name)

    n_unique_perms = compute_n_unique_assignments(treatment_times, length(all_states))
    if nperm > n_unique_perms
        @warn "'nperm' was set to $nperm but only $n_unique_perms unique permutations exist. \n 
        Setting 'nperm' to $n_unique_perms."
        nperm = n_unique_perms
    end 
    if nperm < 500
        @warn "'nperm' is less than 500!"
    end 

    randomized_diff_df = diff_df
    Random.seed!(seed)

    i = 1
    seen = Set{String}()
    pairs = zip(original_treated.silo_name, original_treated.gvar)
    key = join(sort([string(s, "-", t) for (s, t) in pairs]), "")
    push!(seen, key)
    while i < nperm
        shuffled_states = shuffle(all_states)
        new_treated_states = shuffled_states[1:k]
        new_treated_times = treatment_times[randperm(k)]
        pairs = zip(new_treated_states, new_treated_times)
        key = join(sort([string(s, "-", t) for (s, t) in pairs]), "")
        if key in seen
            continue
        end 
        assigned_tt = Dict(new_treated_states[j] => new_treated_times[j] for j in 1:k)
        new_treat = Vector{Int}(undef, nrow(diff_df))
        for row_idx in 1:nrow(diff_df)
            state = diff_df[row_idx, :silo_name]
            trt_time = diff_df[row_idx, :gvar]
            if haskey(assigned_tt, state)
                if trt_time == assigned_tt[state]
                    new_treat[row_idx] = 1
                else 
                new_treat[row_idx] = -1
                end
            else
            new_treat[row_idx] = 0
            end
        end
        randomized_diff_df[!, Symbol("treat_random_", i)] = new_treat
        i += 1
    end

    # PART TWO: COMPUTE RI_ATT & RI_ATT_SUBGROUP
    if length(unique(treatment_times)) == 1
        agg = "none"
    end 
    att_ri = Vector{Float64}(undef, nperm - 1)
    if agg == "g" 
        att_ri_cohort = Matrix{Float64}(undef, nperm - 1, length(treatment_times))
        for j in 1:nperm - 1 
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, length(treatment_times))
            for i in eachindex(treatment_times)
                # Compute sub aggregate ATT
                trt = treatment_times[i]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.gvar .== trt), :]
                att_ri_cohort[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)

                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(W))
            end 
            att_ri[j] = dot(W, att_ri_cohort[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "silo"
        att_ri_state = Matrix{Float64}(undef, nperm - 1, length(treatment_times))
        for j in 1:nperm - 1 
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, length(treatment_states))
            for i in eachindex(treatment_states)
                # Compute sub agg ATT
                trt = treatment_times[i]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.gvar .== trt), :]
                temp_treated_silos = unique(temp[temp[!, colname] .== 1, "silo_name"])
                temp_treated_silo = shuffle(temp_treated_silos)[1]
                temp = temp[(temp.silo_name .== temp_treated_silo) .|| (temp[!, colname] .== 0), :]
                att_ri_state[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)

                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 

            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(W))
            end 
            att_ri[j] = dot(W, att_ri_state[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end 
        end 
    elseif agg == "gt"
        unique_diffs = unique(select(diff_df[diff_df.treat .== 1,:], :gvar, :t))
        att_ri_simple = Matrix{Float64}(undef, nperm - 1, nrow(unique_diffs))
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, nrow(unique_diffs))
            for i in 1:nrow(unique_diffs)
                # Compute sub agg ATT
                t = unique_diffs[i,"t"]
                gvar = unique_diffs[i,"gvar"]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.t .== t) .&& (diff_df.gvar .== gvar), :]
                att_ri_simple[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)

                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(W))
            end 
            att_ri[j] = dot(W, att_ri_simple[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "sgt"
        unique_diffs = unique(select(diff_df[diff_df.treat .== 1,:], :silo_name, :t, :gvar))
        att_ri_sgt = Matrix{Float64}(undef, nperm - 1, nrow(unique_diffs))
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            W = Vector{Float64}(undef, nrow(unique_diffs))
            for i in 1:nrow(unique_diffs)
                # Compute sub agg ATT
                t = unique_diffs[i,"t"]
                gvar = unique_diffs[i,"gvar"]
                temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.t .== t) .&& (diff_df.gvar .== gvar), :]
                temp_treated_silos = unique(temp[temp[!, colname] .== 1, "silo_name"])
                temp_treated_silo = shuffle(temp_treated_silos)[1]
                temp = temp[(temp.silo_name .== temp_treated_silo) .|| (temp[!, colname] .== 0), :]
                att_ri_sgt[j,i] = compute_ri_sub_agg_att(temp, weighting, colname)

                # Compute weights
                if in(weighting, ["att", "both"])
                    W[i] = sum(temp[temp[!, colname] .== 1, "n_t"])
                end 
            end

            # Compute aggregate ATT
            if in(weighting, ["att", "both"])
                W ./= sum(W)
            elseif in(weighting, ["none", "diff"])
                W .= (1 / length(W))
            end 
            att_ri[j] = dot(W, att_ri_sgt[j,:])
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    elseif agg == "none"
        for j in 1:nperm - 1
            colname = Symbol("treat_random_$(j)")
            temp = diff_df[diff_df[!, colname] .!= -1,:]
            att_ri[j] = compute_ri_sub_agg_att(temp, weighting, colname)
            if verbose && j % 100 == 0
                println("Completed $(j) of $(nperm - 1) permutations")
            end
        end 
    end

    # PART THREE: COMPUTE P-VALS BASED ON RI_ATT & RI_ATT_SUBGROUP
    agg_att = results.agg_att[1]
    if agg == "g"
        results.ri_pval_att_g = Vector{Union{Missing, Float64}}(missing, nrow(results))
        for i in eachindex(treatment_times)
            sub_agg_att = results[results.gvar .== treatment_times[i], "att_g"][1]
            results[results.gvar .== treatment_times[i], "ri_pval_att_g"] .= (sum(abs.(att_ri_cohort[:,i]) .> abs(sub_agg_att))) / length(att_ri_cohort[:,i])
        end
    elseif agg == "silo"
        results.ri_pval_att_s = Vector{Union{Missing, Float64}}(missing, nrow(results))
        for i in eachindex(treatment_states)
            sub_agg_att = results[results.silos .== treatment_states[i], "att_s"][1]
            results[results.silos .== treatment_states[i], "ri_pval_att_s"] .= (sum(abs.(att_ri_state[:,i]) .> abs(sub_agg_att))) / length(att_ri_state[:,i])
        end 
    elseif agg == "gt"
        results.ri_pval_att_gt = Vector{Union{Missing, Float64}}(missing, nrow(results))
        for i in 1:nrow(unique_diffs)
            t = unique_diffs[i,"t"]
            gvar = unique_diffs[i,"gvar"]
            sub_agg_att = results[(results.t .== t) .&& (results.g .== gvar), "att_gt"][1]
            results[(results.t .== t) .&& (results.g .== gvar), "ri_pval_att_gt"] .= (sum(abs.(att_ri_simple[:,i]) .> abs(sub_agg_att))) / length(att_ri_simple[:,i])
        end
    elseif agg == "sgt"
        for i in 1:nrow(unique_diffs)
            t = unique_diffs[i,"t"]
            gvar = unique_diffs[i, "gvar"]
            state = unique_diffs[i, "silo_name"]
            sub_agg_att = results[(results.t .== t) .&& (results.g .== gvar) .&& (results.silo_name .== state), "att_sgt"][1]
            results[(results.t .== t) .&& (results.g .== gvar) .&& (results.silo_name .== state), "ri_pval_att_sgt"] .= (sum(abs.(att_ri_sgt[:,i]) .> abs(sub_agg_att))) / length(att_ri_sgt[:,i])
        end 
    end
    
    results.ri_pval_agg_att[1] = ((sum(abs.(att_ri) .> abs(agg_att))) / length(att_ri))
    results.nperm[1] = nperm - 1
    return results
end

function compute_n_unique_assignments(treatment_times::Vector, total_n_states::Number)
    # This computes the combinations formula * the multinomial coefficient
    n_assignments = length(treatment_times)
    unique_assignments = unique(treatment_times)
    num = factorial(big(total_n_states))           
    den = factorial(big(total_n_states - n_assignments))
    for m in unique_assignments
        n_m = sum(treatment_times .== m)
        den *= factorial(big(n_m))
    end
    return num ÷ den                  
end

function calculate_agg_att_g(diff_df::DataFrame,
                             weighting::AbstractString)

    # Grab all of the gvars and create preallocation dataframe for the results
    gvars = sort!(unique(diff_df.gvar))
    n_rows = length(gvars)
    results = DataFrame(gvar = Vector{Date}(undef, n_rows),
                        att_g = Vector{Float64}(undef, n_rows),
                        agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jackknife_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        nperm = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_g_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_g_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_g_se_jackknife = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_g_jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_att_g = Vector{Union{Missing, Float64}}(missing, n_rows),
                        weights = Vector{Union{Nothing, Float64}}(nothing, n_rows))
    
    # Compute all of the sub aggregate level ATTs
    for i in eachindex(gvars)
        trt = gvars[i]
        temp = diff_df[diff_df.gvar .== trt, :]
        X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
        Y = convert(Vector{Float64}, temp.y)
        if in(weighting, ["diff", "both"])
            W = convert(Vector{Float64}, temp.n)
            W ./= sum(W)
        else
            W = fill(nothing, nrow(temp))
        end            
        results.gvar[i] = trt
        result_dict = regression_results(X, Y, W = W)
        results.att_g[i] = result_dict["beta_hat"]
        results.att_g_se[i] = result_dict["beta_hat_se"]
        results.att_g_pval[i] = result_dict["pval_att"] 
        results.att_g_se_jackknife[i] = result_dict["beta_hat_se_jknife"]
        results.att_g_jknife_pval[i] = result_dict["pval_att_jknife"]
        if weighting in ["att", "both"]
            results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
        end
    end

    # Scale weights and convert to either Float64 or Nothing
    results = scale_weights_final(results, weighting)

    # Compute the aggregate level ATT
    Y = convert(Vector{Float64}, results.att_g)
    results = compute_aggregate_att_from_sub_aggregates(results, Y)

    return results
end 

function calculate_agg_att_gt(diff_df::DataFrame, weighting::AbstractString)
    
    # Grab gt combinations and create preallocation df
    unique_diffs = unique(select(diff_df, :t, :gvar))
    date_format = diff_df.date_format[1]
    n_rows = nrow(unique_diffs)

    results = DataFrame(gt = Vector{String}(undef, n_rows),
                        t = Vector{Date}(undef, n_rows),
                        g = Vector{Date}(undef, n_rows),
                        att_gt = Vector{Float64}(undef, n_rows),
                        agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jackknife_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        nperm = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_gt_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_gt_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_gt_se_jackknife = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_gt_jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_att_gt = Vector{Union{Missing, Float64}}(missing, n_rows),
                        weights = Vector{Union{Nothing, Float64}}(nothing, n_rows))

    # Compute all of the sub aggregate level ATTs
    for i in 1:n_rows
        t = unique_diffs[i,"t"]
        gvar = unique_diffs[i, "gvar"]
        temp = diff_df[(diff_df.t .== t) .& (diff_df.gvar .== gvar), :]
        X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
        Y = convert(Vector{Float64}, temp.y)
        if in(weighting, ["diff", "both"])
            W = convert(Vector{Float64}, temp.n)
            W ./= sum(W)
        else
            W = fill(nothing, nrow(temp))
        end 
        results.t[i] = t
        results.gt[i] = string(parse_date_to_string(gvar, date_format), ";", parse_date_to_string(t, date_format))
        results.g[i] = gvar
        result_dict = regression_results(X, Y, W = W)
        results.att_gt[i] = result_dict["beta_hat"]
        results.att_gt_se[i] = result_dict["beta_hat_se"]
        results.att_gt_pval[i] = result_dict["pval_att"] 
        results.att_gt_se_jackknife[i] = result_dict["beta_hat_se_jknife"]
        results.att_gt_jknife_pval[i] = result_dict["pval_att_jknife"]
        if weighting in ["att", "both"]
            results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
        end
    end
    
    # Scale weights and convert to either Float64 or Nothing
    results = scale_weights_final(results, weighting)

    # Sort results df by gt
    sort!(results, [:t])
    sort!(results, [:g])

    # Compute the aggregate level ATT
    Y = convert(Vector{Float64}, results.att_gt)
    results = compute_aggregate_att_from_sub_aggregates(results, Y)
            
end 

function calculate_agg_att_silo(diff_df::DataFrame, weighting::AbstractString)

    # Grab treated states are create preallocation dataframe for results
    treated_states = unique(select(diff_df[diff_df.treat .== 1,:], :silo_name, :gvar))
    n_rows = nrow(treated_states)
    results = DataFrame(silos = Vector{String}(undef, n_rows),
                        att_s = Vector{Float64}(undef, n_rows),
                        agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jackknife_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        nperm = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_s_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_s_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_s_se_jackknife = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_s_jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_att_s = Vector{Union{Missing, Float64}}(missing, n_rows),
                        weights = Vector{Union{Nothing, Float64}}(nothing, n_rows))

    # Compute all of the sub aggregate level ATTs
    for i in 1:n_rows
        silo = treated_states[i,"silo_name"]
        gvar = treated_states[i, "gvar"]
        temp = diff_df[diff_df.silo_name .== silo.|| (diff_df.treat .== 0 .&& diff_df.gvar .== gvar),:]
        X = convert(Matrix{Float64},(hcat(ones(nrow(temp)), temp.treat)))
        Y = convert(Vector{Float64}, temp.y)
        if in(weighting, ["diff", "both"])
            W = convert(Vector{Float64}, temp.n)
            W ./= sum(W)
        else
            W = fill(nothing, nrow(temp))
        end 
        results.silos[i] = silo
        result_dict = regression_results(X, Y, W = W)
        results.att_s[i] = result_dict["beta_hat"]
        results.att_s_se[i] = result_dict["beta_hat_se"]
        results.att_s_pval[i] = result_dict["pval_att"] 
        results.att_s_se_jackknife[i] = result_dict["beta_hat_se_jknife"]
        results.att_s_jknife_pval[i] = result_dict["pval_att_jknife"]
        if weighting in ["att", "both"]
            results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
        end
    end

    # Scale weights and convert to either Float64 or Nothing
    results = scale_weights_final(results, weighting)

    # Sort by silo_name 
    results.tuple_state = custom_sort_order.(results.silos)
    sort!(results,[order(:tuple_state)])
    select!(results, Not([:tuple_state]))

    # Compute the aggregate ATT based on the sub-aggregate computations
    Y = convert(Vector{Float64}, results.att_s)
    results = compute_aggregate_att_from_sub_aggregates(results, Y)

    return results
end 

function calculate_agg_att_sgt(diff_df::DataFrame, weighting::AbstractString)

    # Grab sgt combinations and create preallocation df
    unique_sgt = unique(select(diff_df[diff_df.treat .== 1,:], :silo_name, :t, :gvar))
    date_format = diff_df.date_format[1]
    n_rows = nrow(unique_sgt)
    results = DataFrame(silo_name = Vector{String}(undef, n_rows),
                        g = Vector{Date}(undef, n_rows),
                        t = Vector{Date}(undef, n_rows),
                        att_sgt = Vector{Float64}(undef, n_rows),
                        agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        agg_att_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jackknife_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, n_rows),
                        nperm = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_sgt_se = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_sgt_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_sgt_se_jackknife = Vector{Union{Missing, Float64}}(missing, n_rows),
                        att_sgt_jknife_pval = Vector{Union{Missing, Float64}}(missing, n_rows),
                        ri_pval_att_sgt = Vector{Union{Missing, Float64}}(missing, n_rows),
                        weights = Vector{Union{Nothing, Float64}}(nothing, n_rows),
                        sgt = Vector{String}(undef, n_rows))

    # Calculate sub aggregate level ATTs
    for i in 1:n_rows
        silo = unique_sgt[i, "silo_name"]
        gvar = unique_sgt[i, "gvar"]
        t = unique_sgt[i, "t"]
        temp = diff_df[((diff_df.silo_name .== silo) .&& (diff_df.gvar .== gvar) .&& (diff_df.t .== t)) .|| 
                       ((diff_df.treat .== 0) .&& (diff_df.gvar .== gvar) .&& (diff_df.t .== t)),:]
        X = convert(Matrix{Float64}, hcat(ones(nrow(temp)), temp.treat))
        Y = convert(Vector{Float64}, temp.y)
        if in(weighting, ["diff", "both"])
            W = convert(Vector{Float64}, temp.n)
            W ./= sum(W)
        else
            W = fill(nothing, nrow(temp))
        end 
        results.silo_name[i] = silo
        results.g[i] = gvar
        results.t[i] = t
        results.sgt[i] = string(silo, ";", parse_date_to_string(gvar, date_format), ";", parse_date_to_string(t, date_format))
        result_dict = regression_results(X, Y, W = W)
        results.att_sgt[i] = result_dict["beta_hat"]
        results.att_sgt_se[i] = result_dict["beta_hat_se"]
        results.att_sgt_pval[i] = result_dict["pval_att"] 
        results.att_sgt_se_jackknife[i] = result_dict["beta_hat_se_jknife"]
        results.att_sgt_jknife_pval[i] = result_dict["pval_att_jknife"]
        if weighting in ["att", "both"]
            results.weights[i] = sum(temp[temp.treat .== 1, "n_t"])
        end
    end

    # Scale weights and convert to either Float64 or Nothing
    results = scale_weights_final(results, weighting)

    # Sort dataframe of results
    sort!(results, [:t])
    sort!(results, [:g])

    # Compute the aggregate ATT based on the sub-aggregate computations
    Y = convert(Vector{Float64}, results.att_sgt)
    results = compute_aggregate_att_from_sub_aggregates(results, Y)
    
    return results
end 

function calculate_agg_att_none(diff_df::DataFrame, weighting::AbstractString,
                                scenario::AbstractString = "staggered_adoption",
                                covariates::Bool = false)
    
    # Run final regrsesion -- each diff weighted equally
    results = DataFrame(agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                        agg_att_se = Vector{Union{Missing, Float64}}(missing, 1),
                        agg_att_pval = Vector{Union{Missing, Float64}}(missing, 1),
                        jackknife_se = Vector{Union{Missing, Float64}}(missing, 1),
                        jknife_pval = Vector{Union{Missing, Float64}}(missing, 1),
                        ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, 1),
                        nperm = Vector{Union{Missing, Float64}}(missing, 1))
    
    # Get X matrix(intercept + treat) and Y vector (diffs)
    X = convert(Matrix{Float64}, hcat(ones(nrow(diff_df)), diff_df.treat))
    Y = convert(Vector{Float64}, diff_df.y)

    # Assign weighting
    if in(weighting, ["diff", "both"])
        W = convert(Vector{Float64}, diff_df.n)
        W ./= sum(W)
    else
        W = fill(nothing, nrow(diff_df))
    end 

    # Compute aggregate ATT directly from diffs 
    result_dict = regression_results(X, Y, W = W)
    results.agg_att[1] = result_dict["beta_hat"]

    # Deal with corner case for staggered adoption
    if scenario == "common_adoption" && nrow(diff_df) == 2
        if covariates == false 
            se = sqrt(sum(convert(Vector{Float64}, diff_df.diff_var)))
        elseif covariates == true 
            se = sqrt(sum(convert(Vector{Float64}, diff_df.diff_var_covariates)))
        end
        results.agg_att_se[1] = se
    else
        results.agg_att_se[1] = result_dict["beta_hat_se"]
    end

    results.agg_att_pval[1] = result_dict["pval_att"]
    results.jackknife_se[1] = result_dict["beta_hat_se_jknife"]
    results.jknife_pval[1] = result_dict["pval_att_jknife"]

    return results

end 

function regression_results(X::Matrix{<:Number}, Y::Vector{<:Number};
                                  W::Vector{T} where T <: Union{Nothing, Number} = [nothing])
    beta_hat = nothing
    beta_hat_cov = nothing
    ncolx = size(X, 2)

    # Run OLS (normally if weights aren't provided, and scale (X,Y) -> (Xw,Yw) otherwise)
    if eltype(W) <: Nothing
        try
            beta_hat = (X \ Y) 
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat = pinv(X' * X) * X' * Y
        end
        resid = Y - X * beta_hat
        omega = Diagonal(resid .^ 2)
        try
            beta_hat_cov = inv(X' * X) * (X' * omega * X) * inv(X' * X)
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat_cov = pinv(X' * X) * (X' * omega * X) * pinv(X' * X)
        end 
        beta_hat_se_jknife = compute_jknife_se(X, Y, beta_hat[ncolx]) 
    elseif eltype(W) <: Number
        sw = sqrt.(W)           
        Xw = X .* sw            
        Yw = Y .* sw
        try
            beta_hat = (Xw) \ (Yw) 
        catch e
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat = pinv(Xw' * Xw) * Xw' * Yw
        end
        resid_w = Yw - Xw * beta_hat
        Ωw = Diagonal(resid_w.^2)
        try
            beta_hat_cov = inv(Xw'Xw) * (Xw' * Ωw * Xw) * inv(Xw'Xw)
        catch e 
            @warn "Direct solve failed, using pseudoinverse: $e"
            beta_hat_cov = pinv(Xw'Xw) * (Xw' * Ωw * Xw) * pinv(Xw'Xw)
        end 
        beta_hat_se_jknife = compute_jknife_se(Xw, Yw, beta_hat[ncolx]) 
    end 
    beta_hat_var = diag(beta_hat_cov)
    beta_hat_se = sqrt(beta_hat_var[ncolx]) 
    dof = length(Y) - ncolx
    pval_att = dof > 0 ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se))) : missing 
    pval_att_jknife = dof > 0 && !ismissing(beta_hat_se_jknife) ? 2 * (1 - cdf(TDist(dof), abs(beta_hat[ncolx] / beta_hat_se_jknife))) : missing
    result_dict = Dict("beta_hat" => beta_hat[ncolx], "beta_hat_se" => beta_hat_se, "pval_att" => pval_att,
                       "beta_hat_se_jknife" => beta_hat_se_jknife, "pval_att_jknife" => pval_att_jknife)
    return result_dict
end 

function compute_jknife_se(X::Matrix{<:Number}, Y::Vector{<:Number}, original_att::Number)
    
    n = length(Y)
    if n == 1 
        return missing
    end 
    jknife_beta = Vector{Float64}(undef, n)
    ncolx = size(X,2)
    treat_count = sum(X[:,ncolx] .!= 0)
    control_count = sum(X[:,ncolx] .== 0)
    if (treat_count < 2 || control_count < 2) & (ncolx > 1)
        return missing
    end
    for i in eachindex(Y)
        idx = [1:i-1; i+1:size(X, 1)]
        X_sub = X[idx, :]
        Y_sub = Y[idx]
        jknife_beta[i] = (X_sub \ Y_sub)[ncolx]
    end 
    jknife_se = sqrt(sum((jknife_beta .- original_att).^2) * ((n - 1) / n))
    return jknife_se
end 

function scale_weights_final(results::DataFrame, weighting::AbstractString)
    
    # Scale weights and convert to either Float64 or Nothing
    if weighting in ["att", "both"]
        results.weights ./= sum(results.weights)
        results.weights = convert(Vector{Float64}, results.weights)
    elseif weighting in ["diff", "none"]
        results.weights = convert(Vector{Nothing}, results.weights)
    end
    return results
end 

function custom_sort_order(s)
    parsed = tryparse(Int, s)
    if isnothing(parsed)
        (0, lowercase(s))
    else
        (1, parsed)
    end 
end 

function compute_aggregate_att_from_sub_aggregates(results::DataFrame, Y::Vector{<:Number})
    
    # Compute the aggregate level ATT
    X = ones(nrow(results), 1)
    W = results.weights
    result_dict = regression_results(X, Y, W = W)
    results.agg_att[1] = result_dict["beta_hat"]
    results.agg_att_se[1] = result_dict["beta_hat_se"]
    results.agg_att_pval[1] = result_dict["pval_att"]
    results.jackknife_se[1] = result_dict["beta_hat_se_jknife"]
    results.jknife_pval[1] = result_dict["pval_att_jknife"]

    return results
end 

function compute_ri_sub_agg_att(temp::DataFrame, weighting::AbstractString, colname::Symbol)

    X = convert(Vector{Float64}, temp[!, colname])
    Y = convert(Vector{Float64}, temp.y)
    if in(weighting, ["both", "diff"])
        W_diff  = convert(Vector{Float64}, temp.n)
        W_diff ./= sum(W_diff)
        sub_agg_att = (dot(W_diff[X .== 1], Y[X .== 1]) / sum(W_diff[X .== 1])) - 
                             (dot(W_diff[X .== 0], Y[X .== 0]) / sum(W_diff[X .== 0]))
    elseif in(weighting, ["att", "none"])
        sub_agg_att = mean(Y[X .== 1]) - mean(Y[X .== 0])
    end 
    return sub_agg_att
end 