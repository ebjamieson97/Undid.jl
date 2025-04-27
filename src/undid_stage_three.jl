function undid_stage_three(dir_path::AbstractString; agg::AbstractString = "silo", covariates::Bool = false, save_all_csvs::Bool = false, interpolation = false, weights::Bool = true, nperm = 1000, verbose::Bool = true, seed::Number = rand(1:1000000))

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
        combined_diff_data.t = [gt[2] for gt in combined_diff_data.gt]
        diff_matrix_no_RI = combined_diff_data[combined_diff_data.RI .== 0,:]
        
        # Calculate att by specified aggregation method 
        results = calculate_agg_att_df(diff_matrix_no_RI; agg = agg, covariates = covariates, save_all_csvs = save_all_csvs, printinfo = true, weights = weights, nperm = nperm)  
        results = randomization_inference_v2_undid(combined_diff_data, nperm, results, agg, verbose, seed, covariates)
        return results
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
            results.tuple_silos = custom_sort_order.(results.silos)
            sort!(results,[order(:tuple_silos)])
            select!(results, Not([:tuple_silos]))         
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
                g_vec = [gt[1] for gt in gt_vec]
                t_vec = [gt[2] for gt in gt_vec]
                results = DataFrame(g = g_vec, t = t_vec, gt = gt_vec, gt_n = Float64.(ATT_n), jack_n = Float64.(jack_n), att_gt = Float64.(ATT_vec), att_gt_se = Float64.(ATT_vec_se), att_gt_se_jackknife = Float64.(ATT_vec_se_jack), agg_att = vcat([ATT_gt], fill(NaN, length(ATT_vec) - 1)), agg_ATT_se = vcat([ATT_gt_se], fill(NaN, length(ATT_vec) - 1)), jackknife_se = vcat([jackknife_SE], fill(NaN, length(ATT_vec) - 1)))
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
                results = DataFrame(gvar = gvars, g_n = Float64.(ATT_g_n), jack_n = Float64.(jack_n), att_g = Float64.(ATT_g), att_g_se = Float64.(ATT_g_se), att_g_se_jackknife = Float64.(ATT_g_se_jack),
                                    agg_att = vcat([agg_ATT], fill(NaN, length(ATT_g) - 1)), agg_ATT_se = vcat([agg_ATT_se], fill(NaN, length(ATT_g) - 1)), jackknife_se = vcat([jackknife_SE], fill(NaN, length(ATT_g) - 1)))
            end
        end 
        return results
    end   

end   

function randomization_inference_v2_undid(diff_df::DataFrame, nperm::Int, results::DataFrame,
    agg::AbstractString, verbose::Bool, seed::Number, covariates::Bool)

# PART ONE: CREATE RANDOMIZED TREATMENT COLUMNS
original_treated = unique(diff_df[diff_df.treat .== 1, [:silo_name, :gvar]])
k = nrow(original_treated)  
treatment_times = original_treated.gvar
treatment_states = original_treated.silo_name
all_states = unique(diff_df.silo_name)

# Select whether to use diff_estimate or diff_estimate_covariates
if covariates == true
    diff_df.diff = diff_df.diff_estimate_covariates
elseif covariates == false
    diff_df.diff = diff_df.diff_estimate
end 

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
if agg in ["g", "gt"]
agg = "unweighted"
end
end 
att_ri = Vector{Float64}(undef, nperm - 1)
if agg == "g" 
att_ri_cohort = Matrix{Float64}(undef, nperm - 1, length(treatment_times))
for j in 1:nperm - 1 
colname = Symbol("treat_random_$(j)")
for i in eachindex(treatment_times)
trt = treatment_times[i]
temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.gvar .== trt), :]
X = convert(Vector{Float64}, temp[!, colname])
Y = convert(Vector{Float64}, temp.diff)
att_ri_cohort[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])
end
att_ri[j] = mean(att_ri_cohort[j,:])
if verbose && j % 100 == 0
println("Completed $(j) of $(nperm - 1) permutations")
end
end 
elseif agg == "silo"
att_ri_state = Matrix{Float64}(undef, nperm - 1, length(treatment_times))
for j in 1:nperm - 1 
colname = Symbol("treat_random_$(j)")
for i in eachindex(treatment_states)
trt = treatment_times[i]
temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.gvar .== trt), :]
temp_treated_silos = unique(temp[temp[!, colname] .== 1, "silo_name"])
temp_treated_silo = shuffle(temp_treated_silos)[1]
temp = temp[(temp.silo_name .== temp_treated_silo) .|| (temp[!, colname] .== 0), :]
X = convert(Vector{Float64}, temp[!, colname])
Y = convert(Vector{Float64}, temp.diff)
att_ri_state[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])
end
att_ri[j] = mean(att_ri_state[j,:])
if verbose && j % 100 == 0
println("Completed $(j) of $(nperm - 1) permutations")
end 
end 
elseif agg == "gt"
unique_diffs = unique(select(diff_df[diff_df.treat .== 1,:], :gvar, :t))
att_ri_simple = Matrix{Float64}(undef, nperm - 1, nrow(unique_diffs))
for j in 1:nperm - 1
colname = Symbol("treat_random_$(j)")
for i in 1:nrow(unique_diffs)
t = unique_diffs[i,"t"]
gvar = unique_diffs[i,"gvar"]
temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.t .== t) .&& (diff_df.gvar .== gvar), :]
X = convert(Vector{Float64}, temp[!, colname])
Y = convert(Vector{Float64}, temp.diff)
att_ri_simple[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])
end
att_ri[j] = mean(att_ri_simple[j,:])
if verbose && j % 100 == 0
println("Completed $(j) of $(nperm - 1) permutations")
end
end 
elseif agg == "sgt"
unique_diffs = unique(select(diff_df[diff_df.treat .== 1,:], :state, :t, :r1, :treated_time))
att_ri_sgt = Matrix{Float64}(undef, nperm - 1, nrow(unique_diffs))
for j in 1:nperm - 1
colname = Symbol("treat_random_$(j)")
for i in 1:nrow(unique_diffs)
t = unique_diffs[i,"t"]
r1 = unique_diffs[i,"r1"]
temp = diff_df[(diff_df[!, colname] .!= -1) .&& (diff_df.t .== t) .&& (diff_df.r1 .== r1), :]
temp_treated_silos = unique(temp[temp[!, colname] .== 1, "state"])
temp_treated_silo = shuffle(temp_treated_silos)[1]
temp = temp[(temp.state .== temp_treated_silo) .|| (temp[!, colname] .== 0), :]
X = convert(Vector{Float64}, temp[!, colname])
Y = convert(Vector{Float64}, temp.diff)
att_ri_sgt[j,i] = mean(Y[X .== 1]) - mean(Y[X .== 0])
end
att_ri[j] = mean(att_ri_sgt[j,:])
if verbose && j % 100 == 0
println("Completed $(j) of $(nperm - 1) permutations")
end
end 
elseif agg == "unweighted"
for j in 1:nperm - 1
colname = Symbol("treat_random_$(j)")
temp = diff_df[diff_df[!, colname] .!= -1,:]
X = convert(Vector{Float64}, temp[!, colname])
Y = convert(Vector{Float64}, temp.diff)
att_ri[j] = mean(Y[X .== 1]) - mean(Y[X .== 0])
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
gvar = unique_diffs[i, "treated_time"]
state = unique_diffs[i, "state"]
sub_agg_att = results[(results.t .== t) .&& (results.gvar .== gvar) .&& (results.state .== state), "att_sgt"][1]
results[(results.t .== t) .&& (results.gvar .== gvar) .&& (results.state .== state), "ri_pval_att_sgt"] .= (sum(abs.(att_ri_sgt[:,i]) .> abs(sub_agg_att))) / length(att_ri_sgt[:,i])
end 
end
results.ri_pval_agg_att = Vector{Union{Missing, Float64}}(missing, nrow(results))
results.ri_pval_agg_att[1] = ((sum(abs.(att_ri) .> abs(agg_att))) / length(att_ri))
results.nperm = Vector{Union{Missing, Float64}}(missing, nrow(results))
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


function custom_sort_order(s)
    parsed = tryparse(Int, s)
    if isnothing(parsed)
        (0, lowercase(s))
    else
        (1, parsed)
    end 
end 