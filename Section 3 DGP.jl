# First generating synthetic data to work with
# Using a similar DGP as presented in section 3.2 of the UNDID paper

using Statistics # This is just to make sure that functions like mean() work
using Random # For setting seeds

rows = 50
columns = 10

ATTs = Vector{Any}(undef, 1000)
Random.seed!(123)
# Computing ATT with 1000 different random samples with the same treatment effect
for i in 1:length(ATTs)

    case1 = randn(rows, columns)

    case2 = copy(case1)
    case2[end-16:end, end-4:end] .+= 0.1

    # If we assume poolability of data, we can compute ATT as normal
    ATT = (mean(case2[end-16:end, end-4:end]) - mean(case2[end-16:end, 1:5])) - (mean(case2[1:end-17, end-4: end]) - mean(case2[1:end-17, 1:5]))
    ATTs[i] = ATT
end

mean(ATTs)
using Plots
# Plotting results of ATT
histogram(ATTs, nbins = 100)

# We can also rearrange the matrix and compute difference in differences via OLS
# I think this will be most straightforward way to compute UNDID going forward
ATTs_2 = Vector{Any}(undef, 1000)
ATTs_var_2 = Vector{Any}(undef, 1000)
Random.seed!(125)
for i in 1:length(ATTs_2)

    case1 = randn(rows, columns)

    case2 = copy(case1)
    case2[end-16:end, end-4:end] .+= 0.1
    Y = vec(copy(case2))

    # Create matrices indicating before/after treatment period and treatment group status
    time_indicator = ones(Int, size(case2))  
    time_indicator[:, 1:5] .= 0 

    treated_group_indicator = zeros(Int, size(case2))
    treated_group_indicator[end-16:end, :] .= 1
    
    X = hcat(fill(1, 500), vec(time_indicator), vec(treated_group_indicator), vec(time_indicator .* treated_group_indicator))
    beta_hat = X\Y
    ATTs_2[i] = beta_hat[4]

    # Also want to calculate SE(DID) ## 
    # Calculate the residuals
    residuals = Y - X * beta_hat
    # Estimate the variance of the residuals
    sigma_sq = sum(residuals.^2) / (length(residuals) - length(beta_hat))
    # Compute the covariance matrix of the coefficients
    cov_beta_hat = sigma_sq * inv(X' * X)
    
    
    # Could write a little for loop at this point to extract variance of each coefficient
    # But for now just grabbing SE of beta_4, i.e. coefficient on DID 
    ATTs_var_2[i] = sqrt(cov_beta_hat[4,4])
end 

mean(ATTs_2)
mean(ATTs_var_2)

ATTs_UNDID = Vector{Any}(undef, 1000)
ATTs_UNDID_var = Vector{Any}(undef, 1000)
Random.seed!(125)
for  i in 1:length(ATTs_UNDID)
    ## Now compute ATT via UNDID method ##
    # Generate data
    case1 = randn(rows, columns)
    case2 = copy(case1)
    case2[end-16:end, end-4:end] .+= 0.1

    # Seperate into 'silos' of treated and control
    treated = copy(case2[end-16:end, :])
    control = copy(case2[1:end-17,:])

    # Run regression locally on treated silo
    treated_post_indicator = [zeros(17, 5) ones(17, 5)]
    treated_pre_indicator = [ones(17, 5) zeros(17, 5)]
    X = hcat(fill(1,170), vec(treated_post_indicator))
    Y = vec(copy(treated))
    beta_hat = X\Y
    treated_gamma = beta_hat[2]

    # Also want to calculate SE(gamma) ## 
    # Calculate the residuals
    residuals = Y - X * beta_hat
    # Estimate the variance of the residuals
    sigma_sq = sum(residuals.^2) / (length(residuals) - length(beta_hat))
    # Compute the covariance matrix of the coefficients
    cov_beta_hat = sigma_sq * inv(X' * X)

    treated_gamma_var = cov_beta_hat[2,2]

    # Run regression locally on control silo
    control_post_indicator = [zeros(33, 5) ones(33, 5)]
    control_pre_indicator = [ones(33, 5) zeros(33, 5)]
    X = hcat(fill(1,330), vec(control_post_indicator))
    Y = vec(copy(control))
    beta_hat = X\Y 
    control_gamma = beta_hat[2] 

    # Also want to calculate SE(gamma) ## 
    # Calculate the residuals
    residuals = Y - X * beta_hat
    # Estimate the variance of the residuals
    sigma_sq = sum(residuals.^2) / (length(residuals) - length(beta_hat))
    # Compute the covariance matrix of the coefficients
    cov_beta_hat = sigma_sq * inv(X' * X)

    control_gamma_var = cov_beta_hat[2,2]

    ATTs_UNDID[i] = treated_gamma - control_gamma

    ATTs_UNDID_var[i] = sqrt(treated_gamma_var + control_gamma_var)
end

# Here we can see the ATT computed via UNDID yields
# the same ATT using DID
mean(ATTs_UNDID)
mean(ATTs_UNDID_var)
