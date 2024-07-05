# This script is a test for implementing the Callaway and Sant'Anna procedure (CSDID)
# And UNDID as reported in section 6.2 of the UNDID paper

# Basically I just want to make sure I can replicate the reported ATTs in the UNDID paper section 6.2
# Once I have this done, I'll have the math-y portion of the coding done and I'll just need
# to work on optimizing code efficiency and handling file types, timeseries and repeated cross sectional data types etc.

using StatFiles # For reading .dta
using DataFrames # For dataframes
using Statistics # Part of base Julia


df = DataFrame(load("C:/Users/Eric Bruce Jamieson/Documents/Dalhousie Work/Merit Data Sent from Matt/merit_data_changed.dta"))

# Note that the 'merit' column indicates 1 once receiving treatment and 0 for pre-treatment periods
# First I'd like to try results shown in Panel A (Table 2, Page 36)... 
panel_A = df[((df.state .== 88) .|| (df.state .== 11)),:]

# Indicate before/after treatment period
# Replace 1990 with 0
panel_A.year[panel_A.year .< 2000] .= 0
# Replace 1991 with 1
panel_A.year[panel_A.year .>= 2000] .= 1

# Indicate treatment group versus control
# Replace 73 with 0
panel_A.state[panel_A.state .== 11] .= 0
# Replace 71 with 1
panel_A.state[panel_A.state .== 88] .= 1

X = hcat(fill(1,nrow(panel_A)), panel_A.year, panel_A.state, (panel_A.year.*panel_A.state))
Y = copy(panel_A.coll)

function DID(X,Y)

    conventional_DID  = X\Y
    X = convert(Matrix{Float64}, X) # NOTE THAT X HAS TO BE A Float64 MATRIX IN ORDER FOR INV() TO WORK!
    # Calculate the residuals
    residuals = Y - X * conventional_DID
    # Estimate the variance of the residuals
    sigma_sq = sum(residuals.^2) / (length(residuals) - length(conventional_DID))
    # Compute the covariance matrix of the coefficients
    cov_beta_hat = sigma_sq * inv(X' * X)
    conventional_DID_SE = sqrt(cov_beta_hat[4,4])

    return [conventional_DID[4], conventional_DID_SE]
end 

Panel_A_results = DID(X,Y)

## OK NOW SUPPOSE THE DATA IS SILOED 
silo0 = copy(panel_A[panel_A.state .== 0,:])
silo1 = copy(panel_A[panel_A.state .== 1,:])

function local_silo_regression(treated_or_control, columns=nothing)
    
    X = convert(Matrix{Float64},hcat(fill(1,nrow(treated_or_control)),treated_or_control.year)) # nrow comes from DataFrames package
    if columns !== nothing 
        columns = Symbol.(columns)
        X = hcat(X, Matrix(treated_or_control[:, columns]))
    end 
    Y = treated_or_control.coll
    beta_hat = X\Y
    gamma = beta_hat[2]

    
    residuals = Y - X * beta_hat
    sigma_sq = sum(residuals.^2) / (length(residuals) - length(beta_hat))

    
    X = convert(Matrix{Float64}, X)
    cov_beta_hat = sigma_sq * inv(X' * X)
    gamma_var = cov_beta_hat[2,2]
    
    return [gamma, gamma_var]
end 

untreated = local_silo_regression(silo0, "year", "coll")
treated = local_silo_regression(silo1, "year", "coll")

function undid(treated, untreated)

    UNDID_ATT = treated[1] - untreated[1]

    UNDID_ATT_SE = sqrt(treated[2] + untreated[2])

    return [UNDID_ATT, UNDID_ATT_SE]
end

undid_result = undid(treated, untreated)

# OK that works fine, might need to figure compatability issues when porting to R, Stata, or Python
# Mostly with respect to the use of the DataFrames packages within Julia


########################### Now I'm going to re-produce results from Panel B #####
X = hcat(fill(1,nrow(panel_A)), panel_A.year, panel_A.state, (panel_A.year.*panel_A.state), panel_A.asian, panel_A.black, panel_A.male)
Y = copy(panel_A.coll)

Panel_B_results = DID(X,Y)

untreated = local_silo_regression(silo0, "year", "coll", columns = ["black", "asian", "male"], intercept = false)
treated = local_silo_regression(silo1, "year", "coll", columns = ["black", "asian", "male"], intercept = false)

undid_result_no_intercept = undid(treated, untreated)
undid_result_with_intercept = undid(treated, untreated)

undid_result_no_intercept
undid_result_with_intercept