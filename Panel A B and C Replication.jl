using StatFiles # For reading .dta
using DataFrames # For dataframes
using Statistics # Part of base Julia
using LinearAlgebra # Also part of base Julia


df = DataFrame(load("C:/Users/Eric Bruce Jamieson/Documents/Dalhousie Work/Merit Data Sent from Matt/merit_data_changed.dta"))

# Note that the 'merit' column indicates 1 once receiving treatment and 0 for pre-treatment periods
df_subset = df[((df.state .== 71) .|| (df.state .== 73)),:]

# Indicate before/after treatment period
df_subset.year[df_subset.year .< 1991] .= 0
df_subset.year[df_subset.year .>= 1991] .= 1

# Indicate treatment group versus control
df_subset.state[df_subset.state .== 73] .= 0
df_subset.state[df_subset.state .== 71] .= 1


########## REPLICATE PANEL A RESULT ##########

X = hcat(fill(1,nrow(df_subset)), df_subset.year, df_subset.state, (df_subset.year.*df_subset.state))
Y = copy(df_subset.coll)

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

df_subset_results = DID(X,Y)

## OK NOW SUPPOSE THE DATA IS SILOED 
silo0 = copy(df_subset[df_subset.state .== 0,:])
silo1 = copy(df_subset[df_subset.state .== 1,:])

untreated = local_silo_regression(silo0, "year", "coll")
treated = local_silo_regression(silo1, "year", "coll")

undid_result = undid(treated, untreated)

########## REPLICATE PANEL B RESULT ##########

X = hcat(fill(1,nrow(df_subset)), df_subset.year, df_subset.state, (df_subset.year.*df_subset.state), df_subset.asian, df_subset.black, df_subset.male)
Y = copy(df_subset.coll)

Panel_B_results = DID(X,Y)

untreated = local_silo_regression(silo0, "year", "coll", columns = ["black", "asian", "male"], intercept = false)
treated = local_silo_regression(silo1, "year", "coll", columns =  ["asian", "black", "male"], intercept = false)

undid_result_no_intercept = undid(treated, untreated)

########## REPLICATE PANEL C RESULT ##########

