module UNDID

using Statistics
using LinearAlgebra
using DataFrames

export local_silo_regression, undid

function local_silo_regression(silo, post_col, outcome_col; columns = 0, intercept = 0)

    post = silo[!, Symbol.(post_col)]
    outcome = silo[!, Symbol.(outcome_col)]


    if intercept == 0
        X = convert(Matrix{Float64},hcat(1 .- post, post))
    else
        X = convert(Matrix{Float64},hcat(fill(1,nrow(silo)), post))
    end 

    if columns !== 0 
        columns = Symbol.(columns)
        X = hcat(X, Matrix(silo[:, columns]))
    end 

    Y = outcome
    beta_hat = X\Y 

    if intercept == false
        lambda_pre = beta_hat[1]
        lambda_post = beta_hat[2]
        gamma = lambda_post - lambda_pre
    else
        gamma = beta_hat[2]
    end

    residuals = Y - X * beta_hat
    sigma_sq = sum(residuals.^2) / (length(residuals) - length(beta_hat))
    X = convert(Matrix{Float64}, X)
    

    cov_beta_hat = zeros(size(X, 2), size(X, 2))
    try
        cov_beta_hat = sigma_sq * inv(X' * X)
    catch ex 
        if isa(ex, SingularException)
            det_Gram = det(X' * X)
            println("Warning!! Gram matrix (X' * X) is singular (det = $det_Gram), using pseudoinverse instead.")
            cov_beta_hat = sigma_sq * pinv(X' * X)
        else
            println("Unexpected error occurred:", ex)
            rethrow(ex)
        end
    end 

    if intercept == false
        gamma_var = cov_beta_hat[1,1] + cov_beta_hat[2,2] - 2*cov_beta_hat[1,2]
    else
        gamma_var = gamma_var = cov_beta_hat[2,2]
    end 

    return DataFrame(gamma1 = [gamma], gamma1_var = [gamma_var])
end 

function undid(treated, untreated)

    UNDID_ATT = treated.gamma1 - untreated.gamma1

    UNDID_ATT_SE = sqrt((treated.gamma1_var + untreated.gamma1_var)[1])

    return DataFrame(UNDID_ATT = UNDID_ATT, SE_UNDID_ATT = [UNDID_ATT_SE])
end


end 