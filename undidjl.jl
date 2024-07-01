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
ATTs_2 = Vector{Any}(undef, 1000)
Random.seed!(123)
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
    ATTs_2[i] = (X\Y)[4]
end 

mean(ATTs_2)
