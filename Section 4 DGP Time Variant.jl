# First generating synthetic data to work with
# Using a similar DGP as presented in section 4.2 of the UNDID paper

using Statistics # This is just to make sure that functions like mean() work
using Random # For setting seeds

rows = 50
columns = 10
ATTs = Vector{Any}(undef, 1000)
Random.seed!(127)
for i in 1:length(ATTs)

    age = repeat(rand(18:45,rows), 1, 10) + repeat(range(0, stop=9, step=1), 1, 50)'
    case1 = randn(rows, columns) + 0.5*age

    case2 = copy(case1)
    case2[end-16:end, end-4:end] .+= 0.1
    Y = vec(copy(case2))

    # Create matrices indicating before/after treatment period and treatment group status
    time_indicator_post = ones(Int, size(case2))  
    time_indicator_post[:, 1:5] .= 0 

    treated_group_indicator = zeros(Int, size(case2))
    treated_group_indicator[end-16:end, :] .= 1
    
    X = hcat(fill(1, 500), vec(time_indicator_post), vec(treated_group_indicator), vec(time_indicator_post .* treated_group_indicator), vec(age))
    ATTs[i] = (X\Y)[4]
end 

# ATT computed via regular DID
mean(ATTs)

ATTs_UNDID = Vector{Any}(undef, 1000)
Random.seed!(127)
for i in 1:length(ATTs_UNDID)
    
    age = repeat(rand(18:45,rows), 1, 10) + repeat(range(0, stop=9, step=1), 1, 50)'
    case1 = randn(rows, columns) + 0.5*age

    case2 = copy(case1)
    case2[end-16:end, end-4:end] .+= 0.1
    Y = vec(copy(case2))

     # Seperate into 'silos' of treated and control
     treated = copy(case2[end-16:end, :])
     control = copy(case2[1:end-17,:])

     treated_age = copy(age[end-16:end,:])
     control_age = copy(age[1:end-17,:])

     treated_post_indicator = [zeros(17, 5) ones(17, 5)]

    # Run regression locally on treated silo
    X = hcat(fill(1,170), vec(treated_post_indicator),vec(treated_age))
    Y = vec(copy(treated))

    treated_gamma = (X\Y)[2]

    # Run regression locally on control silo
    control_post_indicator = [zeros(33, 5) ones(33, 5)]
    
    X = hcat(fill(1,330), vec(control_post_indicator), vec(control_age))
    Y = vec(copy(control))

    control_gamma = (X\Y)[2] 

    ATTs_UNDID[i] = treated_gamma - control_gamma
end

mean(ATTs_UNDID)