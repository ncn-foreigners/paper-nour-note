## functions
### generate counts according to bivariate Bernoulli distribution
function generate_counts(
        nsims::Integer,
        N::Integer, 
        p₁::Real, 
        p₂::Real, 
        α::Real, 
        dependence = "+") 
    
    if dependence == "+"
        p₀₀ = α * (1 - p₁) + (1 - α) * (1 - p₁)* (1 - p₂)
        p₁₀ = (1 - α) * p₁ * (1 - p₂)
        p₀₁ = (1 - α) * (1 - p₁) * p₂
        p₁₁ = α * p₁ + (1 - α) * p₁ * p₂
    elseif dependence == "-"
        p₀₀ = (1 - α) * (1 - p₁) * (1 - p₂)
        p₁₀ = α * p₁ + (1 - α) * p₁ * (1 - p₂)
        p₀₁ = α * (1 - p₁) + (1 - α) * (1 - p₁) * p₂
        p₁₁ = (1 - α) * p₁ * p₂
    end 
    
    res = rand(Multinomial(N, [p₁₁, p₁₀, p₀₁, p₀₀]), nsims)
    
    return res
end

### bias approximation (equations 28-30)
function bias_approximation(;N::Integer, α::Real, 
                            p₁::Real, p₂::Real, 
                            dependence = "+", 
                            bias_type = "Nour")
    if dependence == "+"
        p₀₀ = α * (1 - p₁) + (1 - α) * (1 - p₁)* (1 - p₂)
        p₁₀ = (1 - α) * p₁ * (1 - p₂)
        p₀₁ = (1 - α) * (1 - p₁) * p₂
        p₁₁ = α * p₁ + (1 - α) * p₁ * p₂
    elseif dependence == "-"
        p₀₀ = (1 - α) * (1 - p₁) * (1 - p₂)
        p₁₀ = α * p₁ + (1 - α) * p₁ * (1 - p₂)
        p₀₁ = α * (1 - p₁) + (1 - α) * (1 - p₁) * p₂
        p₁₁ = (1 - α) * p₁ * p₂
    end #end if
    
    # These expressions do not change when dependence type change **
    # only inclusion probablities change.
    if bias_type == "Nour"
        bias_lower = N * (- p₁₁ * (p₁₁ * p₀₀ - p₁₀ * p₀₁) + (p₁₁ - p₀₀) * p₁₀ * p₀₁) / (p₁₁ ^ 2 + p₁₀ * p₀₁)
        bias_mean  = "Not claculated"
        bias_upper = "Not claculated"
    else 
        bias_upper  = N * (p₁₁ + p₁₀ + p₀₁ + sqrt(p₁₀ * p₀₁) - 1)
        bias_upper -= (sqrt(p₁₀ * p₀₁) / 2 + sqrt(p₀₁ / p₁₀) * (1 - p₁₀)/ 4 + sqrt(p₁₀ / p₀₁) * (1 - p₀₁)/ 4) / 2
        # constant term
        bias_lower  = (p₁₁ + p₁₀ + p₀₁ + 2 * p₁₁ * p₁₀ * p₀₁ / (p₁₁ ^ 2 + p₁₀ * p₀₁) - 1) * N
        # variance terms of order N^-2
        bias_lower += 2 * (p₁₁ * (1 - p₁₁) * p₁₁ * p₁₀ * p₀₁ * 
        ((p₁₁ ^ 2 - 3 * p₁₀ * p₀₁) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 4))
        - p₁₀ * (1 - p₁₀) * (p₁₁ ^ 3) * (p₀₁ ^ 2) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 4)
        - p₀₁ * (1 - p₀₁) * (p₁₁ ^ 3) * (p₁₀ ^ 2) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 4)) / (N ^ 2)
        # covariance terms of order N^0
        bias_lower -= 2 * (p₁₁ * p₁₀ * (3 * (p₁₁ ^ 2) * p₀₁ * ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2) - 
        4 * (p₁₁ ^ 4) * p₀₁ * (p₁₁ ^ 2 + p₁₀ * p₀₁)) + 
        p₁₁ * p₀₁ * (3 * (p₁₁ ^ 2) * p₁₀ * ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2) - 
        4 * (p₁₁ ^ 4) * p₁₀ * ((p₁₁ ^ 2) + p₁₀ * p₀₁)) + 
        p₁₀ * p₀₁ * (p₁₁ ^ 3) * (p₁₁ ^ 3 - p₁₀ * p₀₁) * (p₁₁ ^ 2 + p₁₀ * p₀₁)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 4)
        bias_mean = (bias_upper + bias_lower) / 2 # thankfully this is linear
    end #end if
    
    return (bias_lower, bias_mean, bias_upper)
end


### variance approximation (equations 28-30)
function variance_approximation(;N::Integer, α::Real, 
                                 p₁::Real, p₂::Real, 
                                 dependence = "+")
    if dependence == "+"
        p₀₀ = α * (1 - p₁) + (1 - α) * (1 - p₁)* (1 - p₂)
        p₁₀ = (1 - α) * p₁ * (1 - p₂)
        p₀₁ = (1 - α) * (1 - p₁) * p₂
        p₁₁ = α * p₁ + (1 - α) * p₁ * p₂
    elseif dependence == "-"
        p₀₀ = (1 - α) * (1 - p₁) * (1 - p₂)
        p₁₀ = α * p₁ + (1 - α) * p₁ * (1 - p₂)
        p₀₁ = α * (1 - p₁) + (1 - α) * (1 - p₁) * p₂
        p₁₁ = (1 - α) * p₁ * p₂
    end #end if
    
    # These expressions do not change when dependence type changes **
    # only inclusion probablities change.
    
    # this is from variance terms
    var_upper  = p₁₁ * (1 - p₁₁) + p₀₁ * (1 - p₀₁) + p₁₀ * (1 - p₁₀)
    var_upper += p₁₀ * (1 - p₁₀) * ((p₀₁ / p₁₀) ^ .5 + (p₀₁ / p₁₀) / 4)
    var_upper += p₀₁ * (1 - p₀₁) * ((p₁₀ / p₀₁) ^ .5 + (p₁₀ / p₀₁) / 4)
    # this is from covariance temrs
    var_upper -= (2 * p₁₁ * p₀₁ * (1 + ((p₁₀ / p₀₁) ^ .5) / 2) + 
    2 * p₁₁ * p₁₀ * (1 + ((p₀₁ / p₁₀) ^ .5) / 2) +
    2 * p₁₀ * p₀₁ * (5 / 4 + ((p₀₁ / p₁₀) ^ .5 + (p₁₀ / p₀₁) ^ .5) / 2))
    var_upper *= N
    # variance terms
    var_lower  = (p₁₁ * (1 − p₁₁) * (((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) ^ 2) + 
    p₀₁ * (1 − p₀₁) * ((1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) ^ 2) +
    p₁₀ * (1 − p₁₀) * ((1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) ^ 2))
    # covariance terms
    var_lower -= 2 * (p₁₀ * p₀₁ * (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) *
    (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) +
    p₁₁ * p₁₀ * (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    ((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) +
    p₁₁ * p₀₁ * (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    ((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)))
    var_lower *= N
    
    var_mean = (2 * p₁₁ * (1 − p₁₁) * (p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2) + 
    2 * p₀₁ * (1 − p₀₁) * (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * (1 + ((p₀₁ / p₁₀) ^ .5) / 2) +
    2 * p₁₀ * (1 − p₁₀) * (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * (1 + ((p₁₀ / p₀₁) ^ .5) / 2))
    var_mean -= 2 * p₁₀ * p₀₁ * ((1 + ((p₁₀ / p₀₁) ^ .5) / 2) * 
    (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) + 
    (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    (1 + ((p₀₁ / p₁₀) ^ .5) / 2))
    var_mean -= 2 * p₁₁ * p₁₀ * (((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    (1 + ((p₁₀ / p₀₁) ^ .5) / 2) + (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)))
    var_mean -= 2 * p₁₁ * p₀₁ * (((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    (1 + ((p₀₁ / p₁₀) ^ .5) / 2) + (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)))
    
    var_mean *= N / 4
    var_mean += var_lower / 4 + var_upper / 4
    
    return (var_lower, var_mean, var_upper)
end

### variance estimation based on estimated N and known counts: n_11, n_10, n_01 (equations 28-30)

function variance_estimation(;N_est::Real, n₁₁::Real, 
                             n₁₀::Real, n₀₁::Real)
    p₀₀ = (N_est - (n₁₀ + n₀₁ + n₁₁)) / N_est
    p₁₀ = n₁₀ / N_est
    p₀₁ = n₀₁ / N_est
    p₁₁ = n₁₁ / N_est
    N   = N_est
    
    
    var_upper  = p₁₁ * (1 - p₁₁) + p₀₁ * (1 - p₀₁) + p₁₀ * (1 - p₁₀)
    var_upper += p₁₀ * (1 - p₁₀) * ((p₀₁ / p₁₀) ^ .5 + (p₀₁ / p₁₀) / 4)
    var_upper += p₀₁ * (1 - p₀₁) * ((p₁₀ / p₀₁) ^ .5 + (p₁₀ / p₀₁) / 4)
    # this is from covariance temrs
    var_upper -= (2 * p₁₁ * p₀₁ * (1 + ((p₁₀ / p₀₁) ^ .5) / 2) + 
    2 * p₁₁ * p₁₀ * (1 + ((p₀₁ / p₁₀) ^ .5) / 2) +
    2 * p₁₀ * p₀₁ * (5 / 4 + ((p₀₁ / p₁₀) ^ .5 + (p₁₀ / p₀₁) ^ .5) / 2))
    var_upper *= N
    # variance terms
    var_lower  = (p₁₁ * (1 − p₁₁) * (((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) ^ 2) + 
    p₀₁ * (1 − p₀₁) * ((1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) ^ 2) +
    p₁₀ * (1 − p₁₀) * ((1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) ^ 2))
    # covariance terms
    var_lower -= 2 * (p₁₀ * p₀₁ * (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) *
    (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) +
    p₁₁ * p₁₀ * (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    ((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) +
    p₁₁ * p₀₁ * (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    ((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)))
    var_lower *= N
    
    var_mean = (2 * p₁₁ * (1 − p₁₁) * (p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2) + 
    2 * p₀₁ * (1 − p₀₁) * (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * (1 + ((p₀₁ / p₁₀) ^ .5) / 2) +
    2 * p₁₀ * (1 − p₁₀) * (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * (1 + ((p₁₀ / p₀₁) ^ .5) / 2))
    var_mean -= 2 * p₁₀ * p₀₁ * ((1 + ((p₁₀ / p₀₁) ^ .5) / 2) * 
    (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) + 
    (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    (1 + ((p₀₁ / p₁₀) ^ .5) / 2))
    var_mean -= 2 * p₁₁ * p₁₀ * (((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    (1 + ((p₁₀ / p₀₁) ^ .5) / 2) + (1 + 2 * p₀₁ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)))
    var_mean -= 2 * p₁₁ * p₀₁ * (((p₁₁ ^ 4 + 3 * (p₁₀ ^ 2) * (p₀₁ ^ 2)) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)) * 
    (1 + ((p₀₁ / p₁₀) ^ .5) / 2) + (1 + 2 * p₁₀ * (p₁₁ ^ 3) / ((p₁₁ ^ 2 + p₁₀ * p₀₁) ^ 2)))
    
    var_mean *= N / 4
    var_mean += var_lower / 4 + var_upper / 4
    
    return (var_lower, var_mean, var_upper)
end

## simulation study with 100 000 replicates
function sim_study(N, p1, p2, alpha, nsims=100_000, bt = "New")
    tab = generate_counts(nsims, N, p1, p2, alpha, "-")
    lb = (mapslices(x -> 2*x[2]*x[3]*x[1]/(x[2]*x[3]+x[1]^2), tab, dims = 1) .+ sum(tab[1:3,:], dims=1))
    ub = (mapslices(x -> sqrt(x[2]*x[3]), tab, dims = 1) .+ sum(tab[1:3,:], dims=1))
    me = (lb .+ ub) ./2
    res = hcat(lb', me', ub')

     table1 = vcat(mean(res, dims=1) .- N, var(res, dims=1))'

    table2=hcat(
        bias_approximation(N = N, p₁ = p1, p₂ = p2, α = alpha, dependence = "-", bias_type = bt) |> collect,
        variance_approximation(N = N, p₁ = p1, p₂ = p2, α = alpha, dependence = "-") |> collect
        )

    res1=hcat(table1,table2)
end
