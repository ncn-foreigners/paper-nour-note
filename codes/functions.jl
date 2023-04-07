## functions
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