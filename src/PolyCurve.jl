function poly_exponent(qs, N)
    M = length(qs)
    il, jl = N, (N*(M-1)+1)
    dp = rand(0:0, il, jl)

    dp[1,1:M] = qs
    dp[:,1] = qs[1] .^ (1:N)
    
    for j in 2:jl
        qt = qs[min(j, M):-1:2]
        s = sum(dp[1:il-1,max(1,j-M+1):j-1] .* qt', dims=2)
        q = qs[1]
        s[1] += q*dp[1,j]
        dp[2:il, j] = accumulate((x,y)-> q*x+y, s)
    end
    return dp
end

function poly_compose(ps, qs)
    N = length(ps)
    coeffs = poly_exponent(qs, N)
    return ps'*coeffs
end

# Some test code
# const N = 4
# const M = 3
# qs = [rand(0:5, M-1); 1]
# ps = [rand(0:5, N-1); 1]

