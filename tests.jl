include("distributions.jl")
include("distances_1.jl")
include("methods.jl")


function mean_bound_unif(n::Vector{Int})
    # this bound uses triangular inequality and only can be used for two uniform(0,1) measures
    bounds = zeros(length(n))
    for i in 1:length(n)
        bounds[i] = π/(4*sqrt(n[i]))
    end
    return bounds
end

function thresh_bobkov(p::Betarv, n_rv::Vector{Int}, θ::Float64)
    # p must be Uniform(0,1)
    @assert p.a==1.0 "a must be 1.0"
    @assert p.b==1.0 "b must be 1.0"
    thresholds = zeros(length(n_rv))
    for i in 1:length(n_rv)
        n = n_rv[i]
        c1 = π/(8*sqrt(n)) # before instead of this we had radamacher complexity
        c2 = sqrt(2*log(1/θ)/n)
        c3 = 4*log(1/θ)/(3*n)
        c4 = 2*sqrt((π/8)*2*log(1/θ))/n^0.75
        thresholds[i] = c1+c2+c3+c4
    end
    return thresholds
end


function J_1(p::PM, M, n_rv, a,b)
    # p is pm on [a,b]
    # firstly we simulate r.vs from p 
    # 
    atoms, weights = emp_distr(pm, n_rv)
    weights = weights .* n_rv

    j = 0.0
    for i in 1:M
        u = unif(a,b)
        j += sqrt(F_n(atoms,weights)*(1-F_n(atoms,weights)))
    end
    j = j * (b-a)/M

end





t = time()
p_1,p_2 = Betarv(1.0,1.0), Betarv(1.0,1.0)
S = 15
n_rv = collect(100:100:1000)
emp_c = emp_threshold("wass", p_1,p_2,n_rv,S,0.05)
c = 2*(thresh_bobkov(p_1,n_rv,0.05/2) + thresh_bobkov(p_1,n_rv,0.05/2))
plot(n_rv,c,label = "theoretical threshold", xlabel = "n", title="Theoretical and Empirical thresholds,θ=0.05")
plot!(n_rv, emp_c, label = "empirical threshold")
time_diff = time() - t
# t = time()
# n_rv = collect(10:100:1000)
# S = 10
# p_1, p_2 = Betarv(1.0,1.0), Betarv(1.0,1.0)
# distances = convergence_dist("wass",p_1,p_2,n_rv,S,1)
# bounds = mean_bound_unif(n_rv)
# time_diff = time()-t
# p = plot(n_rv,distances, xlabel = "n",ylabel = "distance",
#         title = "mean distances for Uniform(0,1) measures",label = "empirical mean")
# plot!(p,n_rv,bounds,label="Bound on mean")


