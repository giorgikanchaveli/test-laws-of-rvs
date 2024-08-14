
include("methods.jl")

using QuadGK
using Interpolations


discr_1 = Betarv(1.0,1.0)
discr_2 = Betarv(5.0,13.0)

n = 3000

atoms_1,weights_1 = emp_distr_cdf(discr_1,n)
atoms_2,weights_2 = emp_distr_cdf(discr_2,n)

proj_atoms_1,proj_weights_1 = project_grid(atoms_1)
proj_atoms_2,proj_weights_2 = project_grid(atoms_2)

function F_n(x::Float64, atoms::Vector{Float64}, weights::Vector{Float64}, n::Int)
    # returns empirical cdf function
    counts = 0.0
    weights = weights*n
    for i in 1:length(atoms)
        if atoms[i] <= x
            counts += weights[i]
        end
    end
    return counts/sum(weights)
end

f(x::Float64) = F_n(x, atoms_1,weights_1,n)
g(x::Float64) = F_n(x, atoms_2,weights_2,n)
h(x::Float64) = abs(f(x)-g(x))

time_wass = time()
cost = compute_cost(proj_atoms_1,proj_atoms_2,1)
gamma = ExactOptimalTransport.emd(proj_weights_1, proj_weights_2, cost, Tulip.Optimizer())
dist = sum(cost.*gamma)
time_wass = time() - time_wass

time_cdf = time()
int,er = quadgk(h,-1.0,1.0)
time_cdf = time()-time_cdf

