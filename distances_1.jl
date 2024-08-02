using ExactOptimalTransport
using Tulip
using StatsBase
using Plots
using IPMeasures:mmd, GaussianKernel
using KernelFunctions


include("distributions.jl")


function compute_cost(atoms1::Vector{Float64}, atoms2::Vector{Float64}, p::Int)
    # builds cost matrix for atoms. Returns a matrix c of size length(atoms1)*length(atoms2) 
    # where c[i,j] is p-euclidian distance between atoms[i] and atoms[j]

    c = zeros(length(atoms1),length(atoms2))
    for i in 1:length(atoms1)
        for j in 1:length(atoms2)
            c[i,j] = (sum(abs.((atoms1[i] .- atoms2[j]))).^p)^(1/p)
        end
    end
    return c
end

function compute_distance(dist::String, emp_dist_1::Function, emp_dist_2::Function, p::Int)
    # this function is used for empirical distributions
    # if dist == wass returns wasserstein distance and optimal coupling given atoms and associated weights
    # if dist == mmd returns mmd distance
    # p is power for the associated p-euclidian distance

    #emp_dist_i: function s.t. when we call it returns unique atoms and weights associated to it
    
    # simulate emp distributions
    atoms_1, weights_1 = emp_dist_1()
    atoms_2, weights_2 = emp_dist_2()
    #compute disatnce
    if dist == "wass"
        cost = compute_cost(atoms_1,atoms_2,p)
        gamma = ExactOptimalTransport.emd(weights_1, weights_2, cost, Tulip.Optimizer())
        return sum(cost.*gamma)
    elseif dist == "mmd"
        atoms_1,atoms_2 = reshape(atoms_1,(1,length(atoms_1))), reshape(atoms_2,(1,length(atoms_2)))#mmd requries
        return mmd(GaussianKernel(1.0),atoms_1,atoms_2)
    end
end




function compute_distance(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int}, S::Int, p::Int)
    # dist: distance used for prob. measures
    # pm_i: probability measure
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance

    # for each n in n_rv computes S times the distance between prob measures
    # returns length(n_rv)xS matrix
    nrows = length(n_rv)
    d = zeros(nrows,S)
    for i in 1:nrows
        for s in 1:S
            d[i, s] = compute_distance(dist,() -> emp_distr(pm_1,n_rv[i]),()->emp_distr(pm_2,n_rv[i]), p)
        end
    end
    return d
end



# distances_same,n_rejected_same = hyp_test_discrete(5,10000,50,0.05, true)
# distances_diff,n_rejected_diff = hyp_test_discrete(5,10000,50,0.05, false)
# c = threshold(1000,0.05)

# mean_d_same = mean(distances_same)
# mean_d_diff = mean(distances_diff)

# # results: 
# # we simulated two measures with 5 random atoms on [-1,1] and weights
# # if n is very large (~ 10e6), then we same distributions are rejected 0% of time and different
# # distributions are rejected 100% of time.
# # if n is smaller than that we never reject. In the case n=1000, empirical threshold should be around
# # ... and theoretical threshold is ... (test file page 122)
# # histogram looks like normal
# # variance is 10e-6

# histogram(distances_same,normalize = true)

