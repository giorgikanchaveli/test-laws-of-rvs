include("distances_1.jl")
using Plots
using Statistics

quan = function(arr,α)
    arr = sort(arr)
    index = floor(Int,length(arr)*α)
    if index ==0
        return arr[1]
    else
        return arr[index]
    end
end

convergence_dist_wass = function(pm_1::Discr, pm_2::Discr, n_rv::Vector{Int}, S::Int)
    # for each n_rv computes empirical distance between two p.m

    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance

    distances = []
    for n in n_rv
        d = compute_distances("wass",() -> emp_distr(pm_1,n),()->emp_distr(pm_2,n), S, 1)
        push!(distances, mean(d))
    end
    return distances
end

emp_threshold_wass = function(pm_1::Discr,pm_2::Discr, n_rv::Int,S::Int)
    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance

    # computes empirical distances between pm_1 and om_2 S times and returns
    # 95% quantile
    d = compute_distances("wass",() -> emp_distr(pm_1,n_rv),()->emp_distr(pm_2,n_rv), S, 1)
    return quantile(d, 0.95)
end

emp_thresholds_wass = function(pm_1::Discr, pm_2::Discr, n_rv::Vector{Int},S::Int)
    # for each n_rv computes empirical tnreshold which is 95% quantile
    # of distances between pm_1 and pm_2
    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance
    c = []
    for n in n_rv
        push!(c, emp_threshold_wass(pm_1, pm_2, n, S))
    end
    return c
end


threshold_wass = function(n::Int, θ::Float64, k = 2)
    # n: number of sampled random variables
    # θ: probability level for hypothesis testing
    # k: diameter of space

    # returns the optimal threshold for hypothesis tessting

    c1 = 512*k*sqrt(log(2))/sqrt(n)
    c2 = sqrt((4*(k^2)*log(1/θ))/n)
    c3 = (8*k*log(1/θ))/(3*n)
    c4 = (64*sqrt(2*(k^2)*sqrt(log(2))log(1/θ)))/n^0.75
    return c1+c2+c3+c4
end

hyp_test = function(dist::String, pm_1::Discr, pm_2::Discr, n_rv::Int, k::Float64, θ::Float64, S::Int)
    # S times performs hypothesis testing for two same/diff empirical distributions

    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # k: diameter of the sample space
    # θ: probability level
    # S: number of times to perform hyp. testing

    n_rejected = 0
    distances = compute_distances(dist,() -> emp_distr(pm_1,n_rv),()->emp_distr(pm_2,n_rv), S, 1)
    if dist=="wass"
        c = threshold_wass(n_rv, θ, k)
    else
        c = threshold_wass(n_rv,θ,k)  # needs to be defined
    end

    for s in 1:S
        if distances[s] > c
            n_rejected += 1
        end
    end
    return distances, n_rejected/S,c
end



pm_1 = Discr(10)
pm_2 = Discr(10)
pm_3 = copy(pm_1)
# emp_1 = ()->emp_distr(pm_1,10)
# emp_2 = ()->emp_distr(pm_2,10)
# emp_3 = ()->emp_distr(pm_3,10)
n_rv = [10,100,1000,10000,20000,30000,100000,1000000]
S = 10
thresholds_same = emp_thresholds_wass(pm_1,pm_2,n_rv,S)
thresholds_diff = emp_thresholds_wass(pm_1,pm_3,n_rv,S)
print(thresholds_same)
print(thresholds_diff)

#convergence_dist_wass(pm_1,pm_2,n_rv,S)
# hyp testing:

distances, n_rejected, c= hyp_test("wass",pm_1,pm_2,100,2.0,0.05,10)

n_rejected
distances





