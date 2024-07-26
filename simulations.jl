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

function convergence_dist(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int}, S::Int)
    # for each n_rv computes empirical distance between two p.m

    # dist: distance used for prob. measures
    # pm_i: probability measure
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance

    distances = []
    for n in n_rv
        d = compute_distances(dist,() -> emp_distr(pm_1,n),()->emp_distr(pm_2,n), S, 1)
        push!(distances, mean(d))
    end
    return distances
end

function emp_threshold(dist::String, pm_1::PM,pm_2::PM, n_rv::Int,S::Int)
    # dist: distance used for prob. measures
    # pm_i: probability measure
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance

    # computes empirical distances between pm_1 and om_2 S times and returns
    # 95% quantile
    d = compute_distances(dist,() -> emp_distr(pm_1,n_rv),()->emp_distr(pm_2,n_rv), S, 1)
    return quantile(d, 0.95)
end

function emp_thresholds(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int},S::Int)
    # for each n_rv computes empirical tnreshold which is 95% quantile
    # of distances between pm_1 and pm_2
    
    # dist: distance used for prob. measures 
    # pm_i: probability measure
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance
    c = []
    for n in n_rv
        push!(c, emp_threshold(dist, pm_1, pm_2, n, S))
    end
    return c
end


function threshold_wass(n::Int, θ::Float64, k = 2)
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

function hyp_test(dist::String, pm_1::PM, pm_2::PM, n_rv::Int, k::Float64, θ::Float64, S::Int)
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



pm_discr_1 = Discr(10)
pm_discr_2= Discr(10)
pm_discr_3 = copy(pm_discr_1)

pm_gam_1 = Gammarv(1.0,2.0,5)
pm_gam_2 = Gammarv(2.0,4.0,5)
pm_gam_3 = Gammarv(1.0,2.0,5)

pm_unif_1 = Uniformrv(0.0,1.0)
pm_unif_2 = Uniformrv(0.0,1.0)
S = 10

# emp_1 = ()->emp_distr(pm_1,10)
# emp_2 = ()->emp_distr(pm_2,10)
# emp_3 = ()->emp_distr(pm_3,10)
n_rv = [10,100,100]
dist="mmd"
thresholds_same = emp_thresholds(dist,pm_discr_1,pm_discr_2,n_rv,S)
thresholds_diff = emp_thresholds(dist,pm_discr_1,pm_discr_3,n_rv,S)
# print(thresholds_same)
# print(thresholds_diff)

converg_discr = convergence_dist(dist,pm_discr_1,pm_discr_2,n_rv,10)
converg_gamma = convergence_dist(dist,pm_gam_1,pm_gam_2,n_rv,10)
converg_unif = convergence_dist(dist,pm_unif_1,pm_unif_2,n_rv,10)

# hyp testing:

 distances, n_rejected, c= hyp_test(dist,pm_unif_1,pm_unif_2,100,2.0,0.05,10)

# n_rejected
# distances





