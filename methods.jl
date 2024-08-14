include("distances_1.jl")
using Plots
using Statistics
using Random
using QuadGK

quan = function(arr,α)
    arr = sort(arr)
    index = floor(Int,length(arr)*α)
    if index ==0
        return arr[1]
    else
        return arr[index]
    end
end



function convergence_dist(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int}, S::Int,p)
    # for each n in n_rv computes empirical mean distance between two p.m

    # dist: distance used for prob. measures
    # pm_i: probability measure
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance
    # p: for wasserstein p distance
    distances = wass_distance(dist, pm_1, pm_2, n_rv, S, p)
    return vec(mean(distances, dims=2))
end

function emp_threshold(p::PM, q::PM, S::Int, θ::Float64, n::Int)
    # S: number of simulations of W(p_n,q_n)
    # θ: probability level
    # n: number of r.v's used for constructing empirical p.m.
    
    # given two p.m's returns empirical threshold for fixed θ and n

    d = vec(wass_distance("wass",p,q,[n],S,p)) # distances
    return quantile(d, 1-θ,)
end

function emp_threshold(distances::Matrix{Float64}, θ::Float64)
    #θ: probability level
    
    # given distances matrix where rows denote number of r.v used to compute
    # empirical distance, we return vector of thresholds 
    len = size(distances)[1]
    thresholds = zeros(len)
    d = sort(distances, dims = 2)
    for i in 1:len
        thresholds[i] = quantile(d[i,:], 1-θ, sorted = true)
    end
    return thresholds
end


function emp_threshold(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int},S::Int,θ::Float64)
    # for each n_rv computes empirical tnreshold which is 100*(1-θ)% quantile
    # of distances between pm_1 and pm_2
    
    # dist: distance used for prob. measures 
    # pm_i: probability measure
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance
    # θ: probability level

    distances = wass_distance(dist, pm_1, pm_2, n_rv, S, 1)
    return emp_threshold(distances, θ,false)
end







        
# es funqcia unda shevcvalo threshold given-it
function hyp_test(dist::String, pm_1::PM, pm_2::PM, n_rv::Int, k::Float64, θ::Float64, S::Int)
    # S times performs hypothesis testing for two same/diff empirical distributions

    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # k: diameter of the sample space
    # θ: probability level
    # S: number of times to perform hyp. testing

    n_rejected = 0
    distances = vec(wass_distance(dist,pm_1,pm_2,[n_rv], S, 1))
    c = threshold(dist,n_rv,θ,k)

    for s in 1:S
        if distances[s] > c
            n_rejected += 1
        end
    end
    return distances, n_rejected/S,c
end

function J_1_exact(p::PM)
    # p: pm on [a,b]
    
    # integrates √(F(x)(1-F(x))) on [a,b] assuming that CDF is known
    f = function(x::Float64)
        y  = round(CDF(x, p),digits = 10)
        return sqrt(y*(1-y))
    end
    return quadgk(f,p.a,p.b)[1]
end


function J_1_est(p::PM)
    # p: pm on [a,b]
    
    # integrates √(F(x)(1-F(x))) on [a,b], where F is empirical cdf
    f = function(x::Float64)
        y  = CDF(x, p)
        return sqrt(y*(1-y))
    end
    return quadgk(f,p.a,p.b)[1]
end






function c_n(p::PM, n::Int, θ::Float64)
    j_1 = J_1_exact(p)
    k = p.b-p.a

    c1 = j_1/sqrt(n) # before instead of this we had radamacher complexity
    c2 = sqrt(2*(k^2)*log(1/θ)/n)
    c3 = 4*k*log(1/θ)/(3*n)
    c4 = 2*sqrt((j_1)*2*k*log(1/θ))/n^0.75
    return c1+c2+c3+c4
end

function theor_thresh_bobk(p::PM,q::PM,n::Int,θ::Float64)
    # theoretical threshold given by theorem using bound of bobkhov
    return 2*(c_n(p,n,θ/2) + c_n(q,n,θ/2))
end

function theor_thresh_bobk(p::PM,q::PM,n::Int,θ::Vector{Float64})
    # theoretical threshold given by theorem using bound of bobkhov
    
    m = length(θ)
    t = zeros(m)
    
    for i in 1:m
        t[i] = 2*(c_n(p,n,θ[i]/2) + c_n(q,n,θ[i]/2))
    end
    return t
end
function theor_thresh_bobk(p::PM,q::PM,n_rv::Vector{Int},θ::Float64)
    # theoretical threshold given by theorem using bound of bobkhov
    
    m = length(n_rv)
    t = zeros(m)
    
    for i in 1:m
        t[i] = 2*(c_n(p,n_rv[i],θ/2) + c_n(q,n_rv[i],θ/2))
    end
    return t
end









function thresh_bobkov(p::PM, n_rv::Vector{Int}, θ::Float64)
    # p must be Uniform(0,1)
   
    thresholds = zeros(length(n_rv))
    j_1 = J_1_exact(p)
    for i in 1:length(n_rv)
        n = n_rv[i]
        c1 = j_1/sqrt(n) # before instead of this we had radamacher complexity
        c2 = sqrt(2*log(1/θ)/n)
        c3 = 4*log(1/θ)/(3*n)
        c4 = 2*sqrt((j_1)*2*log(1/θ))/n^0.75
        thresholds[i] = c1+c2+c3+c4
    end
    return thresholds
end

function thresh_bobkov(p::PM, n_rv::Vector{Int}, θ::Vector{Float64})
    # p must be Uniform(0,1)
   
    thresholds = zeros((length(n_rv),length(θ)))
    j_1 = J_1_exact(p)
    for i in 1:length(n_rv)
        for j in 1:length(θ)
            n = n_rv[i]
            c1 = j_1/sqrt(n) # before instead of this we had radamacher complexity
            c2 = sqrt(2*log(1/θ[j])/n)
            c3 = 4*log(1/θ[j])/(3*n)
            c4 = 2*sqrt((j_1)*2*log(1/θ[j]))/n^0.75
            thresholds[i,j] = c1+c2+c3+c4
        end
    end
    return thresholds
end










function check_thresh(thresh::Float64)
    measures = [] # this will keep track of which pair of measures are tested
    rejecteds = [] # this will keep track of amount of time we reject null hypothesis for each 
                   # pair of measures

    sames = [] # denotes whether we test on same measures or not
    for i in 1:10
        
        case = rand(1:4)
        if case == 1
            p_1,p_2 = Discr(rand(5:20)),Discr(rand(5:20))
            push!(sames,false)
        elseif case == 2
            p_1 = Discr(rand(5:20))
            p_2 = copy(p_1)
            push!(sames,true)
        elseif case == 3
            p_1 = Betarv(0.1,1.0)
            p_2 = Betarv(0.5,1.5+rand())
            push!(sames,false)
        else
            p_1,p_2 = Betarv(1.0,1.0),Betarv(1.0,1.0)
            push!(sames,true)
        end
        push!(measures,(p_1,p_2))
        distances = wass_distance("wass" ,p_1, p_2, [10000], 20, 1)
        push!(rejecteds,sum(distances.>thresh)/20)
    end
    return thresh, measures, rejecteds, sames
end


function check_emp_thresh(n_rv::Int, S::Int, θ::Float64)
    p_1 = Discr(15)
    p_2 = copy(p_1)
    thresh = emp_threshold("wass", p_1,p_2,[n_rv],S,θ)[1] # get empirical threshold from p_1 and p_2
    return check_thresh(thresh) 
end

#Random.seed!(1234)
###m,r = check_emp_thresh(100,10000,0.05)




# m =  2
# t = time()
# for i in 1:m
#     p_1,p_2 = Betarv(1.0,1.0),Betarv(1.0,2.0)
#     wass_distance("wass",p_1,p_2,[10000],1,1)
# end
# time_diff = (time()-t)/m
# println("time for computation is $time_diff)")


# time for computing distances between:

#   Discr(20), Discr(20) : 
#               n = 100      -> 0.0037  seconds
#               n = 10000    -> 0.0058  seconds
#               n = 100000   -> 0.0171  seconds
#               n = 1000000  -> 0.133   seconds
#               n = 10000000 -> 1.52    seconds


#  Betarvs: 
#               n = 10       -> 0.0014  seconds
#               n = 100      -> 0.1556  seconds
#               n = 500      -> 10 seconds
#               n = 1000     -> 80 seconds (10 wami axlit)




# aqedan misaxedia 









function mean_bound_unif(n::Vector{Int})
    # this bound uses triangular inequality and only can be used for two uniform(0,1) measures
    bounds = zeros(length(n))
    for i in 1:length(n)
        bounds[i] = π/(4*sqrt(n[i]))
    end
    return bounds
end

function thresh_bobkov(p::Uniformrv, n_rv::Vector{Int}, θ::Float64)
    # p must be Uniform(0,1)

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


function thresh_bobkov(p::PM, n_rv::Vector{Int}, θ::Float64)
    # p must be Uniform(0,1)
   
    thresholds = zeros(length(n_rv))
    j_1 = J_1(p, 2000, 10000, -1.0, 1.0)
    for i in 1:length(n_rv)
        n = n_rv[i]
        c1 = j_1/sqrt(n) # before instead of this we had radamacher complexity
        c2 = sqrt(2*log(1/θ)/n)
        c3 = 4*log(1/θ)/(3*n)
        c4 = 2*sqrt((j_1)*2*log(1/θ))/n^0.75
        thresholds[i] = c1+c2+c3+c4
    end
    return thresholds
end

function thresh_bobkov(p::PM, n_rv::Vector{Int}, θ::Vector{Float64})
    # p must be Uniform(0,1)
   
    thresholds = zeros((length(n_rv),length(θ)))
    j_1 = J_1(p, 2000, 10000, -1.0, 1.0)
    for i in 1:length(n_rv)
        for j in 1:length(θ[j])
            n = n_rv[i]
            c1 = j_1/sqrt(n) # before instead of this we had radamacher complexity
            c2 = sqrt(2*log(1/θ[j])/n)
            c3 = 4*log(1/θ[j])/(3*n)
            c4 = 2*sqrt((j_1)*2*log(1/θ[j]))/n^0.75
            thresholds[i,j] = c1+c2+c3+c4
        end
    end
    return thresholds
end


function J_1(p::PM, M, n_rv, a,b)
    # p is pm on [a,b]
    # firstly we simulate r.vs from p 
    # 
    atoms, weights = emp_distr(p, n_rv)
    weights = weights .* n_rv
    u = Uniform(a, b)
    j = 0.0
    for i in 1:M
        s = rand(u)
        F_n = sum(weights[atoms.<=s])/n_rv
        j += sqrt(F_n*(1-F_n))
    end

    j = j * (b-a)/M
    return j
end

function hyp_test_thresh(dist::String, pm_1::PM, pm_2::PM, n_rv::Int, k::Float64, θ::Float64, S::Int)
    # S times performs hypothesis testing for two same/diff empirical distributions

    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # k: diameter of the sample space
    # θ: probability level
    # S: number of times to perform hyp. testing

    n_rejected = 0
    distances = vec(wass_distance(dist,pm_1,pm_2,[n_rv], S, 1))
    c = 2*(thresh_bobkov(p_1,[n_rv],0.05/2) + thresh_bobkov(p_2,[n_rv],0.05/2))

    for s in 1:S
        if distances[s] > c[1]
            n_rejected += 1
        end
    end
    return distances, n_rejected/S,c
end


# t = time()
# p_1,p_2 = Betarv(1.0,1.0), Betarv(1.5,1.0)
# n_rv=1000000
# d,r, thresh = hyp_test_thresh("wass",p_1,p_2,n_rv,2.0,0.05,50)
# # c = 2*(thresh_bobkov(p_1,n_rv,0.05/2) + thresh_bobkov(p_1,n_rv,0.05/2))
# # thresh,m,r,sames = check_thresh(c[1])
# time_diff = time() - t

# p_unif = Uniformrv(0.0,1.0)
# c_1 = thresh_bobkov(p_1,[n_rv],0.05)
# c_2 = thresh_bobkov(p_unif,[n_rv],0.05)

# p_1 = Betarv(1.0,1.0)
# p_2 = Uniformrv(0.0,1.0)
# n_rv = [1,10,15]
# threshs_1 = thresh_bobkov(p_1,n_rv,0.05)
# threshs_2 =thresh_bobkov(p_2,n_rv,0.05)


# t = time()
# p_1,p_2 = Betarv(1.0,1.0), Betarv(3.0,2.0)
# S = 40
# n_rv = collect(100:1000:100000)
# emp_c = emp_threshold("wass", p_1,p_2,n_rv,S,0.05)
# c = 2*(thresh_bobkov(p_1,n_rv,0.05/2) + thresh_bobkov(p_1,n_rv,0.05/2))
# plot(n_rv,c,label = "theoretical threshold", xlabel = "n", title="Theoretical and Empirical thresholds,θ=0.05")
# plot!(n_rv, emp_c, label = "empirical threshold")
# time_diff = time() - t






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


# ar viyeneb : 
# function threshold_hipm(dist::String, n::Int, θ::Float64, k = 2)
#     # ar viyeneb

#     # dist: distance used for prob measures
#     # n: number of sampled random variables
#     # θ: probability level for hypothesis testing
#     # k: diameter of space

#     # returns the optimal threshold for hypothesis tessting
#     if dist=="wass"
#         c1 = 512*k*sqrt(log(2))/sqrt(n)
#         c2 = sqrt((4*(k^2)*log(1/θ))/n)
#         c3 = (8*k*log(1/θ))/(3*n)
#         c4 = (64*sqrt(2*(k^2)*sqrt(log(2))log(1/θ)))/n^0.75 
#         return c1+c2+c3+c4
#     elseif dist=="mmd"
#         return sqrt(2/n)*(1+sqrt(2*log(1/θ)))
#     end
#     print("Distance must be either wass or mmd")
#     return "Error"
# end




# function thresholds_hipm(dist::String, n::Int, θ::Vector{Float64}, k = 2)
#     # ar viyeneb

#     # dist: distance used for prob measures
#     # n: number of sampled random variables
#     # θ: probability level for hypothesis testing
#     # k: diameter of space

#     # returns the optimal thresholds for hypothesis tessting for each prob level θ
#     s = length(θ)
#     c = zeros(s)
#     for i in 1:s
#         c[i] = threshold_hipm(dist, n, θ[i],k)
#     end
#     return c
# end
