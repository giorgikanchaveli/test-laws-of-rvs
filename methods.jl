include("distances_1.jl")
using Plots
using Statistics
using Random

quan = function(arr,α)
    arr = sort(arr)
    index = floor(Int,length(arr)*α)
    if index ==0
        return arr[1]
    else
        return arr[index]
    end
end





function mean_bound(n_rv::Vector{Int}, k::Float64,p::Int)
    # threoretical bound for the expected distance between two prob. measures on interval of
    # length k
    return k ./ (n_rv.^(1/(2*p)))
end


function convergence_dist(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int}, S::Int,p)
    # for each n in n_rv computes empirical mean distance between two p.m

    # dist: distance used for prob. measures
    # pm_i: probability measure
    # n_rv: vector of number of r.v. you want to simulate from each prob. measure
    # S: number of times to compute distance
    # p: for wasserstein p distance
    distances = compute_distance(dist, pm_1, pm_2, n_rv, S, p)
    return vec(mean(distances, dims=2))
end



function emp_threshold(distances::Matrix{Float64}, θ::Float64,sorted_dist::Bool)
    #θ:probability level
    
    # given distances matrix where rows denote number of r.v used to compute
    # empirical distance, we return vector of thresholds 
    len = size(distances)[1]
    thresholds = zeros(len)
    for i in 1:len
        thresholds[i] = quantile(distances[i,:], 1-θ, sorted = sorted_dist)
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

    distances = compute_distance(dist, pm_1, pm_2, n_rv, S, 1)
    return emp_threshold(distances, θ,false)
end


#thresholdebi ar shemimowmebia

function threshold(dist::String, n::Int, θ::Float64, k = 2)
    # dist: distance used for prob measures
    # n: number of sampled random variables
    # θ: probability level for hypothesis testing
    # k: diameter of space

    # returns the optimal threshold for hypothesis tessting
    if dist=="wass"
        c1 = 512*k*sqrt(log(2))/sqrt(n)
        c2 = sqrt((4*(k^2)*log(1/θ))/n)
        c3 = (8*k*log(1/θ))/(3*n)
        c4 = (64*sqrt(2*(k^2)*sqrt(log(2))log(1/θ)))/n^0.75 
        return c1+c2+c3+c4
    elseif dist=="mmd"
        return sqrt(2/n)*(1+sqrt(2*log(1/θ)))
    end
    print("Distance must be either wass or mmd")
    return "Error"
end


function thresholds(dist::String, n::Int, θ::Vector{Float64}, k = 2)
    # dist: distance used for prob measures
    # n: number of sampled random variables
    # θ: probability level for hypothesis testing
    # k: diameter of space

    # returns the optimal thresholds for hypothesis tessting for each prob level θ
    s = length(θ)
    c = zeros(s)
    for i in 1:s
        c[i] = threshold(dist, n, θ[i],k)
    end
    return c
end


        

function hyp_test(dist::String, pm_1::PM, pm_2::PM, n_rv::Int, k::Float64, θ::Float64, S::Int)
    # S times performs hypothesis testing for two same/diff empirical distributions

    # pm_i: discrete probability measures, i.e. atoms and associated weights
    # n_rv: number of r.v. you want to simulate from each prob. measure
    # k: diameter of the sample space
    # θ: probability level
    # S: number of times to perform hyp. testing

    n_rejected = 0
    distances = vec(compute_distance(dist,pm_1,pm_2,[n_rv], S, 1))
    c = threshold(dist,n_rv,θ,k)

    for s in 1:S
        if distances[s] > c
            n_rejected += 1
        end
    end
    return distances, n_rejected/S,c
end


function emp_thresh_prob_plot_hist(dist::String, pm_1::PM, pm_2::PM, n_rv::Vector{Int}, S::Int)
    # returns plots for empirical thresholds per probability level for each n_atoms
    # in n_rv
    distances = sort(compute_distance(dist, pm_1, pm_2, n_rv, S, 1),dims = 2) # length(n_rv)xS matrix
    probs = collect(1.0:-0.01:0.0)
    hists = [] # list for histograms
    plots = [] # list for threshold per probability plots
    for i in 1:length(n_rv)
        thresh_emp = quantile(distances[i,:], probs,sorted=true)
        
        push!(plots, plot(0.0:0.01:1.0, thresh_emp, label = "empirical thresholds",
                             title = "thresholds per probability level, n = $(n_rv[i])", 
                             xlabel = "probability", ylabel = "threshold"))
        mean_dist = mean(distances[i,:])
        vline!(plots[i],[0.05],label="θ=0.05->$(round(quantile(distances[i,:],0.95,sorted=true),digits=3))",
                            color = :green, linestyle = :dash)

        hline!(plots[i],[mean_dist], label="mean distance = $(round(mean_dist,digits=3))",
                 linestyle=:dash, color=:red)

        # histogram
        push!(hists,histogram(distances[i,:],normalize = true, 
                xlabel = "distances",title = "Histogram for distances, n = $(n_rv[i])",legend=false))
        
    end
    return plots,hists
end

function plot_thresh_plot(p::Vector{Any})
    # plots all the plots for different n
    l = @layout [a{0.25h}; b{0.25h}; c{0.25h}; d{0.25h}]
    plot(p...,layout=l,size = (300*length(p),300*length(p)))
end

function plot_thresh_hist(h::Vector{Any})
    # plots all the plots for different n
    l = @layout [a{0.25h}; b{0.25h}; c{0.25h}; d{0.25h}]
    plot(h...,layout=l,size = (300*length(p),300*length(p)))
end

# p_1,p_2 = Discr(5),Discr(5)
# n_rv=[5,10,15]
# S=50
# p,h = emp_thresh_prob_plot_hist("wass",p_1,p_2,n_rv,S)
# plot_thresh_plot(p)

function check_thresh(thresh::Float64)
    measures = [] # this will keep track of which pair of measures are tested
    rejecteds = [] # this will keep track of amount of time we reject null hypothesis for each 
                   # pair of measures
    for i in 1:10
        
        case = rand(1:4)
        if case == 1
            p_1,p_2 = Discr(rand(5:20)),Discr(rand(5:20))
        elseif case == 2
            p_1 = Discr(rand(5:20))
            p_2 = copy(p_1)
        elseif case == 3
            p_1 = Betarv(0.1,1.0)
            p_2 = Betarv(0.5,1.5+rand())
        else
            p_1,p_2 = Betarv(1.0,1.0),Betarv(1.0,1.0)
        end
        push!(measures,(p_1,p_2))
        distances = compute_distance("wass" ,p_1, p_2, [10000], 100, 1)
        push!(rejecteds,sum(distances.>thresh)/100)
    end
    return measures, rejecteds
end


function check_emp_thresh(n_rv::Int, S::Int, θ::Float64)
    p_1 = Discr(15)
    p_2 = copy(p_1)
    thresh = emp_threshold("wass", p_1,p_2,[n_rv],S,θ)[1] # get empirical threshold from p_1 and p_2
    return check_thresh(thresh) 
end

Random.seed!(1234)
m,r = check_emp_thresh(100,10000,0.05)



# m =  2
# t = time()
# for i in 1:m
#     p_1,p_2 = Betarv(1.0,1.0),Betarv(1.0,2.0)
#     compute_distance("wass",p_1,p_2,[10000],1,1)
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