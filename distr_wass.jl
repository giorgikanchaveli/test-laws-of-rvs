# here we explore the distribution of wasserstein distance between two same measures

# √(\frac{n}{2})*W_1(μ_n,v_n) => ∫_a^b abs(B(F(x)))dx for μ=v on X=[a,b]


# explore next: 

#       W_p ?


include("methods.jl")
include("brownian_bridge.jl")
using QuadGK
using Interpolations
using Plots


function F_discr(x::Float64, atoms::Vector{Float64}, weights::Vector{Float64})
    # returns empirical cdf function
    f = 0.0
    for i in 1:length(atoms)
        if atoms[i]<=x
            f+=weights[i]
        end
    end
    return f
end


function W(m::Int, p::PM)
    # we get m sample from integral of brownian bridge applied to CDF of p_1
    # this is limiting distribution of wass distance corrected up to 1/√(n)
    
    values = zeros(m)
    for i in 1:m
        bb = brownian_bridge_approx(10000)
        values[i] = quadgk((x::Float64)->abs(bb(CDF(x,p))),-1.0,1.0)[1]
    end
    return values
end

function W_n(m::Int, p_1::PM, p_2::PM, n_rv::Int)
    # given two same probability measure we get m samples from W(p_n^1,p_n^2)

    values = zeros(m)
    for i in 1:m
        d = wass_distance("wass", p_1,p_2,[n_rv],1,1)[1,1]
        d *= √(n_rv/2)
        values[i] = d
    end
    return values
end

function convergence(p::PM,q::PM, m::Int, n_rv::Int)
    # creates histograms for showing convergence in distribution
    # of wasserstein distance between two same pm'S
    
    # this function is only for one n 

    w_n = zeros(m) # contains samples from W(p_n, q_n)
    w = zeros(m) # contains samples from intergal of BB

    for i in 1:m
        w_n[i] = √(n_rv/2)*(wass_distance("wass", p,q,[n_rv],1,1)[1,1])
        bb = brownian_bridge_approx(10000) # brownian bridge
        w[i] = quadgk((x::Float64)->abs(bb(CDF(x,p))),-1.0,1.0)[1] # Integral of bb applied to cdf of p
    end
    # for histograms
    h1 = histogram(w_n, normalize = true, title = "Histogram with normalization",
                  label = "W_n", size=(1000, 500))
    scatter!(h1,[mean(w_n)], [0], label="mean(W_n) = $(round(mean(w_n),digits=2))", marker=:circle, color=:red)
    h2 = histogram(w, normalize = true,title = "Histogram with normalization",
                  label = "W", size=(1000, 500))
    scatter!(h2,[mean(w)], [0], label="mean(W) = $(round(mean(w),digits=2))", marker=:circle, color=:green)
    h = [h1,h2]  
    
    # for quantile plots
    probss = collect(0.0:0.001:1.0)
    q_w_n = quantile(w_n,probss)
    q_w = quantile(w,probss)
    plots = plot(probss,q_w_n,title = "quantile plots", label = "for W_n",
                xlabel = "probabilities", ylabel = "quantiles",size=(1000, 1000))
    plot!(plots, probss,q_w,title = "quantile plots", label = "for W")
    return h,plots
end


function convergence(p::PM,q::PM, m::Int, n_rv::Vector{Int})
    # returns histograms for each n in n_rv
    plots = []
    for n in n_rv
        push!(plots, convergence(p,q,m,n))
    end
    return plots
end


Random.seed!(1234)


p_1 = Uniformrv(-1.0,1.0)
p_2 = Uniformrv(-1.0,1.0)

# for n = 100 -> 370 seconds

# m = 4000


time_diff = time()

n_rv = [50,100,500,1000,10000]
m = 4000
plots = convergence(p_1,p_2,m,n_rv)

time_diff = time()-time_diff

# I need to plot all hists and all qplots and save them


hists = [l[1] for l in plots]
qplots = [l[2] for l in plots]
hists = [plot(h...) for h in hists]
# I need to plot all hists and all qplots

for i in 1:length(n_rv)
    current_dir = pwd()
    folder_path = joinpath(current_dir, "convergence/quantiles")
    file_path = joinpath(folder_path, "unif_quant_n=$(n_rv[i]).png")
    savefig(qplots[i], file_path)
end
for i in 1:length(n_rv)
    current_dir = pwd()
    folder_path = joinpath(current_dir, "convergence/histograms")
    file_path = joinpath(folder_path, "unif_hist_n=$(n_rv[i]).png")
    savefig(hists[i], file_path)
end

















Random.seed!(1234)
p_1 = Betarv(5.0,15.0)
p_2 = Betarv(5.0,15.0)

# for n = 10 -> 74 seconds

# m = 2000


time_diff = time()
n_rv = [50,100,1000,10000,100000]
m = 2000
plots = convergence(p_1,p_2,m,n_rv)

time_diff = time()-time_diff

# # I need to plot all hists and all qplots and save them


hists = [l[1] for l in plots]
qplots = [l[2] for l in plots]
hists = [plot(h...) for h in hists]
#I need to plot all hists and all qplots

for i in 1:length(n_rv)
    current_dir = pwd()
    folder_path = joinpath(current_dir, "convergence_in_distr/quantiles")
    file_path = joinpath(folder_path, "beta_quant_n=$(n_rv[i]).png")
    savefig(qplots[i], file_path)
end
for i in 1:length(n_rv)
    current_dir = pwd()
    folder_path = joinpath(current_dir, "convergence_in_distr/histograms")
    file_path = joinpath(folder_path, "beta_hist_n=$(n_rv[i]).png")
    savefig(hists[i], file_path)
end

