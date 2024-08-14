include("methods.jl")



function emptheo_theta_plot(p::PM,q::PM, n_rv::Vector{Int}, S::Int,samemeasures::Bool)
    # for each n we plot empirical and theoretical thresholds per θ

    plots = []

    for n in n_rv
        θs = collect(0.0:0.01:1.0)

        distances = vec(sort(wass_distance(p, q, [n], S, 1), dims = 2))
        log_mean_dist = log(mean(distances))

        emp_thresh = quantile(distances, 1 .- θs, sorted = true)
        theor_thresh = theor_thresh_bobk(p, q, n, θs)

        log_emp_thresh = log.(emp_thresh)
        log_theor_thresh = log.(theor_thresh)
        log_diff = log_theor_thresh .- log_emp_thresh

        p_1 = plot(θs, log_emp_thresh, label = "log empirical thresholds", xlabel = "probability",
                    title = "log thresholds per probability level, n = $(n)", 
                    ylabel = "threshold")
        plot!(p_1, θs, log_theor_thresh, label = "log theoretical thresholds")
        hline!(p_1,[log_mean_dist], label="log mean distance = $(round(log_mean_dist,digits=3))",
                 linestyle=:dash, color=:red)
        vline!(p_1,[0.05],label="θ=0.05->$(round(quantile(distances,0.95,sorted=true),digits=3))",
                 color = :green, linestyle = :dash)
        

        if samemeasures
            mean_bound = 2*J_1_exact(p)/sqrt(n)
            log_mean_bound = log(mean_bound)
            hline!(p_1,[log_mean_bound], label="log mean bound = $(round(log_mean_bound,digits=3))",
                 linestyle=:dash, color=:orange)
        end

        p_2 = plot(θs, log_diff, title = "Difference of Log thresholds",xlabel = "probability",legend=false,
                ylabel = "Difference")

        scatter!(p_2, [0.05, 0.25], [log_diff[6],log_diff[26]])
        annotate!(p_2, (0.05+0.015, log_diff[6], text("$(round(log_diff[6],digits=3))", :black, :left)))
        annotate!(p_2, (0.25+0.015, log_diff[26], text("$(round(log_diff[26],digits=3))", :black, :left)))

        push!(plots, (p_1,p_2))
    end
    return plots
end

function emptheo_n_plot(p::PM,q::PM, n_rv::Vector{Int}, S::Int, θ::Float64)
    # here instead we plot emp and theor thresh per n for fixed θ
    
    distances = sort(wass_distance(p, q, n_rv, S, 1), dims = 2)
    emp_thresh = [quantile(distances[i,:], θ, sorted = true) for i in 1:length(n_rv)]
    theor_thresh = theor_thresh_bobk(p, q, n_rv, θ)

    log_emp_thresh = log.(emp_thresh)
    log_theor_thresh = log.(theor_thresh)

    log_diff = log_theor_thresh .- log_emp_thresh

    pl_1 = plot(n_rv, log_emp_thresh, title = "thresholds, θ=$θ", label = "empirical",
             xlabel = "n", ylabel = "threshold")
    plot!(pl_1, n_rv, log_theor_thresh, title = "thresholds, θ=$θ", label = "theoretical")
    
    pl_2 = plot(n_rv, log_diff, title = "Difference of log thresholds",
             xlabel = "n")

    return (pl_1,pl_2)
end

function save_plots_n(S::Int, p::PM, q::PM, n_rv::Vector{Int}, samemeasure::Bool)
    # plot theoretical and empirical threshodls per n for θ = 0.05 and θ = 0.25


    pl_1 = emptheo_n_plot(p,q,n_rv,S,0.05)
    pl_2 = emptheo_n_plot(p,q,n_rv,S,0.25)

    current_dir = pwd()
    if samemeasure
        folder_path = joinpath(current_dir, "empvs_theor_thresh/pern/same")
        file_path_1 = joinpath(folder_path, "same_$(typeof(p))_0.05.png")
        file_path_2 = joinpath(folder_path, "same_$(typeof(p))_0.25.png")
    else
        folder_path = joinpath(current_dir, "empvs_theor_thresh/pern/diff")
        file_path_1 = joinpath(folder_path, "diff_$(typeof(p))_0.05.png")
        file_path_2 = joinpath(folder_path, "diff_$(typeof(p))_0.25.png")
    end
   
    pl_1 = plot(pl_1..., layout = (1,2),size = (1000,900))
    pl_2 = plot(pl_2..., layout = (1,2),size = (1000,900))
    
    savefig(pl_1,file_path_1)
    savefig(pl_2,file_path_2)
end




Random.seed!(1234)
p = Discr(15, -1.0,1.0)
q = copy(p)
n_rv = collect(100:1000:100000)
t_1 = time()
S = 5
save_plots_n(S,p,q, n_rv, true)
t_1 = time()-t_1



p = Uniformrv(-1.0,1.0)
q = Uniformrv(-1.0,1.0)
n_rv = collect(100:500:1000)
t_2 = time()
S = 4
save_plots_n(S,p,q, n_rv,true)
t_2 = time()-t_2

p = Betarv(1.0,1.0)
q = Betarv(4.0,6.0)
n_rv = collect(100:500:1000)
t_3 = time()
S = 3
save_plots_n(S,p,q, n_rv, false)
t_3 = time()-t_3









function save_plots_theta(S::Int, p::PM,q::PM, same::Bool)
    # 

    n_rv = [10,100,500, 1000,500,10000,100000]

    t = time()
    plots = emptheo_theta_plot(p, q, n_rv, S, same)
    t = time()-t

    for i in 1:length(plots)

        pl= plot(plots[i]...,layout = (1,2), size = (1000,900))
        current_dir = pwd()
        if same
            folder_path = joinpath(current_dir, "empvs_theor_thresh/pertheta/same")
            file_path = joinpath(folder_path, "same_$(typeof(p))_n=$(n_rv[i]).png")
        else
            folder_path = joinpath(current_dir, "empvs_theor_thresh/pertheta/diff")
            file_path = joinpath(folder_path, "diff_$(typeof(p))_n=$(n_rv[i]).png")
        end
        savefig(pl,file_path)
    end
end

# Random.seed!(1234)

# t_1 = time()
# p = Discr(15,-1.0,1.0)
# q = copy(p)
# save_plots_theta(100,p,q,true)
# t_1 = time()-t_1
# # 100 ->

# t_2 = time()
# p = Uniformrv(-1.0,1.0)
# q = Uniformrv(-1.0,1.0)
# save_plots_theta(100,p,q,true)
# t_2 = time()-t_2
# # 100 ->  600 seconds

# t_3 = time()
# p = Betarv(1.0,1.0)
# q_2 = Betarv(4.0,6.0)
# save_plots_theta(100,p,q_2,false)
# t_3 = time()-t_3 
# # 100 -> 300 seconds




















# ar viyeneb




# function plot_thresh_plot(p::Vector{Any})
#     # plots all the plots for 4 different n
#     l = grid(length(p), 1)
#     #l = @layout [a{0.25h}; b{0.25h}; c{0.25h}; d{0.25h}]
#     plot(p...,layout=l,size = (300*length(p),300*length(p)))
# end

# function plot_thresh_hist(h::Vector{Any})
#     # plots all the plots for different n
#     #l = @layout [a; b; c; d; e]
#     l = grid(length(h), 1)

#     plot(h...,layout=l,size = (300*length(h),300*length(h)))
# end

