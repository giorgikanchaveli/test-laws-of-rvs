# here we given threshold we plot ROC curve in order to measure
# how effective is method to recognize whether two p.meaures are same or not

include("methods.jl")


emp = function(n, S::Int)
    # n: #rv to construct empirical measures

    # from two same discrete measures with 15 atoms, we construct function
    # that takes as an input θ and outputs threshold
    p = Discr(15,-1.0,1.0)
    q = copy(p)
    distances = vec(wass_distance(p, q, [n], S, 1))
    f = function(θ::Float64)
        return quantile(distances, 1-θ)
    end
    return f
end

theor = function(p::PM, q::PM, n::Int, θ::Float64)
    # later this should be changed since p and q are unknown
    return theor_thresh_bobk(p,q, n, θ)
end





Random.seed!(1234)
function roc_discr(n::Int,distr::String)
    # firstly we get empirical threshold
    # then from two same and two different measures
    # we get false positive and true positive rates

    S = 10000

    emp_thresh = emp(n,S) #
    θs = collect(0.0:0.01:1.0)

    p_same = Discr(15,-1.0,1.0)
    q_same = copy(p_same)
    if distr == "discr"
        p_diff = Discr(15, -1.0,1.0)
        q_diff = Discr(15,-1.0,1.0)

    else
        p_diff = Betarv(1.0,1.0)
        q_diff = Betarv(1.3,1.5)
    end

    distances_same = vec(wass_distance(p_same, q_same, [n], S, 1))
    distances_diff = vec(wass_distance(p_diff, q_diff, [n], S, 1))
    α_emp, α_theor = [], []
    β_emp, β_theor = [], []
    for θ in θs
        push!(α_emp, sum(distances_same .> emp_thresh(θ))/S) # false positive for empirical threshold
        push!(α_theor, sum(distances_same .> theor(p_same, q_same, n, θ))/S) # false positive for theoretical
        push!(β_emp, sum(distances_diff .> emp_thresh(θ))/S) # true positive for empirical
        push!(β_theor, sum(distances_diff .> theor(p_diff, q_diff, n, θ))/S) # true positive for theoretical
    end    

    α_β_plot = [plot(α_emp,β_emp,lw = 5,xlims = (0.0,1.0),ylims = (0.0,1.0),label = "empirical",
                title = "ROC for $(distr),n = $n", xlabel = "False positive rate", ylabel = "True positive rate") 
            plot(α_theor,β_theor,lw=5,xlims = (0.0,1.0),ylims = (0.0,1.0),label = "theoretical",
                title = "ROC for $(distr),n = $n", xlabel = "False positive rate", ylabel = "True positive rate")]

    α_plot = [plot(θs, α_emp, lw = 5, xlims = (0.0,1.0),ylims = (0.0,1.0),label = "empirical",
                title = "$(distr),n = $n",ylabel = "False positive rate", xlabel = "θ"),
                plot(θs, α_theor, lw = 5,xlims = (0.0,1.0),ylims = (0.0,1.0),label = "theoretica",
                title = "$(distr),n = $n",ylabel = "False positive rate", xlabel = "θ")]

    β_plot = [plot(θs, β_emp, lw = 4, xlims = (0.0,1.0),ylims = (0.0,1.0),label = "empirical",
                  title = "$(distr),n = $n",ylabel = "True Positive rate", xlabel = "θ"),
                plot(θs, β_theor, lw = 4,xlims = (0.0,1.0),ylims = (0.0,1.0),label = "theoretical",
                title = "$(distr),n = $n",ylabel = "True Positive rate", xlabel = "θ")]

    return α_β_plot,α_plot,β_plot
end


# n = 50
# S = 5
# α_β_plot,α_plot,β_plot = roc_discr(n,"discr")

# plot(α_β_plot...,layout = (1,2))
# plot(α_plot...,layout = (1,2))
# plot(β_plot...,layout = (1,2))



n_rv = [5, 6]
S = 3
for n in n_rv
    α_β_plot,α_plot,β_plot = roc_discr(n,"discr")

    α_β_plot = plot(α_β_plot...,layout = (1,2)) 
    α_plot = plot(α_plot...,layout = (1,2))
    β_plot = plot(β_plot..., layout = (1,2))

    folder_path = joinpath(pwd(), "effectiveness_thresh")
    file_path_1 = joinpath(folder_path, "alpha_beta_discr_n=$(n)")
    file_path_2 = joinpath(folder_path, "alpha_discr_n=$(n)")
    file_path_3 = joinpath(folder_path, "beta_discr_n=$(n)")

    savefig(α_β_plot,file_path_1)
    savefig(α_plot,file_path_2)
    savefig(β_plot,file_path_3)
end


for n in n_rv
    α_β_plot,α_plot,β_plot = roc_discr(n,"Beta")

    α_β_plot = plot(α_β_plot...,layout = (1,2)) 
    α_plot = plot(α_plot...,layout = (1,2))
    β_plot = plot(β_plot..., layout = (1,2))

    folder_path = joinpath(pwd(), "effectiveness_thresh")
    file_path_1 = joinpath(folder_path, "alpha_beta_Betarv_n=$(n)")
    file_path_2 = joinpath(folder_path, "alpha_Betarv_n=$(n)")
    file_path_3 = joinpath(folder_path, "beta_Betarv_n=$(n)")

    savefig(α_β_plot,file_path_1)
    savefig(α_plot,file_path_2)
    savefig(β_plot,file_path_3)
end