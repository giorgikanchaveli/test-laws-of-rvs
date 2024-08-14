using Distributions
using Plots
using Random
using QuadGK
using Interpolations

function BM(t::Vector{Float64})

    # given times vector t = t[1]<t[2]<...<t[k]<t[k] = 1, we simulate brownian motion
    # at those times

    # returns vector w where w[i] is bm at t[i]
    k = length(t)
    @assert (t[1]==0) &&(t[k]==1) 
    @assert k>1

    w = zeros(k) # w[i] corresponds to w(t[i])
    w[1] = 0.0
    for i in 2:k
        z = rand(Normal())
        w[i] = sqrt(t[i]-t[i-1])*z + w[i-1]
    end
    return w
end

function BM(t::Vector{Float64},s::Int)

    # given times vector t = t[1]<t[2]<...<t[k]<t[k] = 1, we simulate s brownian motion
    # at those times

    # returns matrix w where w[j,i] is j-th bm at t[i]
    k = length(t)
    @assert (t[1]==0) &&(t[k]==1) 
    @assert k>1

    w = zeros(s,k) 
    w[:,1] .= 0.0
    for j in 1:s
        w[j,:] = BM(t)
    end
    return w
end



function brownian_bridge(t::Vector{Float64}, w::Vector{Float64})
    # t is vector of times
    # w is brownian motion evaluated at times tessting

    # returns standard brownian bridge over [0,1] evaluated at times t
    b = zeros(length(t))
    for i in 1:length(t)
        b[i] = w[i] - t[i]*w[end]
    end
    return b
end


brownian_bridge(t::Vector{Float64}) = brownian_bridge(t, BM(t)) # without given BM

function brownian_bridge(t::Vector{Float64}, w::Matrix{Float64})
    # t is vector of times
    # w matrix of brownian motions


    # returns s standard brownian bridge over [0,1] evaluated at times t
    s = size(w,1) 
    b = zeros(s,length(t))
    for i in 1:s
        b[i,:] = brownian_bridge(t,w[i,:])
    end
    return b
end

brownian_bridge(t::Vector{Float64},s::Int) = brownian_bridge(t,BM(t,s))

function brownian_bridge_series(t::Vector{Float64},k::Int)
    # brownian bridge using series approximation (not good)
    z = randn(k)
    bb = zeros(length(t))
    for i in 1:length(t)
        bb[i] = 0.0
        for j in 1:k
            bb[i] = bb[i] + z[j]* sqrt(2)*sin(j*π*t[i])/(j*π)
        end
    end
    return bb
end


function brownian_bridge_series(t::Vector{Float64},k::Int, s::Int)
    # t is vector of times
    # k is threshold where we stop series
    # s: number of brownian bridge to generate

    # returns s brownian bridges

    bb = zeros(s,length(t))
    for i in 1:s
        bb[i,:] = brownian_bridge_series(t, k)
    end
    return bb
end
    

function brownian_bridge_approx(n::Int)
    # n : number of points for interpolation
    t = 0.0:1/n:1.0
    y = brownian_bridge(collect(t))
    return cubic_spline_interpolation(t,y)
end


function plot_paths(t, x::Matrix{Float64})
    # plots brownian motions
    p = plot(title = "Sample paths of Stochastic process")
    for i in 1:size(x,1)
        plot!(t, x[i,:], xlabel = "times", ylabel = "values",legend = false)
    end
    return p
end


# Random.seed!(1234)

# t = collect(0.0:0.01:1.0)
# w = BM(t,5)


# Random.seed!(1234)
# bb = brownian_bridge(t)
# Random.seed!(1234)
# bbser = brownian_bridge_series(t,100)
# plot(t,bb)
# plot!(t,bbser)



# Random.seed!(1234)
# t = [0.0, 0.5,1.0]
# bb = brownian_bridge(t, )





#integral, error = quadgk(brownian_bridge, 0, 1)




# Random.seed!(1234)
# n = 1000
# t = collect(0.0:1/n:1.0)
# bb_approx = brownian_bridge_approx(n)

# plot(t, bb_approx(t))
# Random.seed!(1234)
# plot!(t,brownian_bridge(t))




# final_integral = 0
# for i in 1:1000 
#     n = 10000
#     bb_approx = brownian_bridge_approx(n)
#     integral,error = quadgk(bb_approx,0.0,1.0)
#     final_integral+=integral
# end
# final_integral = final_integral/1000


# t = time()
# m = 20000
# integrals = []
# for i in 1:m
#     n = 20000
#     bb_approx = brownian_bridge_approx(n)
#     integral,error = quadgk(bb_approx,0.0,1.0)
#     push!(integrals,integral)
# end
# exact = 1/2 * randn(m)
# time_diff = time()-t

# probss = collect(0.0:0.001:1.0)
# y1 = quantile(integrals,probss)
# y2 = quantile(exact,probss)
# plot(probss,y1)
# plot!(probss,y2)



# histogram(integrals,normalize = true)
