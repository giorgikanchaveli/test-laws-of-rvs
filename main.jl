
include("simulations.jl")


# we check how well does empirical threshold work for hypothesis testing for discrete
# measures with fixed number of atoms
# start_time = time()
# Random.seed!(1234)
# n_atoms = 15
# S = 5
# n_rv = 10
# discr_1 = Discr(n_atoms)
# discr_2 = copy(discr_1)
# θ=0.05
# same = true
# rejecteds = check_emp_thresh_discr("wass",discr_1,discr_2,n_rv,θ,S,same)
# if same
#     println("for probability level $θ we reject same prob m $rejecteds times")
# else
#     println("for probability level $θ we reject diff prob m $rejecteds times")
# end
# time_diff = time() - start_time


# start_time = time()
# Random.seed!(1234)
# n_atoms = 15
# S = 5
# M = 5
# n_rv = 10
# discr_1 = Discr(n_atoms)
# discr_2 = copy(discr_1)
# θ = 0.05
# thresh = emp_threshold("wass",discr_1,discr_2,n_rv,θ, S)
# time_diff = time()-start_time
# measures,rejecteds = check_emp_threshold("wass",thresh, n_rv, S, M)


# we check theoretical bound on the mean distance with n simulated r.v. for each p.m
start_time = time()
Random.seed!(1234)
n_atoms = 15
k = 2.0 # length of interval of sample space
S = 1
n_rv = [5,5,10,10,20,20,20,20]
dist = "wass"
discr_1 = Discr(n_atoms)
discr_2 = Discr(n_atoms)
distances = convergence_dist(dist,discr_1,discr_2,n_rv,S)
bound = mean_bound(n_rv, k, 1)










p_1,p_2 = Discr(15),Discr(15)
d = compute_distance("wass",p_1,p_2,[10],10,1)

p_1,p_2 = Betarv(1.0,1.0)
d = compute_distance("wass",p_1,p_2,[10],10,1)
