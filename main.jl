

# t = time()
# n_rv = [100,1000,10000,100000,1000000]
# S = 1 # 100000 shegvidzlia 5 saatshi
# p_discr,h_discr = emp_thresh_prob_plot_hist("wass",discr_1,discr_2,n_rv,S)
# time_diff = time()-t

# plot_thresh_plot(p_discr)


# t = time()
# unif_1,unif_2 = Betarv(1.0,1.0), Betarv(1.0,1.0)
# S = 000
# p_unif,h_unif = emp_thresh_prob_plot_hist("wass",unif_1,unif_2,n_rv,S)
# time_diff = time()-t






# we check how well does empirical threshold work for hypothesis testing for discrete
# measures with fixed number of atoms (15)

# Random.seed!(1234)
# t = time()
# thresh,m,r,sames = check_emp_thresh(10000,1000,0.05)
# time_diff = time()-t




# start_time = time()
# Random.seed!(1234)
# n_atoms = 15
# S = 5m

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












# p_1,p_2 = Betarv(1.0,1.0), Betarv(1.0,1.0)
# n_rv=[100]
# c = 2*(thresh_bobkov(p_1,n_rv,0.05/2) + thresh_bobkov(p_1,n_rv,0.05/2))
