load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage



# For GLP
# Find optimal attack parameters
# input: ring dimension n, modulus q, list of possible r's to check (multiples of 6)
# output: map with different r's, corresponding hermite deltas and log_2(runtime) (in total)
#
# EXPECTED ERROR NORM: sqrt(2*(m-r)/3)

	
def find_r_GLP_under(n, q, k_vals = [6], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}
	m = 2*n + 1
	det = q**(n)
	
	for r in k_vals:
		error_norm = sqrt(2*(m-r)/3)
		dim = m-r
		c = r/6
		pr_term = prob_terminate_ter(r,c)
		size_S = 2
		
		# x = delta
		#inner function with fixed r, c, det, dim, q, error_norm
		def f(x):
			f_x = RR(bkzcosts_one_round(dim,x) - log((nr_loops_new(r,c,dim, det, x, error_norm,q, size_S)*rt_NP_under(dim)) / (pr_term * prob_NP(dim, det, error_norm, x,q)), 2))
			return f_x
		
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkzcosts_one_round(dim,delta_r) + 1)
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map
 
