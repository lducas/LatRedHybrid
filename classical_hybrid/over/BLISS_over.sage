load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage



# For BLISS
# Find optimal attack parameters
# input: ring dimension n, modulus q, parameters d_1 and d_2,list of possible r's to check
# output: map with different r's, correspunding hermite deltas and runtime (in total)
#
# EXPECTED ERROR NORM: sqrt(d_1 + 4*d_2 + ((n-r)/n)*(d_1 + 4*d_2))

	
def find_r_BLISS_over(n, q, d_1, d_2, k_vals = [6], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}
	m = 2*n
	det = q**(n)
	
	for r in k_vals:
		error_norm = sqrt(d_1 + 4*d_2 + ((n-r)/n)*(d_1 + 4*d_2))
		dim = m-r
		c_2 = round((r*d_1)/(4*n))
		c_4 = round((r*d_2)/(4*n))
		pr_term = prob_terminate_BLISS_new(n,r,d_1,d_2) # WITHOUT p_NP
		
		# x = delta
		#inner function with fixed r, c, det, dim, q, error_norm
		def f(x):
			size_S = 1
			f_x = RR(bkzcosts_mult_rounds(dim,x) - log((nr_loops_BLISS(r,c_2,c_4,size_S,dim, det, x, error_norm,q)*rt_NP_over(dim)), 2))
			return f_x
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		p_S = prob_S_BLISS(pr_term, dim, det, error_norm, delta_r,q) # p_NP is included here
		p_succ = p_S
		bit_hardness_r = RR(bkzcosts_mult_rounds(dim,delta_r) + 1 - log(p_succ, 2))
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map

	