load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage



# For NTRUprime
# Find optimal attack parameters
# input: ring dimension n = p, modulus q, parameters t,list of possible r's to check
# output: map with different r's, corresponding hermite deltas and runtime (in total)
#
# EXPECTED ERROR NORM: sqrt(n*2/3 + ((n-r)/n)*(2*t))

	
def find_r_NTRUprime_over(p, q, t, k_vals = [6], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}
	n = p
	m = 2*n
	det = q**(n)
	d_1 = 2*t
	
	for r in k_vals:
		error_norm = sqrt(n*2/3 + ((n-r)/n)*(d_1))
		dim = m-r
		c_1 = round((r*d_1)/(4*n))
		pr_term = prob_terminate_NTRUprime(n,r,d_1) # WITHOUT p_NP
		
		# x = delta
		#inner function with fixed r, c, det, dim, q, error_norm
		def f(x):
			p_S = prob_S_BLISS(pr_term, dim, det, error_norm, x,q) # p_NP is included here
			# because they assume n-t good rotations in the paper
			size_S = 2 + 2 * (n - t) * p_S
			f_x = RR(bkzcosts_mult_rounds(dim,x) - log((nr_loops_NTRUprime(r,c_1,c_1,size_S,dim, det, x, error_norm,q)*rt_NP_over(dim)), 2))
			return f_x
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		p_S = prob_S_BLISS(pr_term, dim, det, error_norm, delta_r,q) # p_NP is included here
		p_succ = max(1 - (1 - p_S)**(n-t), p_S)
		bit_hardness_r = RR(bkzcosts_mult_rounds(dim,delta_r) + 1 - log(p_succ, 2))
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map

	
	
	




