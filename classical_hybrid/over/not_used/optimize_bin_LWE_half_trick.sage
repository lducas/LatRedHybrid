load get_runtime_ternary.sage
load binary_search.sage
load BKZ_cost.sage

    
	
	
def find_r_bin_LWE_half_trick(n, q, m, k_vals = [4], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}
	det = q**(m-n-1)

	for r in k_vals:
		c = r/4
		dim = m-r
		error_norm = sqrt(dim/4)
		### NOT INCLUDING THE -1/2 TRICK!!!!!
		pr_ter = (binomial(r,2*c))/(2**r)
		size_S = 2
		
		# x = delta
		#inner function with fixed r, c, det, dim, q, error_norm
		def f(x):
			f_x = RR(bkzcosts(dim,x) - log(nr_loops_bin(r,c,dim, det, x, error_norm,q, size_S)*2**15 / (pr_ter * prob_NP(dim, det, error_norm, x, q)), 2))
			return f_x
		
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkzcosts(dim,delta_r) + 1)
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map

	
