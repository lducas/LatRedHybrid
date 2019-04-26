load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage
load norm_product_type_NTRU.sage


# For NTRU product type, f = 1 + pF, where F is of the product form
# Find optimal attack parameters
# input: ring dimension n, modulus q, parameters d_1 = nr of one entries and d_minus_1 = nr of -1 entries,list of possible r's to check
# output: map with different r's, corresponding hermite deltas and runtime (in total)
# 
# EXPECTED ERROR NORM: ...

 
 
def find_r_NTRU_product_form_under(n, q, norm_f_square, d_g, k_vals = [6], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}
	m = 2*n
	det = q**(n)
	
	for r in k_vals:
		error_norm = sqrt(norm_f_square + ((n-r)/n)*(2 * d_g + 1))
		dim = m-r
		c_minus_1 = round((r*d_g)/(2*n))
		c_1 = round((r*(d_g + 1))/(2*n))
		pr_term = p_split_product(2*n-r, n, d_g + 1, d_g, 2*c_1, 2*c_minus_1) # WITHOUT p_NP !!!
		
		# x = delta
		#inner function with fixed r, c, det, dim, q, error_norm
		def f(x):
			p_S = prob_S_NTRU(pr_term, dim, det, error_norm, x,q)
			p_succ = p_S
			size_S = 1
			f_x = RR(bkzcosts_one_round(dim,x) - log((nr_loops_NTRU(r,c_1,c_minus_1,size_S,dim, det, x, error_norm,q)*rt_NP_under(dim)) / (p_succ), 2))
			return f_x
		
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkzcosts_one_round(dim,delta_r) + 1)
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map
 