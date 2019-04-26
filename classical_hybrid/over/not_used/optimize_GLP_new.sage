load get_runtime_ternary.sage
load binary_search.sage
load BKZ_cost.sage



# For GLP
# Find optimal attack parameters
# input: ring dimension n, modulus q, list of possible r's to check
# output: map with different r's, correspunding hermite deltas and runtime (in total, with success probability and so on)
#
# EXPECTED ERROR NORM: sqrt(2*(m-r)/3)

	
def find_r_GLP(n, q, k_vals = [6], min_delta = 1.002, max_delta = 1.03):
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
			f_x = RR(bkzcosts(dim,x) - log((nr_loops_new(r,c,dim, det, x, error_norm,q, size_S)*2**15) / (pr_term * prob_NP(dim, det, error_norm, x,q)), 2))
			# PROB_NP MUSS KORRIGIERT WERDEN!!!! DA DIM DET UND ERROR NORM SICH ÄNDERN!!! ODER NICHT???!!!
			return f_x
		
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkzcosts(dim,delta_r) + 1)
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map
 





#
# Auxiliary functions
#


 
# for TRINARY case: c=r/6
	
def find_delta(det, m, r, error_norm, q, size_S):
	c = r/6
	dim = m-r
	def inner_fun(delta):
		return log(BKZ2_ops_LP(delta), 2) - log(nr_loops_new(r,c,dim, det, delta, error_norm,q, size_S)*2**15, 2)
	return find_zero(inner_fun, "delta", 1.001, 1.03, value=0, args = {}, eps = 1, print_steps = true)