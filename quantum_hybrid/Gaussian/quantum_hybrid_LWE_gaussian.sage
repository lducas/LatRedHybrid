load functions.sage
load binary_search.sage
load BKZ_cost.sage
load create_p_gauss.sage

#
# input: secret dimension n, nr of samples m, modulus q, standard deviation sigma, bkzmodel, values for guessing dimension r, bounds on Hermite delta
# 
# output: bitsec against the quantum hybrid attack, optimal guessing dimension r, optimal blocksize
#
#### T = (T_red + T_hyb) / p_succ
	
def find_r_LWE_enum(n, m, q, sigma, bkzmodel = 'enum', k_vals = [6], min_delta = 1.0005, max_delta = 1.01):
	global time_map
	time_map = {}
	bitsec = 100000000
	det_ = q**(m)
	d=m+n
	p = create_p_gauss(sigma,tailbound=-1)
	s_gauss = sigma * (sqrt(2*pi))
	
	for r in k_vals:
		dim = d-r
		error_norm = sigma*sqrt(dim)
		
		nr_loops = nr_loops_f(p, r)
		NP_cost = cost_NP(dim)
		guessing_cost = log(nr_loops * NP_cost, 2)
		
		
		if bkzmodel == 'sieve':
			delta_r = find_zero(f_sieve, "x", min_delta, max_delta, value=0, args = {'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'guessing_cost':guessing_cost}, eps = 1, print_steps = false)
		elif bkzmodel == 'qsieve':
			delta_r = find_zero(f_q_sieve, "x", min_delta, max_delta, value=0, args = {'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'guessing_cost':guessing_cost}, eps = 1, print_steps = false)
		elif bkzmodel == 'enum':
			delta_r = find_zero(f_enum, "x", min_delta, max_delta, value=0, args = {'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'guessing_cost':guessing_cost}, eps = 1, print_steps = false)
		elif bkzmodel == 'qenum':
			delta_r = find_zero(f_qenum, "x", min_delta, max_delta, value=0, args = {'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'guessing_cost':guessing_cost}, eps = 1, print_steps = false)
		
		blocksize_r = k_chen(delta_r)
		delta_r = delta_0f(blocksize_r)
		
		
		if bkzmodel == 'sieve':
			bit_hardness_r = RR(bkzcosts_sieve(dim, blocksize_r) + 1 - log(p_NP_LP(dim, det_, s_gauss, delta_r,q), 2))
		elif bkzmodel == 'qsieve':
			bit_hardness_r = RR(bkzcosts_qsieve(dim, blocksize_r) + 1 - log(p_NP_LP(dim, det_, s_gauss, delta_r,q), 2))
		elif bkzmodel == 'enum':
			bit_hardness_r = RR(bkzcosts_enum(dim, blocksize_r) + 1 - log(p_NP_LP(dim, det_, s_gauss, delta_r,q), 2))
		elif bkzmodel == 'qenum':
			bit_hardness_r = RR(bkzcosts_qenum(dim, blocksize_r) + 1 - log(p_NP_LP(dim, det_, s_gauss, delta_r,q), 2))
			
			
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
		if bit_hardness_r < bitsec:
			bitsec = bit_hardness_r
			opt_r = r
			opt_bs = blocksize_r
			
	print "bitsec, opt_r, opt_bs = ", bitsec, opt_r, opt_bs
	
	
	return [bitsec, opt_r, opt_bs]








# x = delta
#inner functions with fixed r, c, det_, dim, q, error_norm
def f_sieve(x, dim, det_, s_gauss,q, guessing_cost):
	blocksize_x = k_chen(x)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_sieve(dim, blocksize_x) - guessing_cost)
	return f_x
	
def f_q_sieve(x, dim, det_, s_gauss,q, guessing_cost):
	blocksize_x = k_chen(x)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_qsieve(dim, blocksize_x) - guessing_cost)
	return f_x

def f_enum(x, dim, det_, s_gauss,q, guessing_cost):
	blocksize_x = k_chen(x)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_enum(dim, blocksize_x) - guessing_cost)
	return f_x

def f_qenum(x, dim, det_, s_gauss,q, guessing_cost):
	blocksize_x = k_chen(x)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_qenum(dim, blocksize_x) - guessing_cost)
	return f_x


### NP costs

def cost_NP(dimension):
	return dimension**2 / (2**(1.06))