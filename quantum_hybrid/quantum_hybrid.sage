load get_runtime_quantum.sage
load binary_search.sage
import numpy as np


# find_r_quantum_improved estimates the cost of the quantum hybrid attack for solving an LWE instance b=As+e mod q with s in Zq^n and e in Zq^m.
# given the LWE parameters and a set of possible attack parameters,
# the function aims at finding the optimal attack parameters and the corresponding complexity of the attack
#####
# using unsuitable parameters as input may result in errors
#####
# we use the embedding of an LWE instance into a uSVP instance that contains the short vector (e,s,1)
#####
# sigma = expected norm per coordinate of e, i.e., sigma such that norm(e) \approx sqrt(m)*sigma
#####
# Possibilities for secret_structure:
# secret_structure = 'random_binary'  --> uniformly random binary secret
# secret_structure = 'hweight_binary'  --> random binary secret  with fixed weight h_weight
# secret_structure = 'hweight_trinary' --> random trinary secret with fixed weight h_weight
# secret_structure = 'random_trinary' --> random trinary secret (weight not fixed)
#####
# bkzmodel in 'sieve', 'qsieve', 'enum', 'qenum', 'NIST_core_enum', 'NIST_core_qsieve'
#####
# nr_rotations = number of rotations of the short vector assumed to be able to be found by the attack
# scaling = re-scaling the secret s to improve the attack
# size_S = size of the search space relative to the maximal search space (0 < size_S <= 1)
# pr_ter = probability of guessing the right number of 1 and -1 entries in the given range
# nr_loops = number of loops in the quantum search
##################

def find_r_quantum_improved(n, q, m, h_weight, sigma, secret_structure = 'random_binary', bkzmodel = 'sieve', nr_rotations = 1, scaling = True, size_S = 1, k_vals = [100], print_results = False, preci = 1):
	global time_map
	time_map = {}
	bitsec, opt_bs, opt_r = 100000000000000000, 0, 0
	det_ = q**(m)
	
	if secret_structure == 'random_binary':
		h_weight = round(n/2)   # expected weight
		
	elif secret_structure == 'random_trinary':
		h_weight = round(n*2/3)   # expected weight
		
	scaling_factor = sigma*sqrt(n / h_weight)

	for r in k_vals:
		nr_ones = round(h_weight/n * r)                     ### expected number of ones
		dim = n+m-r
		d_vec = [1 for j in range(dim)]
		error_norm = sqrt(sigma^2 * m + h_weight-nr_ones)                                             ### without scaling
		
		if scaling == True:                                                                            ### with scaling
			error_norm = sigma * sqrt(dim)
			det_ = q**m * scaling_factor**(n-r)

		s_gauss = error_norm * sqrt(2*pi/dim)
		
		
			
		pr_ter, nr_loops = calc_prter_nrloops(secret_structure, n, r, h_weight, nr_ones, size_S, nr_rotations)
		
		
		
		if bkzmodel == 'sieve':
			delta_r = find_zero(f_sieve, "x", delta_0f(dim), 1.0125, value=0, args = {'d_vec':d_vec, 'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'nr_loops':nr_loops}, eps = 1, print_steps = false)
		elif bkzmodel == 'qsieve':
			delta_r = find_zero(f_q_sieve, "x", delta_0f(dim), 1.0125, value=0, args = {'d_vec':d_vec, 'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'nr_loops':nr_loops}, eps = 1, print_steps = false)
		elif bkzmodel == 'enum':
			delta_r = find_zero(f_enum, "x", delta_0f(dim), 1.0125, value=0, args = {'d_vec':d_vec, 'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'nr_loops':nr_loops}, eps = 1, print_steps = false)
		elif bkzmodel == 'qenum':
			delta_r = find_zero(f_qenum, "x", delta_0f(dim), 1.0125, value=0, args = {'d_vec':d_vec, 'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'nr_loops':nr_loops}, eps = 1, print_steps = false)
		elif bkzmodel == 'NIST_core_enum':
			delta_r = find_zero(f_NIST_core_enum, "x", delta_0f(dim), 1.0125, value=0, args = {'d_vec':d_vec, 'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'nr_loops':nr_loops}, eps = 1, print_steps = false)
		elif bkzmodel == 'NIST_core_qsieve':
			delta_r = find_zero(f_NIST_core_q_sieve, "x", delta_0f(dim), 1.0125, value=0, args = {'d_vec':d_vec, 'dim':dim, 'det_':det_, 's_gauss':s_gauss, 'q':q, 'nr_loops':nr_loops}, eps = 1, print_steps = false)
		
		if delta_r > 0:		
			blocksize_r = k_chen(delta_r)
			delta_r = delta_0f(blocksize_r)
			
			p_NPs_val = p_NPs(d_vec,dim, det_, s_gauss, delta_r,q)
			
			if p_NPs_val * pr_ter <= 0:
				print "p_NPs_val or pr_ter is zero"
				bit_hardness_r = 100000000000000000
			elif bkzmodel == 'sieve':
				bit_hardness_r = RR(bkzcosts_sieve(dim, blocksize_r) + 1 - log(pr_ter * p_NPs_val, 2))
			elif bkzmodel == 'qsieve':
				bit_hardness_r = RR(bkzcosts_qsieve(dim, blocksize_r) + 1 - log(pr_ter * p_NPs_val, 2))
			elif bkzmodel == 'enum':
				bit_hardness_r = RR(bkzcosts_enum(dim, blocksize_r) + 1 - log(pr_ter * p_NPs_val, 2))
			elif bkzmodel == 'qenum':
				bit_hardness_r = RR(bkzcosts_qenum(dim, blocksize_r) + 1 - log(pr_ter * p_NPs_val, 2))
			elif bkzmodel == 'NIST_core_enum':
				bit_hardness_r = RR(bkzcosts_NIST_core_enum(dim, blocksize_r) + 1 - log(pr_ter * p_NPs_val, 2))
			elif bkzmodel == 'NIST_core_qsieve':
				bit_hardness_r = RR(bkzcosts_NIST_core_qsieve(dim, blocksize_r) + 1 - log(pr_ter * p_NPs_val, 2))
				
			time_map[r] = (blocksize_r, bit_hardness_r)
			
			if bit_hardness_r < bitsec:
				bitsec = bit_hardness_r
				opt_r = r
				opt_bs = blocksize_r
				opt_pNPs = p_NPs_val
				opt_pter = pr_ter
			
			if print_results == True:
				print r, time_map[r]
		
#	print "bitsec, opt_bs, opt_r = ", bitsec, opt_bs, opt_r
	return [bitsec, opt_bs, opt_r, m, round(log(size_S, 2))] #, RR(log(opt_pter,2)), RR(log(opt_pNPs,2))



	
# search for optimal parameters using find_r_quantum_improved	
def find_opt_qhybrid(n_, q_, h_, secret_, sig_, bkzm_ = 'NIST_core_enum', m_ = [], S_ = [], r_vals = []):
	best_params = [100000000000]
	for m_S in CartesianProduct(m_, S_):
#		print "m = ", m_S[0], "Size S =2^(-", m_S[1], ")"
		current_params = find_r_quantum_improved(n=n_, q=q_, m=m_S[0], h_weight=h_, sigma=sig_, secret_structure = secret_, bkzmodel = bkzm_, nr_rotations = 1, scaling = True, size_S = 2**(-m_S[1]), k_vals = r_vals, print_results = False, preci = 1)
		
		if current_params[0] < best_params[0]:
			best_params = current_params
			print "currently best parameters found:"
			print "bitsec = ",  best_params[0]
			print "beta = ",  best_params[1]
			print "r = ",  best_params[2]
			print "m = ",  best_params[3]
			print "size_S  = 2^",  best_params[4]
		
	print "best parameters found:"
	print "bitsec = ",  best_params[0]
	print "beta = ",  best_params[1]
	print "r = ",  best_params[2]
	print "m = ",  best_params[3]
	print "size_S  = 2^",  best_params[4]
	
	return best_params
		
		
		
		

# x = delta
#inner functions with fixed r, c, det_, dim, q, error_norm
def f_sieve(x, d_vec, dim, det_, s_gauss,q, nr_loops):
	blocksize_x = min(k_chen(x), dim)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_sieve(dim, blocksize_x) - log((nr_loops * prod(d_vec)*cost_NP(dim)), 2))
	return f_x
	
def f_q_sieve(x, d_vec, dim, det_, s_gauss,q, nr_loops):
	blocksize_x = min(k_chen(x), dim)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_qsieve(dim, blocksize_x) - log((nr_loops * prod(d_vec)*cost_NP(dim)), 2))
	return f_x

def f_enum(x, d_vec, dim, det_, s_gauss,q, nr_loops):
	blocksize_x = min(k_chen(x), dim)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_enum(dim, blocksize_x) - log((nr_loops * prod(d_vec)*cost_NP(dim)), 2))
	return f_x

def f_qenum(x, d_vec, dim, det_, s_gauss,q, nr_loops):
	blocksize_x = min(k_chen(x), dim)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_qenum(dim, blocksize_x) - log((nr_loops * prod(d_vec)*cost_NP(dim)), 2))
	return f_x

def f_NIST_core_enum(x, d_vec, dim, det_, s_gauss,q, nr_loops):
	blocksize_x = min(k_chen(x), dim)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_NIST_core_enum(dim, blocksize_x) - log((nr_loops * prod(d_vec)*cost_NP(dim)), 2))
	return f_x
	
def f_NIST_core_q_sieve(x, d_vec, dim, det_, s_gauss,q, nr_loops):
	blocksize_x = min(k_chen(x), dim)
	delta_x = delta_0f(blocksize_x)
	f_x = RR(bkzcosts_NIST_core_qsieve(dim, blocksize_x) - log((nr_loops * prod(d_vec)*cost_NP(dim)), 2))
	return f_x

### NP costs

def cost_NP(dimension):
#	return 2**(15.1)
	return dimension**2 / (2**(1.06))

	
### BKZ costs

def bkzcosts_sieve(dimension, blocksize):
	return RR(0.292*blocksize + 16.4 + log(8*dimension,2))
	
def bkzcosts_qsieve(dimension, blocksize):
	return RR(0.265*blocksize + 16.4 + log(8*dimension,2))
	
def bkzcosts_enum(dimension, blocksize):
	return RR(0.18728*blocksize*log(blocksize, 2) - 1.0192*blocksize + 16.1 + log(8*dimension,2))
	
def bkzcosts_qenum(dimension, blocksize):
	return RR((0.18728*blocksize*log(blocksize, 2) - 1.0192*blocksize + 16.1)/2  + log(8*dimension,2))

def bkzcosts_NIST_core_enum(dimension, blocksize):
	return RR(0.18728*blocksize*log(blocksize, 2) - 1.0192*blocksize + 16.1)
	
def bkzcosts_NIST_core_qsieve(dimension, blocksize):
	return RR(0.265*blocksize)


def k_chen(delta):
    k = ZZ(40)

    while delta_0f(2*k) > delta:
        k *= 2
    while delta_0f(k+10) > delta:
        k += 10
    while True:
        if delta_0f(k) < delta:
            break
        k += 1

    return k
	
def delta_0f(k):
    return RR(k/(2*pi*e) * (pi*k)**(1/k))**(1/(2*(k-1)))

	

### calculate probability that it is possible to guess the correct vector and the number looos
def calc_prter_nrloops(secret_structure, n, r, h_weight, nr_ones, size_S, nr_rotations):
	if secret_structure == 'random_binary':
		pr_ter_temp = size_S
		nr_loops = sqrt(size_S * 2**r)
		
	elif secret_structure == 'random_trinary': ## random trinary secret (weight not fixed)
		pr_ter_temp = size_S
		nr_loops = sqrt(size_S * 3**r)
		
	elif secret_structure == 'hweight_binary': ## binary secret
		S_i_sizes = [binomial(r, i) for i in [max(0, h_weight-(n-r)) .. min(r, h_weight)]]
		q_i = [binomial(n-r, h_weight-i) / binomial(n, h_weight) for i in [max(0, h_weight-(n-r)) .. min(r, h_weight)]]
		size_M = sum(S_i_sizes)
		pr_ter_temp = RR(0)
		size_S_temp = 0
		sum_temp = RR(0)
		
		while size_S_temp < size_S * size_M:
			index_i_current = np.argmax(q_i)
			size_S_i_current = min(S_i_sizes[index_i_current], ceil(size_S * size_M) - size_S_temp)
			del S_i_sizes[index_i_current]
			q_i_current = q_i[index_i_current]
			del q_i[index_i_current]
			size_S_temp = size_S_temp + size_S_i_current
			pr_ter_temp = pr_ter_temp + RR(q_i_current * size_S_i_current)
			sum_temp = sum_temp + RR(q_i_current**(2/3) * size_S_i_current)
			
		nr_loops = sum_temp**(3/2) / pr_ter_temp
		
	elif secret_structure == 'hweight_trinary': ## random trinary secret with fixed weight h_weight
		S_i_sizes = [binomial(r, i) * 2**i for i in [max(0, h_weight-(n-r)) .. min(r, h_weight)]]
		q_i = [binomial(n-r, h_weight-i) * 2**(-i) / binomial(n, h_weight) for i in [max(0, h_weight-(n-r)) .. min(r, h_weight)]]
		size_M = sum(S_i_sizes)
		pr_ter_temp = RR(0)
		size_S_temp = 0
		sum_temp = RR(0)
		
		while size_S_temp < size_S * size_M:
			index_i_current = np.argmax(q_i)
			size_S_i_current = min(S_i_sizes[index_i_current], ceil(size_S * size_M) - size_S_temp)
			del S_i_sizes[index_i_current]
			q_i_current = q_i[index_i_current]
			del q_i[index_i_current]
			size_S_temp = size_S_temp + size_S_i_current
			pr_ter_temp = pr_ter_temp + RR(q_i_current * size_S_i_current)
			sum_temp = sum_temp + RR(q_i_current**(2/3) * size_S_i_current)
			
		nr_loops = sum_temp**(3/2) / pr_ter_temp
		

	pr_ter = 1-(1-pr_ter_temp)**nr_rotations
		
	return pr_ter, nr_loops
	
	
	
