indef_integral_map = {}



# runtime of NP conservative (under)
# in: dimension
# out: number of operations
def rt_NP_under(dim):
	return dim/(2^(1.06))
	
	
# runtime of NP current (over)
# in: dimension
# out: number of operations
def rt_NP_over(dim):
	return dim^2/(2^(1.06))





def nr_loops_new(r,c,dim, det, delta, error_norm,q, size_S):
	p_s = is_s_admissable(dim, det, error_norm, delta,q)
	if RR(p_s) <= 0:
		print "p_s negativ or zero, p_s, delta = ", p_s, delta
		p_s = (delta-1) * 10 ** (-500)
	return nr_loops_ter(r,c,p_s, size_S)
	
	
	
	
	
def nr_loops_bin(r,c,dim, det, delta, error_norm,q, size_S):
	p_s = is_s_admissable(dim, det, error_norm, delta,q)
	if RR(p_s) <= 0:
		print "p_s negativ or zero, p_s, delta = ", p_s, delta
		p_s = (delta-1) * 10 ** (-500)
	return binomial(r,c)/(sqrt(p_s*size_S*binomial(2*c,c)))
	
	

## for BLISS
	
def nr_loops_BLISS(r,c_2,c_4,size_S,dim, det, delta, error_norm,q):
	p_s = is_s_admissable(dim, det, error_norm, delta,q)
	if RR(p_s) <= 0:
		print "p_s negativ or zero, p_s, delta = ", p_s, delta
		p_s = (delta-1) * 10 ** (-500)
	return (binomial(r,c_4)*binomial(r-c_4,c_2)*binomial(r-c_4-c_2,c_2)*binomial(r-c_4-c_2-c_2,c_4))/(sqrt(p_s * size_S * (binomial(2*c_2,c_2))**2 * (binomial(2*c_4,c_4))**2))
	
	
	
	
## for NTRU	
	
def nr_loops_NTRU(r,c_1,c_minus_1,size_S,dim, det, delta, error_norm,q):
	p_s = is_s_admissable(dim, det, error_norm, delta,q)
	if RR(p_s) <= 0:
		print "p_s negativ or zero, p_s, delta = ", p_s, delta
		p_s = (delta-1) * 10 ** (-500)
	return (binomial(r,c_1)*binomial(r-c_1,c_minus_1))/(sqrt(p_s * size_S * binomial(2*c_1,c_1) * binomial(2*c_minus_1,c_minus_1)))
	
	
	
	
	
	
	
	
	
	
## for NTRUprime
	
def nr_loops_NTRUprime(r,c_1,c_minus_1,size_S,dim, det, delta, error_norm,q):
	p_s = is_s_admissable(dim, det, error_norm, delta,q)
	if RR(p_s) <= 0:
		print "p_s negativ or zero, p_s, delta = ", p_s, delta
		p_s = (delta-1) * 10 ** (-500)
	return (binomial(r,c_1)*binomial(r-c_1,c_minus_1))/(sqrt(p_s * size_S * binomial(2*c_1,c_1) * binomial(2*c_minus_1,c_minus_1)))
	
	
	
	
	
	
	





## Probability for admissible given GS basis and error norm
## Input: GS basis, error norm
## output: probability


def is_s_admissable_gs_basis(B, error_norm):
	y = var('y')
	z = var('z')
	dim=B.nrows()
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)
	b = [norm(B[i]) for i in range(dim)]
	r = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-get_prob(dim, r[i], g) for i in range(dim))

	
	
	
## Probability for admissible given GS lenghts and error norm
## Input: GS lengths, error norm
## output: probability

	
	
def is_s_admissable_gs_lengths(b, error_norm):
	y = var('y')
	z = var('z')
	dim=len(b)
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)
	r = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-get_prob(dim, r[i], g) for i in range(dim))
	
	
	
	
	
## Probability for admissible given dimension, determinant, delta, and error norm
## Input: dimension, determinant, delta, and error norm
## output: probability


	
def is_s_admissable(dim, det, error_norm, delta,q):
	y = var('y')
	z = var('z')
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)

	b = getB_q_ary(delta,det,dim,q)
	r = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-get_prob(dim, r[i], g) for i in range(dim))	
	
	
	
	
	
	

# prob that NP returns e
	
def prob_NP(dim, det, error_norm, delta,q):
	y = var('y')
	z = var('z')
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)

	b = getB_q_ary(delta,det,dim,q)
	r_va = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-2*get_prob_NP(dim, r_va[i], g) for i in range(dim))
	
	
def get_prob_NP(k, r_val, g=None):
	complex_bound = 1e-5
	r_val = QQ(RR(r_val))
	if RR(r_val) >= 1:
		return 0
	
	y = var('y')
	if g == None:
		f = (1-y^2)^((k-3)/2)
		g = integrate(f, y)
	c_k = QQ(RR(1/(beta((k-1)/2, 1/2))))
	
	val = c_k * (g.subs(y = -r_val) - g.subs(y=-1))
	val = CC(val)
	if RR(abs(val.imag_part())) > complex_bound:
		print "error: ", val
		return -1
	return RR(val.real_part())
	
	
	
## c = 2c in our paper
def p_split(y, N, df, c):
	return binomial(y-N,df-c)*binomial(y-N-(df-c),df-c)*binomial(2*N-y,c)*binomial(2*N-y-c,c)/(binomial(N,df)*binomial(N-df,df))
	
	
## ci = 2ci in our paper
def p_split_product(y, N, d_1, d_minus_1, c_1, c_minus_1):
	return binomial(y-N,d_1-c_1)*binomial(y-N-(d_1-c_1),d_minus_1-c_minus_1)*binomial(2*N-y,c_1)*binomial(2*N-y-c_1,c_minus_1)/(binomial(N,d_1)*binomial(N-d_1,d_minus_1))
	
	
	
	
	
## Probability that the algorithm terminates for TERNARY, WITHOUT p_NP
def prob_terminate_ter(r,c):
	return (binomial(r,2*c)*binomial(r-2*c,2*c))/(3**r)
	
	
	
	
	

## Probability that the algorithm terminates for NTRU, WITHOUT p_NP
def prob_terminate_NTRU(n,r,d_1,d_minus_1):
	c_1 = round((r*d_1)/(2*n))
	c_minus_1 = round((r*d_minus_1)/(2*n))
	c_0 = (r-2*(c_1 + c_minus_1))/2
	p_1 = d_1 / (n)
	p_minus_1 = d_minus_1 / (n)
	p_0 = 1- p_1 - p_minus_1
	return (p_0**(2*c_0) * (p_1**(2*c_1)) * (p_minus_1**(2*c_minus_1)) * binomial(r,2*c_0)*binomial(r-2*c_0,2*c_1))
	
	
	
	
## Probability that the algorithm terminates for NTRU, p_S ONLY, that means for one roatation only
def prob_S_NTRU(pr_term, dim, det, error_norm, delta,q):
	p_S = pr_term * prob_NP(dim, det, error_norm, delta,q)
	return p_S
	
	

	
	
	
## Probability that the algorithm terminates for NTRU_1_pF, WITHOUT p_NP
def prob_terminate_NTRU_1_pF(n,r,d_f):
	c_p = round((r*d_f)/(2*n))
	c_0 = (r-4*(c_p))/2
	p_p = d_f / (n)
	p_0 = 1- 2*p_p
	return (p_0**(2*c_0) * (p_p**(2*c_p))**2 * binomial(r,2*c_0)*binomial(r-2*c_0,2*c_p))	
	

	
	
	

	
## Probability that the algorithm terminates for NTRUprime, WITHOUT p_NP
def prob_terminate_NTRUprime(n,r,d_1):
	c_1 = round((r*d_1)/(2*n)) / 2
	c_0 = (r-4*(c_1))/2
	d_0 = n - d_1
	return (multinomial(2*c_0,2*c_1,2*c_1)*multinomial(d_0-2*c_0,d_1-4*c_1)*2**(d_1-4*c_1))/(multinomial(d_0,d_1)*2**(d_1))
	
	
	
	
	
	
	
	


## Probability that the algorithm terminates for BLISS, WITHOUT p_NP
def prob_terminate_BLISS_new(n,r,d_1,d_2):
	c_2 = round((r*d_1)/(2*n)) / 2
	c_4 = round((r*d_2)/(2*n)) / 2
	c_0 = (r-4*(c_2 + c_4))/2
	d_0 = n - d_1 - d_2
	return (multinomial(2*c_0,2*c_2,2*c_2,2*c_4,2*c_4)*multinomial(d_0-2*c_0,d_1-4*c_2,d_2-4*c_4)*2**(d_1-4*c_2+d_2-4*c_4))/(multinomial(d_0,d_1,d_2)*2**(d_1+d_2))





	
	
## Probability that the algorithm terminates for BLISS, WITHOUT p_NP
def prob_terminate_BLISS(n,r,d_1,d_2):
	c_2 = round((r*d_1)/(2*n)) / 2
	c_4 = round((r*d_2)/(2*n)) / 2
	c_0 = (r-4*(c_2 + c_4))/2
	d_0 = n - d_1 - d_2
	return (multinomial(2*c_0,2*c_2,2*c_2,2*c_4,2*c_4)*multinomial(d_0-2*c_0,floor(d_1/2)-2*c_2,ceil(d_1/2)-2*c_2,floor(d_2/2)-2*c_4,ceil(d_2/2)-2*c_4))/(multinomial(d_0,floor(d_1/2),ceil(d_1/2),floor(d_2/2),ceil(d_2/2)))

	
	
	
	
	
	
## Probability that the algorithm terminates for BLISS, WITHOUT p_NP
def prob_terminate_BLISS_old(n,r,d_1,d_2):
	c_2 = round((r*d_1)/(2*n)) / 2
	c_4 = round((r*d_2)/(2*n)) / 2
	c_0 = (r-4*(c_2 + c_4))/2
	p_2 = d_1 / (2*n)
	p_4 = d_2 / (2*n)
	p_0 = 1- 2*p_2 - 2*p_4
	return (p_0**(2*c_0) * (p_2**(2*c_2))**2 * (p_4**(2*c_4))**2 * binomial(r,2*c_0)*binomial(r-2*c_0,2*c_2)*binomial(r-2*c_0-2*c_2,2*c_2)*binomial(r-2*c_0-2*c_2-2*c_2,2*c_4))



	
	
## Probability that the algorithm terminates for BLISS, p_S ONLY, that means for one roatation only
def prob_S_BLISS(pr_term, dim, det, error_norm, delta,q):
	p_S = pr_term * prob_NP(dim, det, error_norm, delta,q)
	return p_S
	
	
	
## Probability that the algorithm terminates for BLISS, TOTAL
def prob_succ_BLISS(n, r, pr_term, dim, det, error_norm, delta,q):
	p_S = pr_term * prob_NP(dim, det, error_norm, delta,q)
	return (1 - (1 - p_S)**(n-1-r))



	
	
	
	
	
## Number of loops in the algorithm for TERNARY
def nr_loops_ter(r,c,prob, size_S):
	return (binomial(r,c)*binomial(r-c,c))/(sqrt(prob* size_S* binomial(2*c,c)*binomial(2*c,c)))





	
	
##
##
## auxiliary functions
##
##






## returns the probability of the following experiment:
## input:
## 		- the interval R = [-r_val, r_val]
## 		- the surface S of the k-dimensional unit sphere S centered around the origin
## experiment:
## 		- sample t in R
## 		- sample s on S
## 		- success iff t+s_1 is contained in R
def get_prob(k, r_val, g=None):
	complex_bound = 1e-5
	r_val = QQ(RR(r_val))
	y = var('y')
	z = var('z')
	if g == None:
		f = (1-y^2)^((k-3)/2)
		g = integrate(f, y)
	c_k = QQ(RR(1/(r_val*beta((k-1)/2, 1/2))))
	
	val = -1
	if RR(r_val) < 1/2:
		int0 = get_indef_integral(k, 0, g)
		int1 = get_indef_integral(k, 1, g)
		h = c_k * (int0.subs(z=r_var-1) - int0.subs(z=-r_var-1) + int1.subs(z=-r_var) - int1.subs(z=r_var-1))
		val = h(r_var=r_val)
	else:
		int2 = get_indef_integral(k, 2, g)
		h = c_k * (int2.subs(z=-r_var) - int2.subs(z=-r_var-1))
		val = h.subs(r_var=r_val)
	val = CC(val)
	if RR(abs(val.imag_part())) > complex_bound:
		print "error: ", val
		return -1
	return RR(val.real_part())
	
## tool function needed by get_prob(...)
def get_indef_integral(k, index, g):
	if not indef_integral_map.has_key((index, k)):
		r_var = var('r_var')
		assume(r_var > 0)
		if index == 0:
			indef_integral_map[(index, k)] = integral(g.subs(y=z+r_var) - g.subs(y=-1), z)
		if index == 1:
			indef_integral_map[(index, k)] = integral(g.subs(y=z+r_var) - g.subs(y=z-r_var), z)
		if index == 2:
			indef_integral_map[(index, k)] = integral(g.subs(y=z+r_var) - g.subs(y=-1), z)
	return indef_integral_map[index, k]
	
## returns a list consisting of the lenghts of the gram-schmidt vectors
def getB(delta,det,dim):
    b = []
    for i in range(dim):
		b_i = delta**(((-2*dim*i) / (dim-1))+dim) * det^(1/dim)
		b.append(RR(b_i))
    return b
	
	



## NEW VERSION	
## returns a list consisting of the lenghts of the gram-schmidt vectors, CONSIDERING Q-ARY STRUCTURE

def getB_q_ary_new(delta,det,dim,q):
	b = []
	k = floor(1/(2*log(delta, q)))
	if k>dim:
		k=dim
	x=log(det,q)-k/2
	if RR(x)<0:
		print "error"
		return 0
	if RR(x+k)>dim:
		return getB_q_ary(delta,det,dim,q)
	if k<RR((dim-log(det,q))):
		print "k too small in getB_q_ary"
	
	for i in range(dim-k):
		b.append(q)
	
	for i in range(k):
		b_i = RR(delta**(((-2*k*i) / (k-1))+k) * det_new**(1/k))
		b.append(b_i)
	return b
	
	
	
	
	
	
	
	
## OLD VERSION
## returns a list consisting of the lenghts of the gram-schmidt vectors, CONSIDERING Q-ARY STRUCTURE

def getB_q_ary(delta,det,dim,q):
	b = []
	k = floor(sqrt((dim-log(det,q))/(log(delta,q))))
	if k>dim:
		k=dim
	det_new = q**(k-(dim-log(det,q)))
	if k<RR((dim-log(det,q))):
		print "k too small in getB_q_ary"
	
	for i in range(dim-k):
		b.append(q)
	
	for i in range(k):
		b_i = RR(delta**(((-2*k*i) / (k-1))+k) * det_new**(1/k))
		b.append(b_i)
	return b
	

	
