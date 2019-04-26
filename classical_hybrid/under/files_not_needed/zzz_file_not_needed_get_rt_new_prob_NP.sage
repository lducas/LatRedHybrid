indef_integral_map = {}

## returns the probability of the following experiment:
## input:
## 		- the interval R = [-r_val, r_val]
## 		- the surface S of the k-dimensional unit sphere S centered around the origin
##  experiment:
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
	if r_val < 1/2:
		int0 = get_indef_integral(k, 0, g)
		int1 = get_indef_integral(k, 1, g)
		h = c_k * (int0.subs(z=r_var-1) - int0.subs(z=-r_var-1) + int1.subs(z=-r_var) - int1.subs(z=r_var-1))
		val = h(r_var=r_val)
	else:
		int2 = get_indef_integral(k, 2, g)
		h = c_k * (int2.subs(z=-r_var) - int2.subs(z=-r_var-1))
		val = h.subs(r_var=r_val)
	val = CC(val)
	if abs(val.imag_part()) > complex_bound:
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
	
## returns a list consisting of the lenghts of the gram-schmidt vectors ##
def getB(delta,det,dim):
#    print RR(dim), RR(det), RR(alpha), b_i
    b = []
    for i in range(dim):
		b_i = delta^(-2*i+dim) * det^(1/dim)
		b.append(b_i)
    return b
	
	
	

# prob NP returns e
	
def prob_NP(dim, det, error_norm, delta):
	y = var('y')
#	z = var('z')
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)

	b = getB(delta,det,dim)
	r_va = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-2*get_prob_NP(dim, r_va[i], g) for i in range(dim))
	
	
def get_prob_NP(k, r_val, g=None):
	complex_bound = 1e-5
	r_val = QQ(RR(r_val))
	if r_val > 1:
		return 0
	
	y = var('y')
	if g == None:
		f = (1-y^2)^((k-3)/2)
		g = integrate(f, y)
	c_k = QQ(RR(1/(beta((k-1)/2, 1/2))))
	
	val = c_k * (g.subs(y = -r_val) - g.subs(y=-1))
	val = CC(val)
	if abs(val.imag_part()) > complex_bound:
		print "error: ", val
		return -1
	return RR(val.real_part())
	

def do_it(dim, det, error_norm, delta):
	r = [delta**(-2*i+dim) * det**(1/dim) / (2*error_norm) for i in range(dim)]
	my_prod = 1
	for i in range(dim):
		if r[i] > 1:
			continue
		y = var('y')
		f = (1-y^2)^((dim-3)/2)
		g = integrate(f, y)
		val = 1 - 2/beta((dim-1)/2, 1/2)*(g.subs(y = -r[i]) - g.subs(y=-1))
		my_prod = my_prod * val
	return RR(my_prod)
	
	
# prob s-admissible	
	
def is_s_admissable(dim, det, error_norm, delta):
	y = var('y')
	z = var('z')
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)

	b = getB(delta,det,dim)
	r_va = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-get_prob(dim, r_va[i], g) for i in range(dim))
	
## returns the number of operatoins BKZ needs to achieve	## 
##   a given hermite delta according to Lindner/Peikert		##
def BKZ2_ops_LP(delta):
	return RR(2**(1.8/log(delta, base=2) - 110)*2.3e9)
	
def is_s_admissable_gs_basis(B):
	y = var('y')
	z = var('z')
	dim=B.nrows()
	error_norm=sqrt(dim/2)
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)
	b = [norm(B[i]) for i in range(dim)]
	r_va = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-2*get_prob(dim, r_va[i], g) for i in range(dim))

	
def is_s_admissable_gs_lengths(b, error_norm):
	y = var('y')
	z = var('z')
	dim=len(b)
	f = (1-y^2)^((dim-3)/2)
	g = integrate(f, y)
	r_va = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-2*get_prob(dim, r_va[i], g) for i in range(dim))
		
		
## Probability that the algorithm terminates
def prob_terminate(r,c):
	return (binomial(r,2*c))/(2**r)


## Probability of a certain error norm
def prob_error_norm(dim,t):
	return (binomial(dim,t))/(2**dim)


	
## Number of loops in the algorithm OLD VERSION
def nr_loops(r,c,prob):
	return binomial(r,c)/(sqrt(prob*binomial(2*c,c)))



## Number of loops in the algorithm NEW VERSION
## r = number of entries in the guessing part
## 2c = number of ones we expect in the last r entries
## dim = dimension of our lattice
## delta = hermite delta
## error_norm = sqrt(dim/4) assumes that we transformed the binary error into something in {+-1/2}^n
def nr_loops_new(r,c,dim, det, delta):
#	print r,c,dim,RR(log(det, 2)),delta
	p_s = is_s_admissable(dim, det, sqrt(dim/4), delta)
#	print "r = ", r, " p_s = ", p_s
	return nr_loops(r,c,p_s)
	
	
def nr_loops_new_print(r,c,dim, det, delta):
#	print r,c,dim,RR(log(det, 2)),delta
	p_s = is_s_admissable(dim, det, sqrt(dim/4), delta)
	print "r = ", r, " p_s = ", RR(p_s)
	return nr_loops(r,c,p_s)
	
	
	
	
def get_bit_hardness(n, q, m, r, delta):
	det = q**(m-n+r)
	c = r/4
	dim = m
	error_norm = sqrt(dim/4)
	p_NP = prob_NP(dim, det, error_norm, delta)
	print "r = ", r, " p_NP = ", RR(p_NP)
	return RR(log((BKZ2_ops_LP(delta) + nr_loops_new_print(r,c,dim, det, delta)*2**16) / (prob_terminate(r,c)*p_NP), 2))
	
	
	
	
