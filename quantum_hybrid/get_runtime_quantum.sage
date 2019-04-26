import numpy as np
indef_integral_map = {}




# number of loops
# in: probability vector p, dimension r
# out: minimal (optimal) number of loops
def nr_loops(p, r):
	return ((sum((p[k])**(2/3) for k in range(len(p))))**(3/2))**r
	


	
def p_NP_LP(dim, det, s, delta,q):
	b = getB(delta,det,dim)
	return RR(getSuccesP(dim,b,s))
	
def p_NPs(d, dim, det, s, delta,q):
	b = getB(delta,det,dim)
	return RR(getSuccesP_NPs(d,dim,b,s))
	

def myErf(x):
    if x > 100:
        return 1
    return erf(x)

## returns the success probability of nearest plane for LWE with gaussian error ##
## input: dimension m, array of GS lengths b, s
def getSuccesP(m,b,s):
    return RR(prod(myErf(b[i]*sqrt(pi)/2/s) for i in range(m)))
	
## returns the success probability of nearest planes for LWE with gaussian error ##
## input: dimension m, array of GS lengths b, s
def getSuccesP_NPs(d,m,b,s):
	bd = [d[i] * b[i] for i in range(m)]
	factor_NPs = RR(sqrt(pi)/2/s)
	probs_bd = [RDF((bd[i]  * factor_NPs)).erf() for i in range(m)]
	return prod(probs_bd)
	
	

# prob that NP returns e
	
def prob_NP(dim, det, error_norm, delta,q):
	y = var('y')
	z = var('z')
	f = (1-y^2)**((dim-3)/2)
	g = integrate(f, y)

	b = getB(delta,det,dim)
	r_va = [b[i] / (2*error_norm) for i in range(dim)]
	return prod(1-2*get_prob_NP(dim, r_va[i], g) for i in range(dim))
	
	
def get_prob_NP(k, r_val, g=None):
	complex_bound = 1e-5
	r_val = QQ(RR(r_val))
	if r_val >= 1:
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
def getB_old(delta,det,dim):
	b = []
	for i in range(dim):
		b_i = delta**(((-2*dim*i) / (dim-1))+dim) * det^(1/dim)
		b.append(RR(b_i))
	return b
	
def getB(delta,det,dim):
	b = []
	constant_factor = delta**(((-2*dim) / (dim-1)))
	b.append(RR(delta**(dim) * det**(1/dim)))
	for i in range(dim-1):
		b.append(RR(b[len(b)-1]*constant_factor))
	return b
	
	



## NEW VERSION	## USES OLD VERSION
## returns a list consisting of the lenghts of the gram-schmidt vectors, !!!!!!!!CONSIDERING Q-ARY STRUCTURE!!!!!!! ##

def getB_q_ary_new(delta,det,dim,q):
	b = []
	k = floor(1/(2*log(delta, q)))
#	print "k = ", k
	if k>dim:
		k=dim
#	print "k = ", k
	x=log(det,q)-k/2
	if x<0:
		print "error"
		return 0
	if x+k>dim:
		return getB_q_ary(delta,det,dim,q)
	if k<(dim-log(det,q)):
		print "k too small in getB_q_ary"
	
	for i in range(dim-k):
		b.append(q)
	
	for i in range(k):
		b_i = RR(delta**(((-2*k*i) / (k-1))+k) * det_new**(1/k))
		b.append(b_i)
	return b
	
	
	
	
	
	
	
	
## OLD VERSION
## returns a list consisting of the lenghts of the gram-schmidt vectors, !!!!!!!!CONSIDERING Q-ARY STRUCTURE!!!!!!! ##

def getB_q_ary(delta,det,dim,q):
#    print RR(dim), RR(det), RR(alpha), b_i
	b = []
	k = floor(sqrt((dim-log(det,q))/(log(delta,q))))
#	print "k = ", k
	if k>dim:
		k=dim
#	print "k = ", k
	det_new = q**(k-(dim-log(det,q)))
	if k<(dim-log(det,q)):
		print "k too small in getB_q_ary"
	
	for i in range(dim-k):
		b.append(q)
	
	for i in range(k):
		b_i = RR(delta**(((-2*k*i) / (k-1))+k) * det_new**(1/k))
		b.append(b_i)
	return b
	

	
## returns the number of operatoins BKZ needs to achieve	## 
##   a given hermite delta according to Lindner/Peikert		##
def BKZ2_ops_LP(delta):
	return RR(2**(1.8/log(delta, base=2) - 110)*2.3e9)
	


