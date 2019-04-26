indef_integral_map = {}




# number of loops
# in: probability vector p, dimension r
# out: minimal (optimal) number of loops
def nr_loops_f(p, r):
	return ((sum((p[k])**(2/3) for k in range(len(p))))**(3/2))**r
	


	

	
	

# prob that NP returns e
	
def prob_NP(dim, det, error_norm, delta,q):
	y = var('y')
	z = var('z')
	f = (1-y^2)**((dim-3)/2)
	g = integrate(f, y)

	b = getB_q_ary(delta,det,dim,q)
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
	

# for Gaussian according to Lindner and Peikert
# s = gaussian parameter, sigma = std deviation
	
	
def p_NP_LP(dim, det, s, delta,q):
	b = getB_q_ary(delta,det,dim,q)
	return RR(getSuccesP(dim,b,s))
	

def myErf(x):
    if x > 100:
        return 1
    return erf(x)

	
## returns the success probability of nearest plane for LWE with gaussian error ##
## input: dimension m, array of GS lengths b, s
def getSuccesP(m,b,s):
    return RR(prod(myErf(b[i]*sqrt(pi)/2/s) for i in range(m)))
	
	



## NEW VERSION	## USES OLD VERSION
## returns a list consisting of the lenghts of the gram-schmidt vectors, considering q-ary structure

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
## returns a list consisting of the lenghts of the gram-schmidt vectors, considering q-ary structure

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
	



