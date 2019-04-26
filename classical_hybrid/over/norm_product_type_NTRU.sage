###
### Experimental norm of product typt
###




def sample_ter_vector(n, d_1, d_minus_1):
	v = vector([0 for i in range(n)])
	count = 0
	while count < d_1:
		i = randint(0, n-1)
		if v[i] == 0:
			v[i] = 1
			count = count+1
	count = 0
	while count < d_minus_1:
		i = randint(0, n-1)
		if v[i] == 0:
			v[i] = -1
			count = count+1
	return v
	
	
	
	

	
def norm_square_experiments(n_samples, n, q, p, d_1, d_2, d_3):
	list = []
	S = PolynomialRing(IntegerModRing(q),'y'); y = S.gen()
	R = S.quotient(y^n - 1, 'x'); x = R.gen()
	for k in range(n_samples):
		a1 = sample_ter_vector(n, d_1, d_1)
		A1 = sum(a1[i] * x**i for i in range(n))
		a2 = sample_ter_vector(n, d_2, d_2)
		A2 = sum(a2[i] * x**i for i in range(n))
		a3 = sample_ter_vector(n, d_3, d_3)
		A3 = sum(a3[i] * x**i for i in range(n))
		f = p * (A1 * A2 + A3)
		#print A1
		#print A2
		#print A3, 
		#print f
		list.append(sum((min(abs(ZZ(f[i])), abs(ZZ(f[i])-q)))**2 for i in range(n)))
	return list
		
	
	
		
	
