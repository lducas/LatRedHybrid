load attack_implementation/NP.sage
load get_runtime.sage

def get_instance(n,m,q):
	error = vector(Zmod(q), [mod(i,2) for i in range(m)])
	A = matrix(Zmod(q), m, n, [randint(0,q-1) for i in range(m*n)])
	s = vector(Zmod(q), n, [randint(0,q-1) for i in range(n)])
	b = A*s+error
	return A, b

def get_bases(A, block_size=2):
	n=A.ncols()
	m=A.nrows()
	q = A.base_ring().order()
	
	B = A*A[0:n].inverse()
	B = B.change_ring(ZZ)
	B = B.augment(matrix(ZZ, n, m-n).stack(q*identity_matrix(ZZ, m-n)))
	
	B_reduced = B.transpose().BKZ(block_size = block_size).transpose()
	B_reduced_transposed = B_reduced.transpose()
	hash_map = {}
	B_reduced_gs = B_reduced_transposed.gram_schmidt()[0]
	return B_reduced, B_reduced_gs
	
def count_coefficients(B, B_gs, q, n_samples = 1):
	counts = {}
	for i in range(n_samples):
#		print i
		t = vector(ZZ, [randint(0,q-1) for i in range(B.nrows())])
		x = t - NP(B.transpose(), t, B_gs)
		for i in x:
			if counts.keys().count(i) == 0:
				counts[i] = 1
			else:
				counts[i] = counts[i] + 1
	return counts
	
	
	
## simulation part ##
def sample_gs_basis(dim, b):
	M = matrix(dim, dim, [gauss(0,1) for i in range(dim**2)])
	B = M.gram_schmidt()[0]
	basis_vectors = [B[i]*b[i] for i in range(dim)]
	return basis_vectors
	
def sample_point(basis_vectors):
	dim = len(basis_vectors)
	v = sum( (random()-0.5)*vector for vector in basis_vectors ).list()
	return [round(i) for i in v]
	
def simulate_B(basis_vectors,n_samples):
	counts = {}
	for i in range(n_samples):
		x = sample_point(basis_vectors)
#		print i, x
		for i in x:
			if counts.keys().count(i) == 0:
				counts[i] = 1
			else:
				counts[i] = counts[i] + 1
	return counts
def simulate(n,m,q,r,delta,n_samples):
	det = q**(m-n+r)
	dim = m
	global b
	b = getB(delta,det,dim)
	basis_vectors = sample_gs_basis(dim, b)
	return simulate_B(basis_vectors, n_samples)
	
##comparison##
def compare(n,m,q,r, block_size = 2):
	B, B_gs = get_bases(A, block_size)
	A, b = get_instance(n-r,m,q)
	counts_experiments = count_coefficients(B, B_gs, q, n_samples = 100)
	basis_vectors = [x for x in B]
	counts_simulation = simulate_B(basis_vectors, 100)
	print counts_experiments[0], counts_simulation[0]
	return counts_experiments, counts_simulation