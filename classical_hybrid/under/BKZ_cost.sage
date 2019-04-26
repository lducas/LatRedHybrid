

#computes the estimated cost of bkz as a function of the blocksize, number of rounds and dimension of the lattice
def bkzoperations(bs,dim,rounds):
	return RR(0.187281*bs*log(bs, 2)-1.0192*bs+ log((dim+1-bs)*rounds,2) +16.1)





# finds min block size according to chen's thesis. sets nr of rounds = 1, and computes runtime
# use for under
# returns log_2(operations)
def bkzcosts_one_round(dim, x):
	k =36
	while RR(((((pi*k)**(1/k))*k)/(2*pi*e))**(1/(2*(k-1))))>RR(x):
		k=k+1
	return bkzoperations(min(k,dim),dim,1)






def minimal_delta(dim):
	return ((((pi*dim)**(1/dim))*dim)/(2*pi*e))**(1/(2*(dim-1)))
	

	