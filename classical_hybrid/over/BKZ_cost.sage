gp.read("bkzsim.gp")



#computes the estimated cost of bkz as a function of the blocksize, number of rounds and dimension of the lattice
def bkzoperations_over(bs,dim,rounds):
	return RR(0.187281*bs*log(bs, 2)-1.0192*bs+ log((dim)*rounds,2) +16.1)

	


# finds min block size according to chen's thesis. finds nr of rounds using simulator, and computes runtime
# use for over
# returns log_2(operations)
def bkzcosts_mult_rounds(dim, x):
	k =36
	while RR(((((pi*k)**(1/k))*k)/(2*pi*e))**(1/(2*(k-1))))>RR(x):
		k=k+1
	iterations = ZZ(gp.simulate(dim,min(k,dim),x)[3]) 
	return bkzoperations_over(min(k,dim),dim,iterations)












def minimal_delta(dim):
	return ((((pi*dim)**(1/dim))*dim)/(2*pi*e))**(1/(2*(dim-1)))
	

	