gp.read("bkzsim.gp")




#Do we want the current or conservative parameter set
conserv = True


#computes the estimated cost of bkz as a function of the blocksize, number of rounds and dimension of the lattice
def bkzoperations(bs,dim,rounds):
	LogNodes = 0.00405892*bs^2 - 0.337913*bs + 34.9018
	if conserv: LogNodes = 0.000784314*bs^2 + 0.366078*bs - 6.125
	if rounds != 0: return ZZ(round(LogNodes + log(dim*rounds,2) + 7))
	else: return ZZ(round(LogNodes + log(dim,2) + 7))

def bkzoperations2(bs,dim,rounds):
	LogNodes = 0.292*bs
	return ZZ(round(LogNodes + log(bs,2)))

#binary search tree ala schanck to speed up finding optimal BKZ params
def binarysearch(dim,val,low,high):
	if(low >= high -1): return low
	temp = ceil((low + high)/2)
	if(-gp.simulate(dim,temp,val)[2] <=val): return binarysearch(dim,val,temp,high)
	else: return binarysearch(dim,val,low,temp)


def bkzcosts(dim,hermite):
	b = binarysearch(dim,-hermite,60,dim)+1
	beta = gp.simulate(dim,b-1,hermite)[2]
	iterations = ZZ(gp.simulate(dim,b,beta)[3])
	return bkzoperations(b,dim,iterations)
