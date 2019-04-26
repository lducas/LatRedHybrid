# sigma = standard deviation (error_norm=sqrt(dim)*sigma)

def create_p_gauss(sigma,tailbound=-1):
	if tailbound == -1:
		tailbound = 14*sigma
	tailbound = floor(tailbound)
	gauss = [exp(-1/2 * (k-tailbound)**2 / ((sigma)**2)) for k in range(2*tailbound+1)]
	s = sum(gauss)
	for k in range(len(gauss)):
		gauss[k] = RR(gauss[k]/s)
	return gauss
