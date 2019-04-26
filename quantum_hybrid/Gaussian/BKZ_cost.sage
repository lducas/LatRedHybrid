### BKZ costs

def bkzcosts_sieve(dimension, blocksize):
	return RR(0.292*blocksize + 16.4 + log(8*dimension,2))
	
def bkzcosts_qsieve(dimension, blocksize):
	return RR(0.265*blocksize + 16.4 + log(8*dimension,2))
	
def bkzcosts_enum(dimension, blocksize):
	return RR(0.18728*blocksize*log(blocksize, 2) - 1.0192*blocksize + 16.1 + log(8*dimension,2))
	
def bkzcosts_qenum(dimension, blocksize):
	return RR((0.18728*blocksize*log(blocksize, 2) - 1.0192*blocksize + 16.1)/2  + log(8*dimension,2))

def k_chen(delta):
    k = ZZ(40)

    while delta_0f(2*k) > delta:
        k *= 2
    while delta_0f(k+10) > delta:
        k += 10
    while True:
        if delta_0f(k) < delta:
            break
        k += 1

    return k
	
def delta_0f(k):
    return RR(k/(2*pi*e) * (pi*k)**(1/k))**(1/(2*(k-1)))