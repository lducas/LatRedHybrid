# This programm calculates the L: 
# L = multinomial(c_(-k),..,c_(k))*(p+|s|+prod(i from -k to k i!=0)(binomial(2 * list_c[i], list_c[i])))^(-1/2)
# which from the Paper 'On the Hardness of LWE with Binary Error:
# Revisiting the Hybrid Lattice-Reduction and Meet-in-the-Middle Attack'
# new Version page 9 and page 10
# and the p is: p = prod_(i=1)^(m-dim_r) (1-my_integral(dim_r, m, r_i)/(r_i * Beta((m-dim_r-1)/2, 0.5))) 
# the input is: c_(-k),..,c_(k), s, dim_r, m, r_i (for all i in {1, ..., m-dim_r})

# input() gives the information which parameter shoube be input
# test_input() tests, if the input parameters are valid
# test() runs an example test
# my_l() is the final calutation


import sys
from sage.all import *
from sage.combinat.q_analogues import *
from sage.symbolic import *
from sage.symbolic.assumptions import GenericDeclaration


# calculate term (8) from the Paper 'On the Hardness of LWE with Binary Error:
# Revisiting the Hybrid Lattice-Reduction and Meet-in-the-Middle Attack'
# If ri < 1/2, J(ri,m) = integral (integral (1 − y^2)^((m−r-3)/2) dy)dz z from -1 to  z+ri y from -ri−1 to ri−1
#                      + integral (integral (1 − y^2)^((m−r-3)/2) dy)dz z from z-ri to  z+ri y from ri−1 to -ri
# if ri > 1/2, J(ri,m) = integral (integral (1 − y^2)^((m−r-3)/2) dy)dz z from -1 to  z+ri y from -ri−1 to -ri
def my_integral(dim_r, m, r_i):
    var('t, z')
    my_exp = QQ((m - dim_r - 3)/2)
    f = (1.0 - t^2)^my_exp
    if r_i < 0.5 :
        temp = integral(f, t)
        low_1 = temp(t = -1)
        up = temp(t = z + r_i)
        low_2 = temp(t = z - r_i)
        res = QQ(integrate(up - low_1, z)(z = r_i - 1) - integrate(up - low_1, z)( z= -r_i - 1))
        res = QQ(res + integrate(up - low_2, z)(z = -r_i) - integrate(up - low_2, z)(z = r_i - 1))
        return res
    else:
        temp = integral(f, t)
        low_lim = temp(t = -1)
        up_lim = temp(t = z + r_i)
        res = QQ(integrate(up_lim - low_lim, z)(z= - r_i) - integrate(up_lim - low_lim, z)(z= -r_i - 1))
        return res

   
# calculate p = prod_(i=1)^(m-dim_r) (1 - my_integrate(r, m, r_i)/(r_i * Beta((m-dim_r-1)/2, 0.5)))
def my_p(dim_r, m, list_r):
    j = 0
    p = 1.0  
    while j < (m - dim_r) :
        r_i = list_r[j].n()
        be = beta((m-dim_r-1)*0.5, 1/2) 
        my_i = my_integral(dim_r, m, r_i)
        #print 'r_i', r_i, 'be', be, 'my_i', my_i
        p = p * (1.0 - my_i / (r_i * be))
        j = j + 1
    return p

# calculate L = multinomial(c_(-k),..,c_(k))*(p+|s|+prod(i from -k to k i!=0)(binomial(2 * list_c[i], list_c[i])))^(-1/2)
def my_l():
    multnomial_c =  gaussian_multinomial(list_c, 1)
    binomial = 1.0
    i = 0    
    while i < len(list_c) - i - 1:
        binomial = binomial * gaussian_binomial(2 * list_c[i], list_c[i], 1)
        i = i + 1
    i = i + 1
    while i < len(list_c):
        binomial = binomial * gaussian_binomial(2 * list_c[i], list_c[i], 1)
        i = i + 1
    my_l = (multnomial_c / sqrt(my_p(dim_r, m, list_r) *  abs(s) * binomial)).n() 
    print 'binomial',binomial
    print 'my_p:',my_p(dim_r, m, list_r)
    print 'my_l:', my_l
    
# integer (m - dim_r) has to be > 3.    
# list_c: (number of -k, ..., number of k), the sum of the elements of the list is dim_r
# list_r is a list with m-dim_r elements and each r_i must not be 0
# the parameters can be redefined locally
def input():
    print 'Please enter the parameter!'
    print 'integer (m - dim_r) has to be > 3'
    print 'list_c: (number of -k, ..., number of k), the sum of the elements of the list is dim_r'
    print 'list_r is a list with m-r elements and each r_i must not be 0\n'
    global dim_r, m, s, list_c, list_r
    dim_r = 6
    m = 10
    s = 1
    list_c = [3, 2, 1]
    list_r = [0.7, 0.3, 0.2, 0.6]

# test_input() tests, if the input parameters are valid
def test_input():
    if  m < 0 or dim_r < 0 or m - dim_r < 3:
        print 'm, dim_r should be > 0 and m - dim_r should be not < 3, please try again'
    elif GenericDeclaration(x, 'odd').contradicts(x == len(list_c)) :
        print 'list_c should be have 2k+1 elements, please try again'
    elif dim_r != sum(list_c) :
        print 'dim_r != sum(list_c), please try again', dim_r, sum(list_c)
    elif len(list_r) != m-dim_r :
        print 'list_r should be have m - dim_r elements, please try again'
    elif 0 in list_r :
        print 'a ri is 0, please try again'
    else :
        print 'valid input, please start my_l()\n'
   
# a test   
def test():
    global dim_r, m, s, list_c, list_r
    dim_r = 5
    m = 30
    s = 1
    list_c = [2, 2, 1]
    list_r = [5.77437174388516*sqrt(5),
 5.69106233564987*sqrt(5),
 5.60895486899513*sqrt(5),
 5.52803200297957*sqrt(5),
 5.44827664684725*sqrt(5),
 5.36967195641809*sqrt(5),
 5.29220133053043*sqrt(5),
 5.21584840753489*sqrt(5),
 5.14059706183882*sqrt(5),
 5.06643140050062*sqrt(5),
 4.99333575987316*sqrt(5),
 4.92129470229566*sqrt(5),
 4.85029301283327*sqrt(5),
 4.78031569606374*sqrt(5),
 4.71134797291037*sqrt(5),
 4.64337527752072*sqrt(5),
 4.57638325419035*sqrt(5),
 4.51035775433086*sqrt(5),
 4.44528483348183*sqrt(5),
 4.38115074836567*sqrt(5),
 4.31794195398515*sqrt(5),
 4.25564510076267*sqrt(5),
 4.19424703172089*sqrt(5),
 4.13373477970399*sqrt(5),
 4.07409556463900*sqrt(5)]
    my_l()

