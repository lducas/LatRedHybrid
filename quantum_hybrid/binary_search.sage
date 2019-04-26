def find_zero(fun, x_name, x_min, x_max, value=0, args = {}, eps = 1, print_steps = false):
	if value != 0:
		def inner_fun(**inner_args):
			return fun(**inner_args) - value
		return find_zero(inner_fun, x_name, x_min, x_max, value=0, args = args, eps = eps, print_steps = print_steps)
	args[x_name] = x_min
	f_x_min = RR(fun(**args))
	args[x_name] = x_max
	f_x_max = RR(fun(**args))
#	print f_x_min, f_x_max
	if f_x_max * f_x_min > 0:
#		print "error, fun(x_min) and fun(x_max) have the same sign", RR(f_x_max * f_x_min)
		return -1
	x_med = (x_max+x_min)/2
	args[x_name] = x_med
	f_x_med = RR(fun(**args))
	while abs(f_x_med) > eps:
		if print_steps:
			print RR(x_min), RR(x_med), RR(x_max)
			print "\t" + str(RR(f_x_min)), str(RR(f_x_med)), str(RR(f_x_max))
		if f_x_min.sign() == f_x_med.sign():
			x_min = x_med
			f_x_min = f_x_med
		else:
			x_max = x_med
			f_x_max = f_x_med
		x_med = (x_max+x_min)/2
		args[x_name] = x_med
		f_x_med = RR(fun(**args))
	return RR(x_med)