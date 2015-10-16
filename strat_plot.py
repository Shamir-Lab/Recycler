import matplotlib.pyplot as plt
import numpy as np


class RefDataType:
	def __init__(self,length,steps,coverage,cv,marker,color,label):
		self.length = length
		self.steps = steps
		self.coverage = coverage
		self.cv = cv
		self.marker = marker
		self.color = color
		self.label = label

	def get_prop(self, prop_str):
		if prop_str == 'length':
			return self.length
		elif prop_str == 'steps':
			return self.steps
		elif prop_str == 'coverage':
			return self.coverage
		elif prop_str == 'cv':
			return self.cv
		elif prop_str == 'marker':
			return self.marker
		elif prop_str == 'color':
			return self.color
		elif prop_str == 'label':
			return self.label
		else:
			return "oh crap!"

def filter_prob_vals(vals, cnts):
	""" replaces cases having no instances
		with last value observed, notes where they occur
		for plotting
	"""
	prob_vals = []
	for ind, v in enumerate(vals):
		if v == 0 and ind != 0 and cnts[ind]==0:
			vals[ind] = vals[ind-1]
			prob_vals.append(ind)
	return vals, prob_vals


def ref_prop_plot(ref, prop, prop_ind, ranges):
	if prop == 'steps':
		x_vals = get_exp_x_mids(ranges[prop_ind],2)
	elif prop == 'coverage':
		x_vals = get_exp_x_mids(ranges[prop_ind],10)
	else:
		x_vals = get_x_mids(ranges[prop_ind])
	num_meas = sum([b[1] for b in ref.get_prop(prop)])
	# s=[a*1000./num_meas for a in [b[1] for b in ref.get_prop(prop)]]
	y_vals = [b[2] for b in ref.get_prop(prop)]
	y_vals, probs = filter_prob_vals(y_vals, [b[1] for b in ref.get_prop(prop)])

	for i in range(len(y_vals)):
		if i == len(y_vals)-1 and prop_ind==3:
			axarr[0][prop_ind].scatter(x_vals[i], y_vals[i], 
			marker=ref.marker, label = ref.get_prop('label'),
			facecolors=ref.color, edgecolors=ref.color, s=100)
		elif i not in probs:
			axarr[0][prop_ind].scatter(x_vals[i], y_vals[i], 
			marker=ref.marker, 
			facecolors=ref.color, edgecolors=ref.color, s=100)	
		else:
			axarr[0][prop_ind].scatter(x_vals[i], y_vals[i], 
				marker=ref.marker, 
				facecolors='none', edgecolors=ref.color, s=1000)

	axarr[0][prop_ind].plot(x_vals, y_vals, c=ref.color)

def ref_counts_plot(ref, prop, prop_ind, ranges):
	x_vals = get_x_mids(ranges[prop_ind])
	row = ref.get_prop(prop)
	instances = np.array([a[1] for a in row])
	instances = np.true_divide(instances, sum(instances))
	axarr[1][prop_ind].scatter(x_vals, instances, 
		c=ref.color, marker=ref.marker, s=100)
	axarr[1][prop_ind].plot(x_vals, instances, c=ref.color)


def get_x_mids(rng):
	return 0.5 * ( np.array(rng[:-1]) + np.array(rng[1:]) ) 

def get_exp_x_mids(rng, base):
	if base == 10:
		vals = np.log10(rng)
	else:
		vals = np.log2(rng)
	return base**get_x_mids(vals)


length_x_rng = [0,4000,8000,12000,16000,20000]
step_x_rng = [1,2,4,8,16,32]
cov_x_rng = [1,10,100,1e3,1e4,1e5]
cv_x_rng = [0, 0.05, 0.1, 0.15, 0.20, 0.25]
ranges = [length_x_rng, step_x_rng, cov_x_rng, cv_x_rng]

# input data - from alignment summary (makefile) output
# rows 0 = length, 1 = steps, 2 = coverage, 3 = CV
# plotting rows 4 = marker, 5 = color
ref_100 = RefDataType(
	[(40, 46, 0.87), (19, 21, 0.9), (3, 6, 0.5), (4, 6, 0.67), (6, 7, 0.86)],
	[(59, 59, 1.0), (6, 11, 0.55), (3, 10, 0.3), (3, 4, 0.75), (1, 2, 0.5)],
	[(4, 4, 1.0), (39, 49, 0.8), (22, 26, 0.85), (7, 7, 1.0), (0, 0, 0)],
	[(67, 73, 0.92), (4, 5, 0.8), (1, 4, 0.25), (0, 3, 0.0), (0, 1, 0.0)], 
	'^', 'y', '100'
)

ref_200 = RefDataType(
	[(59, 73, 0.81), (35, 46, 0.76), (6, 11, 0.55), (3, 6, 0.5), (4, 6, 0.67)],
	[(82, 83, 0.99), (3, 12, 0.25), (12, 24, 0.5), (7, 16, 0.44), (3, 7, 0.43)],
	[(6, 8, 0.75), (55, 71, 0.77), (37, 52, 0.71), (8, 10, 0.8), (1, 1, 1.0)],
	[(90, 95, 0.95), (7, 17, 0.41), (4, 11, 0.36), (3, 7, 0.43), (3, 11, 0.27)],
	'v', 'g', '200'
)

ref_400 = RefDataType(
	[(98, 118, 0.83), (62, 79, 0.78), (22, 27, 0.81), (13, 18, 0.72), (10, 12, 0.83)],
	[(146, 147, 0.99), (24, 35, 0.69), (17, 39, 0.44), (11, 19, 0.58), (7, 13, 0.54)],
	[(17, 22, 0.77), (105, 135, 0.78), (67, 77, 0.87), (11, 14, 0.79), (5, 6, 0.83)],
	[(174, 188, 0.93), (13, 23, 0.57), (7, 14, 0.5), (8, 14, 0.57), (3, 15, 0.2)],
	'h', 'c', '400'
)

ref_800 = RefDataType(
	[(193, 236, 0.82), (107, 145, 0.74), (39, 55, 0.71), (18, 28, 0.64), (8, 13, 0.62)],
	[(271, 278, 0.97), (39, 83, 0.47), (30, 58, 0.52), (10, 24, 0.42), (13, 27, 0.48)],
	[(30, 46, 0.65), (174, 230, 0.76), (127, 162, 0.78), (21, 26, 0.81), (6, 6, 1.0)],
	[(310, 345, 0.9), (25, 53, 0.47), (12, 27, 0.44), (7, 23, 0.3), (11, 27, 0.41)],
	's', 'r', '800'
)

ref_1600 = RefDataType(
	[(325, 404, 0.8), (181, 225, 0.8), (46, 72, 0.64), (35, 53, 0.66), (19, 25, 0.76)],
	[(432, 442, 0.98), (72, 130, 0.55), (48, 95, 0.51), (29, 58, 0.5), (19, 40, 0.47)],
	[(70, 104, 0.67), (253, 328, 0.77), (134, 173, 0.77), (119, 137, 0.87), (20, 23, 0.87)],
	[(500, 548, 0.91), (46, 71, 0.65), (25, 50, 0.5), (23, 62, 0.37), (12, 48, 0.25)], 
	'o', 'b', '1600'
)


# plots
props = ['length', 'steps', 'coverage', 'cv']
f, axarr = plt.subplots(2, len(props), sharey='row', sharex='col')
# axarr[0].set_xticklabels(labels)

for ind, prop in enumerate(props):
	for ref in [ref_100, ref_200, ref_400, ref_800, ref_1600]:
		# print ref, prop, ind
		ref_prop_plot(ref, prop, ind, ranges)
		ref_counts_plot(ref, prop, ind, ranges)
	if prop == 'cv':
		axarr[0][ind].set_xlim(left=0, right=0.25)
		axarr[1][ind].set_xlim(left=0, right=0.25)
	if prop == 'steps':
		axarr[0][ind].set_xscale('log', basex=2)
		axarr[0][ind].set_xlim(1,32)
		axarr[1][ind].set_xscale('log', basex=2)
		axarr[1][ind].set_xlim(1,32)
	if prop == 'coverage':
		axarr[0][ind].set_xscale('log', basex=10)
		axarr[1][ind].set_xscale('log', basex=10)
	axarr[0][ind].set_title(prop.upper())

	# x_vals = get_x_mids(ranges[ind])
	# row = ref_1600.get_prop(prop)
	# instances = np.array([a[1] for a in row])
	# instances = np.true_divide(instances, sum(instances))
	# axarr[1][ind].scatter(x_vals, instances, 
	# 	c=ref_1600.color, marker=ref_1600.marker, s=100)
	# axarr[1][ind].plot(x_vals, instances, c=ref_1600.color)

# legend - put it on the rightmost

axarr[0][ind].legend()
plt.show()
