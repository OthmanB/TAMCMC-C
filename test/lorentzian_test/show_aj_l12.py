from main import eta0_fct
from main import lorentzian_cpp
from main import set_expectations
from acoefs import Snlm
#from function_rot import amplitude_ratio
import matplotlib.pyplot as plt
import numpy as np

def set_arguments_aj_model(l):
	# function that list all of the arguments that can be fed into the 'lorentzian_test' C++ function
	# and set their default values. Names are set identical to their counterpart in the C++ function (see './lorentzian_test get_all' for details)
	# Usual function for the definition for the Lorentzian components in a global fit:
	# -------------------
	Dnu_star=135.1
	# -------------------
	H_l=1    									# Maximum height
	fc_l=1000  									# Central frequency
	a1=2    									# a1 coeficient
	eta0=0 #eta0_fct(Dnu=Dnu_star) 			    # coefficient of distorsion due to the centrifugal force: DO NOT USE JOINTLY WITH a2
	a2=0.250									# a2 coeficient
	a3=0.100			    					# a3 coeficient
	a4=0.050									# a4 coeficient
	a5=0   			    					    # a5 coeficient
	a6=0                					    # a6 coeficient
	asym=0 			        				    # Asymetry of the Lorentzian normalised by width
	gamma_l=0.5              					# Width of the mode
	inclination=40       					    # Stellar inclination
	#V=amplitude_ratio(l, inclination)          # Mode relative height. Varies with l and inclination
	#handled_names=['H_l', 'fc_l', 'a1', 'eta0', 'a2', 'a3', 'a4', 'a5', 'a6', 'asym', 'gamma_l', 'inclination', 'V', 'l']
	#return handled_names, [H_l, fc_l, a1, eta0, a2, a3,a4, a5, a6, asym, gamma_l, inclination, V, l]
	#H_l,  fc_l, a1,  a2, a3, a4, a5, a6, eta0, asym, gamma_l, inclination, l
	handled_names=['H_l', 'fc_l', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'eta0', 'asym', 'gamma_l', 'inclination', 'l']
	return handled_names, [H_l, fc_l, a1, a2, a3,a4, a5, a6, eta0, asym, gamma_l, inclination, l]


def nice_show(freq, power, l, aj_inputs, xlines,ylines, xnames, ynames, text_index=None, file_out='test.jpg'):
	# Function that take the power spectrum that is created with the lorentzian_cpp() function
	# and visualise it along with the expected frequency positions as per
	# defined by the params vector
	unit=1e3
	freq=freq*np.repeat(unit, len(freq))
	xlines=xlines*np.repeat(unit, len(xlines))
	dx=(freq[2]-freq[1])
	dnu=0.4
	f0=xlines[xnames.index('fc_l')]
	propx=0.98
	propy=1.02
	fig, ax=plt.subplots(1, figsize=(12, 6))
	f=freq -np.repeat(f0,len(freq))
	ax.plot(f, power)
	ax.set_ylim([0, max(power)*1.2])
	ax.set_xlim([min(f), max(f)])
	ax.set_xlabel(r'$\nu_{nlm}-\nu_{0}$ (nHz)', fontsize=16)
	ax.set_ylabel("Power", fontsize=16)
	if text_index != None:
		ax.annotate(text_index, xy=(0.05, 0.92), xycoords=ax.transAxes, fontsize=22)

	m=-l
	for i in range(len(xlines)):
		x=xlines[i]
		posOK=np.where(np.bitwise_and(np.asarray(freq) >= (x - dx), np.asarray(freq) <= (x + dx)))
		if len(posOK[0]) > 1:
			px=power[posOK[0][0]]
		else:
			px=power[posOK[0]]
		if  xnames[i] == 'fc_l':
			col='green'
			ax.axvline(x=x-f0, color=col, linestyle='--')
		else:
			if m == 0:
				if l==1:
					fix_x=dnu
				if l==2:
					fix_x=1.5*dnu
			else:
				fix_x=0
			col='blue'	
			ax.plot([x-f0,x-f0], [0, px], color=col, linestyle='--')
			ax.annotate('m=' + str(m), xy=((x-f0)*propx-fix_x, px*propy), fontsize=12, color=col)
			m=m+1
	#
	if l==1:
		ax.text(0.025, 0.80, r"$a_1 =$ {0:0.0f} nHz".format(aj_inputs[0]), fontsize=16, transform=ax.transAxes)
		ax.text(0.025, 0.75, r"$a_2 =$ {0:0.0f} nHz".format(aj_inputs[1]), fontsize=16, transform=ax.transAxes)
		ytext=max(power)*1.1
		ax.annotate(r"$2 a_1$", xy=((xlines[1]+xlines[3])/2 -f0, ytext+0.01), fontsize=16, ha="center")
		ax.annotate("", xy=(xlines[1]-f0, ytext), xytext=(xlines[3]-f0, ytext), va="center", ha="center", arrowprops={"arrowstyle": "<|-|>", "facecolor":"black"}, fontsize=16)
		ytext=max(power)*0.9
		ax.annotate(r"$2 a_2$", xy=((xlines[2]-f0)/2, ytext+0.01), fontsize=16, ha="center")
		ax.annotate("", xy=(xlines[2]-f0, ytext), xytext=(0, ytext), va="center", ha="center", arrowprops={"arrowstyle": "<|-|>", "facecolor":"black"}, fontsize=16)	
	if l==2: 
		ax.text(0.025, 0.80, r"$a_1 =$ {0:0.0f} nHz".format(aj_inputs[0]), fontsize=16 , transform=ax.transAxes)
		ax.text(0.025, 0.75, r"$a_2 =$ {0:0.0f} nHz".format(aj_inputs[1]), fontsize=16 , transform=ax.transAxes)
		ax.text(0.025, 0.70, r"$a_3 =$ {0:0.0f} nHz".format(aj_inputs[2]), fontsize=16 , transform=ax.transAxes)
		ax.text(0.025, 0.65, r"$a_4 =$ {0:0.0f} nHz".format(aj_inputs[3]), fontsize=16 , transform=ax.transAxes)
		S=Snlm(xlines[1:], l) 
		ytext=max(power)*1.1
		ax.annotate(r"$4 (a_1+a_3)$", xy=((xlines[1]+xlines[5])/2 -f0, ytext+0.01), fontsize=16, ha="center")
		ax.annotate("", xy=(xlines[1]-f0, ytext), xytext=(xlines[5]-f0, ytext), va="center", ha="center", arrowprops={"arrowstyle": "<|-|>", "facecolor":"black"}, fontsize=16)
		ax.plot([S[1],S[1]], [ytext*1.01, max(power)*0.9*0.925], color='orange') # S22 indicator
		#
		ytext=max(power)*0.9
		ax.annotate(r"$2 (a_1- 4a_3)$", xy=((xlines[2]+xlines[4])/2 -f0, ytext+0.01), fontsize=16, ha="center")
		ax.annotate("", xy=(xlines[2]-f0, ytext), xytext=(xlines[4]-f0, ytext), va="center", ha="center", arrowprops={"arrowstyle": "<|-|>", "facecolor":"black"}, fontsize=16)
		ax.plot([S[0],S[0]], [ytext*1.01, ytext*0.925], color='orange') # S21 indicator
		# 10a4 + 3a2 = S22 - S21
		ax.annotate("", xy=(S[0], max(power)*0.9*0.925), xytext=(S[1], max(power)*0.9*0.925), va="center", ha="center", arrowprops={"arrowstyle": "<|-|>", "facecolor":"orange"}, fontsize=16)
		ax.annotate(r"$10a_4 + 3a_2$", xy=((S[0]+S[1])/2, max(power)*0.9*0.925*0.98), va="top", ha="center", fontsize=16, color='orange')		

	#ax.axhline(y=0, color='black', linestyle='--')
	fig.savefig(file_out + '.jpg', dpi=300)
	#plt.show()


def do_l(l, file_out):
	unit=1e3 # to get nanoHz
	func_name='optimum_lorentzian_calc_aj'
	# -----  Initial setup -----
	# Frequency range on which we will draw the result
	xmin=996
	xmax=1004
	# Degrees that we want to test:
	# --------------------------
	text_index=None
	if l == 1:
		text_index='(a)'
	if l == 2:
		text_index='(b)'
	if l>1:
		xmin=xmin - (l-1)*2 # For each extra l, we extend by a1=2 microHz the edge
		xmax=xmax + (l-1)*2 # For each extra l, we extend by a1=2 microHz the edge
		
	param_names,params=set_arguments_aj_model(l)
	aj_inputs=[]
	for j in range(1,6):
		aj_inputs.append(unit*params[param_names.index('a'+str(j))])
	#
	freq, power, err=lorentzian_cpp(func_name, params, xmin, xmax)
	xlines, ylines, xnames, ynames=set_expectations(params, param_names, func_name)
	nice_show(freq, power, l, aj_inputs, xlines, ylines, xnames, ynames, text_index=text_index, file_out=file_out)
