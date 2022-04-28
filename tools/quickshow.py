import numpy as np
import matplotlib.pyplot as plt

# Quick show allows you to have a visual on the model atop of data
# The read file is the one created by getmodel. You therefore need to execute getmodel in order
# to get the proper content to use quickshow()
def extract(output_modelfile='output_model.ascii'):
	f=open(output_modelfile, 'r')
	data=f.read()
	f.close()

	txt=data.split('\n')
	nmodels=len(txt[0].split()) - 2
	#
	x=np.zeros(len(txt))
	y=np.zeros(len(txt))
	m=np.zeros((nmodels, len(txt)))
	i=0
	for t in txt:
		line=t.split()
		if len(line)>0:
			x[i]=line[0]
			y[i]=line[1]
			for j in range(nmodels):
				m[j,i]=line[2+j]
			i=i+1
	x=x[0:i-1]
	y=y[0:i-1]
	m=m[:,0:i-1]
	
	return x,y, m

def quickshow(x,y,m, c=['red', 'orange', 'blue', 'cyan', 'purple']):
	nmodels=len(m[:,0])
	plt.plot(x, y, color='gray')
	for j in range(nmodels):
		plt.plot(x,m[j,:], color=c[j])
	plt.show()

output_modelfile='../build/output_model.ascii'
x,y,m=extract(output_modelfile)
quickshow(x,y,m)