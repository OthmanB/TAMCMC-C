import numpy as np

# This is to divide by a given factor the height in the data
# This to avoid very large values in the coviariance matrix
# Ideally this rescaling should ensure that height are below < 100
# Note that the noise background in the model file has to be
# updated accordingly !
def rescale(data_file, factor=1000):
	outfile=data_file + '.rescaled'
	f=open(data_file,'r')
	raw=f.read()
	f.close()

	data=raw.split('\n')
	Ncols=len(data[0].split())
	#x=np.zeros(len(data))
	#y=np.zeros((Ncols-1, len(data)))
	#i=0
	txt='# data taken from {} and divided y-axis with a factor {} \n'.format(data_file, factor)
	for d in data:
		line=d.split()
		if len(line) != 0:
			txt=txt + '{0:20.8f}'.format(float(line[0])) 
			#x[i]=line[0]
			for j in range(1,Ncols):
				txt=txt + '{0:20.8f}'.format(float(line[j])/factor)
				#y[j,i]=line[j]		
			txt=txt + '\n'
	print(txt)
	f=open(outfile, 'w')
	f.write(txt)
	f.close()
	print('Data successfully rescaled by a factor {}'.format(factor))
	print('The new data file is {}'.format(outfile))
	print('Do not also forget to rescale Height values of the nosie and of the modes in the .model file')
