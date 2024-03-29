import numpy as np
from scipy.io import readsav
import os
    
def sav2data(savfile_in, datafile_out, fmin=None, fmax=None):
    '''
        Small program to convert sav file into data files
    '''
    r=readsav(savfile_in)
    freq=r['freq']
    spec_reg=r['spec_reg']
    if fmin == None:
        fmin=min(r['freq'])
    if fmax == None:
        fmax=max(r['freq'])
    print('            Converting sav file : ', savfile_in)
    print('                           Into : ', datafile_out)
    print(' Range of Frequency that is kept: [ ', fmin, ' , ', fmax,' ]')
    pos=np.where(np.bitwise_and(freq>=fmin, freq<=fmax))
    x=freq[pos]
    s=spec_reg[pos]
    f=open(datafile_out, 'w')
    header='# File auto-generated by sav2data.py\n'
    header= header + '# freq_range=[ {} , {} ]\n'.format(fmin, fmax)
    header= header + '! frequency          power\n'
    header=header +  '* (microHz)          (ppm^2/microHz)\n'
    f.write(header)
    for i in range(len(x)):
        chain='{0:20.8f} {1:20.8f}\n'.format(x[i], s[i])
        f.write(chain)
    f.close()
    print('      sav file was converted into a .data file : ' , datafile_out)

def do_all(dir_savfiles='', dir_out='', fmin=2330.0, fmax=3660.0):
    #dir_savfiles='/Users/obenomar/Work/tmp/Sun-data/VIRGO-SPM-SoHo/TESTS/Spec/'
    #dir_out='/Users/obenomar/Work/tmp/Sun-analysis/Pre-fit-outputs/final_setup/'

    files=os.listdir(dir_savfiles)
    for fin in files:
        if fin.endswith(".sav"):
            print('fin : ', fin)
            fout=os.path.join(dir_out, fin + '.data')
            print('      --> fout :', fout)
            sav2data(os.path.join(dir_savfiles, fin), fout, fmin=fmin, fmax=fmax)