import numpy as np
from scipy.io import readsav

def print_s2noise(s2_file, limit_noise=None, exit_after=False):
    '''
    Take a s2 file and format the noise in a way that can be copy-pasted easily in a .model file
    s2_file: File of the s2 phase
    limit_noise: If set (in percent), fixes a minimum uncertainty for the Gaussian error. This prevents too small
            uncertainties for the noise in the .model file (this happens sometimes during the Hessian calculation)
    '''
    passed=False
    r=readsav(s2_file)
    noise_param=r['noise_param']
    noise_param_err=r['noise_param_err']
    print('# Noise parameters: A0/B0/p0, A1/B1/p1, A2/B2/p2, N0')
    cpt=0
    for i in range(3):
        print('{0:18.7f} {1:18.7f} {2:18.7f}'.format(noise_param[cpt], noise_param[cpt+1], noise_param[cpt+2]))
        #print('      ', noise_param[cpt], '  ', noise_param[cpt+1], '  ', noise_param[cpt+2])
        cpt=cpt+3
    print('{0:18.6f}'.format(noise_param[cpt]))
    #
    print('# Noise Info from output_s2 (use as priors): params / err_m / err_p')
    for i in range(0, len(noise_param)):
        if limit_noise == None:
            print('{0:18.7f} {1:18.7f} {2:18.7f}'.format(noise_param[i], noise_param_err[i,0], noise_param_err[i,1]))
            #print('      ', noise_param[i], '  ', noise_param_err[i,0], '  ', noise_param_err[i,1])
        else:
            if noise_param_err[i,0]/noise_param[i] < limit_noise/100 or noise_param_err[i,1]/noise_param[i] < limit_noise/100 or np.isfinite(noise_param_err[i,0]==False or np.isfinite(noise_param_err[i,1])==False):
                print('{0:18.7f} {1:18.7f} {2:18.7f}'.format(noise_param[i], limit_noise*noise_param[i]/100, limit_noise*noise_param[i]/100)) 
                passed=True
            else:  
                print('{0:18.7f} {1:18.7f} {2:18.7f}'.format(noise_param[i], noise_param_err[i,0], noise_param_err[i,1]))            
                #print('      ', noise_param[i], '  ', limit_noise*noise_param[i]/100, '  ', limit_noise*noise_param[i]/100)
    print(' WARNING: The White noise is smaller than 1e-7. IT IS APPROXIMATED TO 0 HERE. MANUALLY CHANGE THIS.')
    print(' WARNING: limit_noise was used. Some error were below the specified threshold. Threshold:', limit_noise, '%')
    if exit_after == True:
        exit()