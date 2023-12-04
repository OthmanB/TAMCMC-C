import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


def read_two_cols(filename):
    x = []
    y = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#') or line.startswith("!") or line.startswith("*"): # Skip some lines
                continue
            else:
                line = line.strip()
                if line:
                    columns = line.split()
                    if len(columns) == 2:
                        x.append(float(columns[0]))
                        y.append(float(columns[1]))
    return np.array(x), np.array(y)

def read_debug_data(file_path):
    x = []
    model = []
    eta_squared = []
    P1 = []
    P2 = []
    N0 = []
    y_in = []
    
    variables = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.startswith('#x'):
                # Skip the header line
                continue
            
            if line.startswith('!'):
                # Extract additional variables
                variable_name, variable_value = line[1:].split('=')
                variables[variable_name.strip()] = float(variable_value.strip())
                continue
            
            values = line.split()
            
            if len(values) == 7:
                x.append( float(values[0]))
                model.append( float(values[1]))
                eta_squared.append( float(values[2]))
                P1.append( float(values[3]))
                P2.append( float(values[4]))
                N0.append( float(values[5]))
                y_in.append( float(values[6]))
    return np.asarray(x), np.asarray(model), np.asarray(eta_squared), np.asarray(P1), np.asarray(P2), np.asarray(N0), np.asarray(y_in), variables


def get_ksinorm(b, c, x):
    integral, _ = integrate.quad(lambda x: 1.0 / (1.0 + (x / b) ** c), x[0], x[-1])
    ksi = b / integral
    return ksi


def eta_squared_kallinger2014(x):
    x_nyquist = np.max(x)
    eta = np.sin(0.5 * np.pi * x / x_nyquist) / (0.5 * np.pi * x / x_nyquist)
    if x[0] == 0:
        eta[0]=1 # This to avoid the Division by 0
    return eta**2

def Kallinger2014(numax, Mass, noise_params, x, y):
    """
    OLD IMPLEMENTATION WITH ONLY 2 HARVEY-LIKE PROFILES
    Using notations from Table 2 of Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
    Note that here we assume the instrumental noise to be Pinstrument(nu) = 0
    The noise_params must have parameters in this order:
        - Noise a : ka, sa, t
        - Noise b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
        - Noise b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
    Such that at the end we have: [ka, sa, t, k1, s1, c1, k2, s2, c2, N0]
    """
    N0=np.zeros(len(y)) + noise_params[-1]
    Power = y + N0
    # Compute the Leakage effect as a sinc function (Eq 1 of Kallinger+2014)
    eta_squared=eta_squared_kallinger2014(x)
    # Compute b1, b2, and a
    a = np.abs(noise_params[0] * numax ** noise_params[1] * Mass ** noise_params[2])
    b1 = np.abs(noise_params[3] * numax ** noise_params[4])
    b2 = np.abs(noise_params[6] * numax ** noise_params[7])
    c1 = np.abs(noise_params[5])
    c2 = np.abs(noise_params[8])

    # Compute the normalization constants ksi1 and ksi2
    ksi1 = get_ksinorm(b1, c1, x)
    ksi2 = get_ksinorm(b2, c2, x)

    # First SuperLorentzian
    P1 = eta_squared * ksi1 * a ** 2 / b1 * ( ( x / b1) ** c1 + 1) ** -1  # Numerator/Denominator
    Power += P1

    # Second SuperLorentzian
    P2 = eta_squared * ksi2 * a ** 2 / b2 * (( x / b2) ** c2 + 1) ** -1  # Numerator/Denominator
    Power += P2

    variables=[a,b1,b2,c1,c2,ksi1,ksi2]
    return Power, eta_squared, P1, P2, N0, variables

def Kallinger2014_V2(numax, noise_params, x, y):
    """
    NEW IMPLEMENTATION WITH 3 HARVEY: Low-Frequency (GRANULATION), Mid-Frequency and High-Frequency (under the modes)
    Using notations from Table 2 of Kallinger+2014 (https://arxiv.org/pdf/1408.0817.pdf)
    Note that here we assume the instrumental noise to be Pinstrument(nu) = 0
    The noise_params must have parameters in this order:
        - Noise a : ka, sa, t
        - Noise b1: k1, s1, ( and c1, the slope of the SuperLorentzian)
        - Noise b2: k2, s2, ( and c2, the slope of the SuperLorentzian)
    Such that at the end we have: [ka, sa, t, k1, s1, c1, k2, s2, c2, N0]
    NOTE THAT THIS IS FOR TESTS AND THAT THE GRANULATION TERM IS HARD-CODED HERE (See Below, parameters p, a0, b0, c0, ksi0)
    """
    N0=np.zeros(len(y)) + noise_params[8]
    Power = y + N0
    # Compute the Leakage effect as a sinc function (Eq 1 of Kallinger+2014)
    eta_squared=eta_squared_kallinger2014(x)
    
    # Compute b1, b2, and a
    a = np.abs(noise_params[0] * numax ** noise_params[1])
    b1 = np.abs(noise_params[2] * numax ** noise_params[3])
    b2 = np.abs(noise_params[5] * numax ** noise_params[6])
    c1 = np.abs(noise_params[4])
    c2 = np.abs(noise_params[7])

    # Compute the normalization constants ksi1 and ksi2
    ksi1 = get_ksinorm(b1, c1, x)
    ksi2 = get_ksinorm(b2, c2, x)

    # Granulation 
    p=[3335, -0.564, 836, -0.886, 2]
    a0=np.abs(p[0]* numax**p[1])
    b0=np.abs(p[2]* numax**p[3])
    c0=p[4]
    ksi0= get_ksinorm(b0, c0, x)
    P0 = eta_squared * ksi0 * a0 ** 2 / b0 * ( ( x / b0) ** c0 + 1) ** -1  # Numerator/Denominator
    Power += P0

    # First SuperLorentzian
    P1 = eta_squared * ksi1 * a ** 2 / b1 * ( ( x / b1) ** c1 + 1) ** -1  # Numerator/Denominator
    Power += P1

    # Second SuperLorentzian
    P2 = eta_squared * ksi2 * a ** 2 / b2 * (( x / b2) ** c2 + 1) ** -1  # Numerator/Denominator
    Power += P2
    variables=[a,b1,b2,c1,c2,ksi1,ksi2]
    return Power, eta_squared, P0, P1, P2, N0, variables

def smooth_data(x, ydata, delta):
    num_bins = int(delta / (x[1] - x[0]))
    smoothed_ydata = np.convolve(ydata, np.ones(num_bins)/num_bins, mode='same')
    return smoothed_ydata


def show_Kallinger2014_V2():
    numax=1500.
    noise_params=[3382, -0.609, 0.317  ,  0.970  ,  2    ,   0.948 ,  0.992,    2.  ,    100]# ,  287.578545 , 200.0000 , 22.544956]
    x=np.linspace(0, 8000, 8000)
    y=np.zeros(len(x))
    Power, eta_squared, P0, P1, P2, N0, variables=Kallinger2014_V2(numax, noise_params, x, y)
    fix, ax = plt.subplots(1)
    ax.plot(x, Power, color="black", label="P1+P2+N0")
    ax.plot(x, P0+N0, color="cyan", label="P0+N0")
    ax.plot(x, P1+N0, color="red", label="P1+N0")
    ax.plot(x, P2+N0, color="blue", label="P2+N0")
    ax.plot(x, N0, color="gray", label="N0")
    ax.legend()
    ax.set_yscale("log")
    ax.set_xscale("log")
    plt.show()
