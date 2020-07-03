"""
Created on Thu Sep 12 17:06:22 2019
@author: Sumeru

"""
#==================================================================================
#This code estimates the transmon parameters using approximate Mathieu function at 
#low temperatures from the given probed resistances and capacitor values.
#==================================================================================
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

h = 6.626e-34                                       #Plank constant
phi0 = 3.291e-16                                    #flux quantum
kb = 1.38e-23                                       #Blotzmann constant
el = 1.6e-19                                        #charge quantum

print ("===========================================================================")
print ("To find the frequency and anharmonicity call the function 'spectra(rp, cp)'.")
print ("rp = probed resistance value in kOHm. cp = qubit capacitance in fF")
print ("===========================================================================")

#==================================================================================
#From the probed resistance, the leakage current equivalent to 40kOhm resistance
#is parallel in cancelled. Then an additional 20% increase in junction resistance
#is added emperically. This value is used to estimate the Josephson energy of the
#juntion. The function prints the low temperature resistance and returns Ej in GHz
#==================================================================================

def josephson_energy(rp):                           #gets josephson energy of junction
    rlow = (rp * 1.2) * 40/(40 - rp)                #low temperature resistance
    i0 = (300e-6)/(rlow * 1e3)                      #ambegaokar baratoff formula
    ej = phi0 * i0 * 1e-9                           #Ej in 1/Lj unit
    ejghz = (ej/h)                                  #Ej in GHz unit
    print("Low temerature resistance",rlow, "kOhm")
    return ejghz

#==================================================================================
#From the estimated josephson energy and given capacitance, Ej/Ec is calculated. 
#This value is used to estimate Transmon spectra using  Mathieu function to calculate 
#the energy of first few levels. However the Mathieu function defined in the scipy 
#library takes only integer input for exponent 'm'. Hence we need to approximate the
#Mathieu functions using the formula: (which holds for small m and large q)
#>>>>>>> mathieu_a(m+0.5,q)=0.5*(mathieu_a(m,q)+mathieu_a(m+1,q)); m = even
#>>>>>>> mathieu_a(m+0.5,q)=0.5*(mathieu_b(m,q)+mathieu_b(m+1,q)); m = odd
#==================================================================================

def energy_level(n, ratio, ec):
    if n % 2 ==0:
        energy_n = ec * (sp.mathieu_a(n, -0.5 * ratio) + sp.mathieu_a(n + 1, -0.5 * ratio))/2
    else:
        energy_n = ec * (sp.mathieu_b(n, -0.5 * ratio) + sp.mathieu_b(n + 1, -0.5 * ratio))/2
    return energy_n

def spectra(rp, cp):
    ec = (1e-9 * el**2)/(2 * h * cp * 1e-15)        #Charging energy in GHz from given C
    print ("===========================================================================")
    print ("charging energy", ec, "GHz.")
    ej = josephson_energy(rp)
    print ("Josephson Energy:", ej, "GHz.")
    ratio = ej/ec
    print ("Ej/Ec ratio:", ratio)
    f01 = energy_level(1, ratio, ec) - energy_level(0, ratio, ec)
    f12 = energy_level(2, ratio, ec) - energy_level(1, ratio, ec)
    print ("Transmon Frequency:", f01, "GHz.")
    print ("Anharmonicity:", f12-f01, "GHz.")
    print ("===========================================================================")
    

#code plots energy levels inside or outside of the cosinusidal potentional well
#All the energy values are divided by Ej 
#code can be changed to get absolute energy levels.    
def plot_energy_levels_in_well(ratio,N):##ratio=Ej/Ec , N= number of levels
    ec=1
    ratio=ratio
    ej=ec*ratio
    e_level=np.array([])
    n1=N
    for i in range(n1):
        e_level=np.append(e_level,energy_level(i,ratio,ec))
    e_level=(e_level/ej) + 1
    U=np.array([])
    for i in range(1001):
        U=np.append(U,1.0 + np.cos(2.0*np.pi*i*(0.001)))
    angle=np.array(range(1001))*(2.0/1000)
    plt.plot(angle,U,c='k')
    plt.xlabel("Phase/Pi")
    plt.ylabel("E/Ej")
    color=['r','g','b']

    print("===============================================================")
    k=0
    for i in range(len(e_level)):
        #print(i%3)
        if e_level[i]<=2.0:
            r=0.00
            temp=np.arccos(e_level[i]-1)/(np.pi)
            temp1=temp/2.0
            temp2=(np.abs(1.0-temp)+1.0)/2.0
            #print(temp1,temp2)
            d=(temp2-temp1)*0.05
            plt.axhline(y=e_level[i],xmin=temp1+d,xmax=temp2-d,c=color[i%len(color)],linestyle="-")
        else:
            #global k
            plt.axhline(y=e_level[i],xmin=0.05,xmax=0.95,c='k',linestyle=":")
            if k==0:
                print("Levels from n=",str(i)+" lie outside of poetential well")
            k+=1
    if k==0:
        print("All the levels are inside the potential well")
    print("===============================================================")