# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 11:14:38 2023

@author: grgil
"""

import os
import numpy as np
import matplotlib.pyplot as plt


# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# Change the current working directory
os.chdir('C:/Users/grgil/Documents/Hapi_Test/linelists')

filename="C3H8.txt"


with open(filename,'r') as f: #open the txt file
            
            # header data. Reads number of characters based on documentation in Hitran
            molecule        = f.read(20).strip() #strip removes leading and trailing whitespaces
            v_min           = float(f.read(10))
            v_max           = float(f.read(10))
            nPoints         = int(f.read(7))
            T               = float(f.read(7))
            P               = float(f.read(6))
            sigma_max       = float(f.read(10))
            instr_res_raw   = f.read(5)
            comm_name       = f.read(15).strip()
            unused          = f.read(4)
            broadener       = f.read(3).strip()
            reference       = int(f.read(3))
            
            # xsec data
            v = np.linspace(v_min,v_max,nPoints)
            xsec_individual = np.zeros(nPoints)

            for i in range(nPoints): #range(nPoints) is just a sequence of numbers from 0 to nPoints-1. Used for lopps
                            # read newline character from previous line after every 10 points
                if i%10==0: #if i is visible by 0
                    new_line = f.read(1) #reads 1 character and assigns it to new_line
                data_point = f.read(10) #next 10 characters
                
                try:
                    xsec_individual[i] = data_point
                except:
                    xsec_individual[i] = 0 #if assigning data_point to xsec_individual doesn't work, make it 0
                   # error_log.append(f'{species_dir} | xsec: {xsec_index+1} | nPoints is too large')
                   
#%%

T=294
R=8.3145 #J/molK
Na=6.02214076e23
P=0.9*101325 #Pa 
molarfrac=0.01
L=1.1 #cm


numden=Na*(P/(R*T))*molarfrac #molec/m^3
numden_cm=numden/1e6 #molec/cm^3

alpha=xsec_individual*numden_cm*L

#%% Plot

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Computer Modern"
    })
plt.rcParams.update({'font.size': 13})

#create functions to make secondary x-axis on top
def wavelength_conversion(x):
    return 10000*(x**(-1))

def wavelength_conversion_rev(x):
    x=np.array(x,float)
    near_zero=np.isclose(x,0) #have to do some stuff so you don't divide by zero
    x[near_zero]=np.inf
    x[~near_zero]=1/(x[~near_zero]*1e-4)
    return x

fig,ax1 = plt.subplots(tight_layout=True)
plt.plot(v,alpha,label=r'C$_3$H$_8$, X=0.01, P=0.9 atm, T=294 K')
ax1.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax1.set_ylabel('Absorbance')
ax1.legend()
secax=ax1.secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength($\mu$m)')
ax1.set(xlim=(2200,3333))

