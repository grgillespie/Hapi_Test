# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 10:05:29 2023

@author: grgil
"""

#beep boop
#boop
import os
import numpy as np
import matplotlib.pyplot as plt
import hapi as hapi

plt.close('all')
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

# Change the current working directory
os.chdir('C:/Users/grgil/Documents/Hapi_Test')

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

def check_file_exists(file_path):
    return os.path.exists(file_path)



#%%
hapi.db_begin('linelists') #create folder to put database
# If data already exists, use db_begin('data') to load the linelists into python I think



#convert from wavelength to wavenumber
lambda1=3 #um
lambda2=5 #um
nu_start=(lambda1*1e-4)**(-1) #cm-1
nu_end=(lambda2*1e-4)**(-1)

file_path='C:/Users/grgil/Documents/Hapi_Test/linelists/CH4.data'
if check_file_exists(file_path):
    print(f"The file '{file_path}' exists.")
else:
    hapi.fetch('CH4',6,1,nu_end,nu_start) #see hapi info 


file_path='C:/Users/grgil/Documents/Hapi_Test/linelists/C2H6.data'
if check_file_exists(file_path):
    print(f"The file '{file_path}' exists.")
else:
    hapi.fetch('C2H6',27,1,nu_end,nu_start)

'''
file_path='C:/Users/grgil/Documents/Hapi_Test/linelists/C3H8.data'
if check_file_exists(file_path):
    print(f"The file '{file_path}' exists.")
else:
    hapi.fetch('C3H8',27,1,nu_end,nu_start)
 '''   
#%% Methane
mol_id=6
iso=1
molefraction=0.95
name='CH4'
pressure=0.9 #atm
temperature=294 #K
pathlength=1.1 #cm
stepsize=0.01 #wavelength step size in cm^-1


nu, coef = hapi.absorptionCoefficient_Voigt(((int(mol_id), int(iso), molefraction),),
            name, WavenumberStep=stepsize, HITRAN_units=False, GammaL='gamma_self',
            Environment={'p':pressure,'T':temperature,'l':pathlength})


#%%
mol_id=27
iso=1
molefraction=0.04
name='C2H6'
pressure=0.9 #atm
temperature=294 #K
pathlength=1.1 #cm
stepsize=0.01 #wavelength step size in cm^-1


nu1, coef1 = hapi.absorptionCoefficient_Voigt(((int(mol_id), int(iso), molefraction),),
            name, WavenumberStep=stepsize, HITRAN_units=False, GammaL='gamma_self',
            Environment={'p':pressure,'T':temperature,'l':pathlength})



#%%
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
plt.plot(nu,coef,label='CH4, X=0.95, P=0.9 atm, T=294 K')
plt.plot(nu1,coef1,label='C2H6, X=0.04, P=0.9 atm, T=294 K')
ax1.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax1.set_ylabel('Absorbance')
ax1.legend()
secax=ax1.secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength($\mu$m)')