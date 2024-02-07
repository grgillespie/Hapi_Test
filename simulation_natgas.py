# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 10:05:29 2023

@author: grgil
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import hapi as hapi
from tqdm import tqdm #this is a progress bar thing

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


#%% Methane
mol_id=6
iso=1
molefraction=0.95
name='CH4'
pressure=1.037 #atm
temperature=298 #K
pathlength=3.56 #cm
stepsize=0.033 #wavelength step size in cm^-1
#see spectra_single in td_support for the Diluent parameter, which is apparently better than the GammaL parameter.


nu, coef = hapi.absorptionCoefficient_Voigt(((int(mol_id), int(iso), molefraction),),
            name, WavenumberStep=stepsize, HITRAN_units=False, GammaL='gamma_self',
            Environment={'p':pressure,'T':temperature})

coef=coef*pathlength


#%% Ethane
mol_id=27
iso=1
molefraction=0.04
name='C2H6'
pressure=1.037 #atm
temperature=298 #K
pathlength=3.56 #cm
stepsize=0.033 #wavelength step size in cm^-1


nu1, coef1 = hapi.absorptionCoefficient_Voigt(((int(mol_id), int(iso), molefraction),),
            name, WavenumberStep=stepsize, HITRAN_units=False, GammaL='gamma_self',
            Environment={'p':pressure,'T':temperature})

coef1=coef1*pathlength
#%% Propane

current_directory=os.getcwd()
filepath_linelists='linelists'
new_path=os.path.join(current_directory,filepath_linelists)
filename="C3H8.txt"


with open(os.path.join(new_path,filename),'r') as f: #open the txt file
            
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
                   

T=298 #K
R=8.3145 #J/molK
Na=6.02214076e23
P=1.037*101325 #Pa 
molarfrac=0.01
L=3.56 #cm


numden=Na*(P/(R*T))*molarfrac #molec/m^3
numden_cm=numden/1e6 #molec/cm^3

alpha=xsec_individual*numden_cm*L
#%% Plotting

plt.close('all')

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

fig,ax1 = plt.subplots(2,1,tight_layout=True)
ax1[0].plot(nu,coef,label=r'CH$_4$, $X$=0.95')
ax1[0].plot(nu1,coef1,label=r'C$_2$H$_6$, $X$=0.04')
ax1[0].plot(v,alpha,label=r'C$_3$H$_8$, $X$=0.01')
ax1[0].set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax1[0].set_ylabel('Absorbance')
ax1[0].legend()
secax=ax1[0].secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength (${\mu}$m)')
ax1[0].set(xlim=(2200, 3333))
ax1[0].set_title(r'$P$=0.9 atm, $T$=294 K, $L$=1.1 cm', pad=20)

#fig,ax2 = plt.subplots(tight_layout=True)
ax1[1].plot(nu,coef,label=r'CH$_4$, $X$=0.95')
ax1[1].plot(nu1,coef1,label=r'C$_2$H$_6$, $X$=0.04')
ax1[1].plot(v,alpha,label=r'C$_3$H$_8$, $X$=0.01')
ax1[1].set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax1[1].set_ylabel('Absorbance')
#ax1[1].legend()
secax=ax1[1].secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength ($\mu$m)')
ax1[1].set(ylim=(0, 1.2))
ax1[1].set(xlim=(2850, 3075))
#ax2.set_title(r'$P$=0.9 atm, $T$=294 K, $L$=1.1 cm', pad=20)


#%% Methane by itself
'''
fig,ax5 = plt.subplots(tight_layout=True)
ax5.plot(nu,coef,label=r'CH$_4$, $X$=0.95')
ax5.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax5.set_ylabel('Absorbance')
ax5.legend()
secax=ax5.secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength (${\mu}$m)')
ax5.set(xlim=(2200, 3333))
#ax1[0].set_title(r'$P$=0.9 atm, $T$=294 K, $L$=1.1 cm', pad=20)
ax5.set(xlim=(2600, 3250))
ax5.set(ylim=(0, 10))

'''
#%% Sum the absorbance and plot transmission

#Pad the ethane data with zeros
numleft=round((nu1[0]-nu[0])/stepsize)
numright=round((nu[-1]-nu1[-1])/stepsize)

padded_coef1=np.pad(coef1, (numleft, numright), mode='constant', constant_values=0)

common_nu= np.linspace(min(min(nu), min(nu1)), max(max(nu), max(nu1)), num=len(nu))

methane_interp=np.interp(common_nu, nu, coef)
ethane_interp=np.interp(common_nu, nu1, coef1)
propane_interp=np.interp(common_nu, v, alpha)

sum_alpha=methane_interp+ethane_interp+propane_interp

trans=np.exp(-sum_alpha)

fig, ax2=plt.subplots(tight_layout=True)
ax2.plot(common_nu, trans)
ax2.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax2.set_ylabel('Transmission')
secax=ax2.secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength (${\mu}$m)')
ax2.set(xlim=(2631, 3236))

fig, ax3=plt.subplots(tight_layout=True)
ax3.plot(common_nu, sum_alpha)
ax3.set_xlabel(r'Wavenumber (cm$^{-1}$)')
ax3.set_ylabel('Absorbance')
secax=ax3.secondary_xaxis('top', functions=(wavelength_conversion, wavelength_conversion_rev))
secax.set_xlabel(r'Wavelength (${\mu}$m)')
ax3.set(xlim=(2631, 3236))
ax3.set(ylim=(-0.5,10))
