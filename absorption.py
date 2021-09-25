#%% Import necessary packages and functions
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.optimize 
from scipy.optimize import curve_fit

def find_nearest(array, value):
    #array is a 1D vector of wavelengths
    #value is the specific wavelength for which want the index
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def trim_data(x, data, limit1, limit2):
    #x is a 1D array of two theta or q values
    #data is an array of x-ray intensities
    #limit1 and limit2 are what you'd like to trime your data to 
    test = np.array(x)
    set1 = find_nearest(test,limit1)
    set2 = find_nearest(test,limit2)
    return x[set1:set2], data[set1:set2]

def back_abs(x, data, limit1, limit2):
    #x is a 1D array of wavelengths
    #data is an array of intensities
    #limit1 and limit2 are what you'd like to use for calculating background
    test = np.array(x)
    set1 = find_nearest(test,limit1)
    set2 = find_nearest(test,limit2)
    back = np.mean(data[set1:set2])
    data_correct=(data-back)
    return data_correct

def Wavelen2eV(wavelen):
    #wavelen: array like, an array of e&m wavelength in nm
    c_light = 2.99792458e8 * 1e9 #
    h_plank = 4.1357e-15 
    eV_tmp=[c_light/each*h_plank for each in wavelen]
    return np.array(eV_tmp)

def linear(x,m,b):
    return m*x + b

def tauc_fit(eV, taucscale, limit1, limit2):
    set1 = find_nearest(eV,limit1)
    set2 = find_nearest(eV,limit2)
    popt, pcov = curve_fit(linear, eV[set1:set2], taucscale[set1:set2])
    return -1*popt[1]/popt[0], popt


#%%     Update values for importing and analysis
path = r'/Volumes/GoogleDrive/My Drive/Wellesley_Courses/Private/PHYS_331/F_2021/Labs/Copper_Experiment/Absorption' #update for your directory
white_light = all_files[-1] #update for your data
all_files = sorted(glob.glob(path + "/*.csv"))
sample = all_files[-7]

Fit_range = [1.7, 1.6] #eV range to do linear fitting of band edge

Thickness = 150*10**-9 #Approximate thickness of sample
Transition_Power = 1/2 #1/2 for direct, 2 for indirect, 3/2 for forbidden etc.

#   Import Data
lightsource = pd.read_csv(white_light) #read in data from ocean optics file
lightsource_data = np.array([lightsource["Wavelength"], lightsource["Intensity"]])
experiment = pd.read_csv(sample) #read in data from ocean optics file
experiment_data = np.array([experiment["Wavelength"], experiment["Intensity"]])


#   Transfrom Data 
absorption = [Wavelen2eV(experiment_data[0,:]), (lightsource_data[1,:]-experiment_data[1,:])/lightsource_data[1,:]]
tauc_scale = [absorption[0], (absorption[0]*absorption[1]/Thickness)**1/Transition_Power]
plt.plot(tauc_scale[0], tauc_scale[1])

#   Fit Data

E_g, fit = tauc_fit(tauc_scale[0],tauc_scale[1],Fit_range[0], Fit_range[1])
print('Bandgap', E_g, 'eV')

#   Plot fit
fig1, ax1 = plt.subplots()
ax1.plot(tauc_scale[0],tauc_scale[1], 'k-', label ='data')
ax1.plot(tauc_scale[0],linear(tauc_scale[0], fit[0], fit[1]), 'r--', label = 'fit')
ax1.set_ylim(0,np.max(tauc_scale[1]))
ax1.set_xlabel('Energy [eV]')
ax1.set_ylabel('Tauc Scale')

# %% Effect of Thickness
plt.plot(nm_150[0],nm_150[1]/np.max(nm_150[1]))
plt.plot(nm_50[0],nm_50[1]/np.max(nm_50[1]))
plt.plot(nm_20[0],nm_20[1]/np.max(nm_20[1]))
plt.plot(nm_10[0],nm_10[1]/np.max(nm_10[1]))
# %% Effect of Temp
plt.plot(nm_150[0],nm_150[1]/np.max(nm_150[1])) #500 C
plt.plot(nm_150_air[0],nm_150_air[1]/np.max(nm_150_air[1]))# 500 C with air
plt.plot(nm_150_250[0],nm_150_250[1]/np.max(nm_150_250[1]))# 250 C  C

# %%
