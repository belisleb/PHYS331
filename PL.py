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


#%%     Update values for importing and analysis
path = r'/Volumes/GoogleDrive/My Drive/Wellesley_Courses/Private/PHYS_331/F_2021/Labs/Copper_Experiment/PL' #update for your directory
all_files = sorted(glob.glob(path + "/*.csv"))
sample = all_files[-7]
background = all_files[1]

experiment = pd.read_csv(sample) #read in data
correction = pd.read_csv(background) #read in background

plt.plot(experiment["Wavelength"], experiment["Intensity"]-correction["Intensity"] )
# %%
# TODO update so I can five it  a row of sampels to look at 
# Make it optional to have a correction