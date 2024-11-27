#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:39:26 2024

@author: lucas
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
c = 3e8
pi = np.pi
def gamma(l,R,lam):
    k=np.pi*2/lam
    return 1/(1-R)*1/(1+(1/(1-R))**2*np.sin(k*l)**2)

def gammaTotal(l,R,lam,sol):
    return 1+(gamma(l,R,lam)-1)*sol*3/(8*np.pi)



from scipy.optimize import curve_fit

# Define a Lorentzian function
def lorentzian(x, amplitude, center, width, offset):
    return offset + amplitude / (1 + ((x - center) / width)**2)
# Define a Gaussian function
def gaussian(x, amplitude, center, width):
    return amplitude * np.exp(-((x - center) ** 2) / (2 * width ** 2))

from scipy.special import wofz

def voigt(x, amplitude, center, sigma, gamma):
    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))
    return amplitude * np.real(wofz(z))
#%%
ECS_E_c = []
ACS_E_c = []
ECS_E_a = []
ACS_E_a = []
T = [78, 100, 125, 150, 175, 200, 225, 250, 275, 300]
T_str = [str(T[i]) + " K" for i in range(len(T))]

# Importo la columna de longitudes de onda y la guardo en lambda_
df = pd.read_excel("2022_May_14_YbYLF_ECS_ACS_simplified_v_1.xls",
                   sheet_name = 0)
df = df.iloc[733:3000].values
lambda_ = df[:,0].astype(np.float64)

# Importo las columnas de absorcion y emision para cada polarizacion
# sheet name es porque cada temperatura estaba en una solapa diferente del
# excel.
for i in range(10):
    
    df = pd.read_excel("2022_May_14_YbYLF_ECS_ACS_simplified_v_1.xls",
                       sheet_name = i)

    df = df.iloc[733:3000].values
    ECS_E_c.append( df[:,1] )
    ACS_E_c.append( df[:,2] )
    ECS_E_a.append( df[:,3] )
    ACS_E_a.append( df[:,4] )

def mean_(x, counts):
    return np.sum(x*counts) / np.sum(counts)

font_size = 10
#%%
def cavity_spont_emi_rate(str_,
                          gamma_plot,
                          spectra_plot,
                          transition_linewidth,
                          lines_plot,
                          save_data):
    # Esto asume siempre que es una cavidad confocal, es decir dos espejos curvos.
    # Por lo tanto si el diametro de apertura de la lente es mayor que su largo,
    # en este caso querria decir que los espejos se atravesarian y no tiene sentido
    
    # L: Largo de la cavidad
    # b: radio de la lente (aperture diameter)
    # R: Reflectividad de los espejos

    
    if str_ == 'Feld 1987':
        L = 50e-3
        b = 0.5 * 4e-3
        R = 0.9771
        lambda_0 = 556e-9 # m

    elif str_ == 'Feld 1987 for 971 nm':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 1
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.9771
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_90':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.90
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_95':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.95
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_97':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.97
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_98':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.98
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_995':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.995
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_999':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.999
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Cavity b_0':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 1.1605 * (L/2)
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Cavity b_1':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 1.1605 * (L/2) * 0.8
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Cavity b_2':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 1.1605 * (L/2) * 0.6
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Cavity b_3':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 1.1605 * (L/2) * 0.4
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Cavity b_4':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 1.1605 * (L/2) * 0.2
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
###############################################################

    elif str_ == 'Reasonable cavity R_96 L_79965':
        print('\n'+str_+'\n')
        L = (0.79965*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_79970':
        print('\n'+str_+'\n')
        L = (0.79970*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_79975':
        print('\n'+str_+'\n')
        L = (0.79975*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_79980':
        print('\n'+str_+'\n')
        L = (0.79980*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_79985':
        print('\n'+str_+'\n')
        L = (0.79985*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_79990':
        print('\n'+str_+'\n')
        L = (0.7999*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_79995':
        print('\n'+str_+'\n')
        L = (0.79995*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80000':
        print('\n'+str_+'\n')
        L = (0.8000*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80005':
        print('\n'+str_+'\n')
        L = (0.80005*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80010':
        print('\n'+str_+'\n')
        L = (0.80010*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80015':
        print('\n'+str_+'\n')
        L = (0.80015*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80020':
        print('\n'+str_+'\n')
        L = (0.80020*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80025':
        print('\n'+str_+'\n')
        L = (0.80025*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80030':
        print('\n'+str_+'\n')
        L = (0.80030*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80035':
        print('\n'+str_+'\n')
        L = (0.80035*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80040':
        print('\n'+str_+'\n')
        L = (0.80040*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80045':
        print('\n'+str_+'\n')
        L = (0.80040*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_96 L_80050':
        print('\n'+str_+'\n')
        L = (0.80040*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.96
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
###############################################################

    elif str_ == 'Reasonable cavity R_99 L_79965':
        print('\n'+str_+'\n')
        L = (0.79965*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_79970':
        print('\n'+str_+'\n')
        L = (0.79970*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_79975':
        print('\n'+str_+'\n')
        L = (0.79975*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_79980':
        print('\n'+str_+'\n')
        L = (0.79980*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_79985':
        print('\n'+str_+'\n')
        L = (0.79985*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_79990':
        print('\n'+str_+'\n')
        L = (0.7999*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_79995':
        print('\n'+str_+'\n')
        L = (0.79995*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80000':
        print('\n'+str_+'\n')
        L = (0.8000*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80005':
        print('\n'+str_+'\n')
        L = (0.80005*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80010':
        print('\n'+str_+'\n')
        L = (0.80010*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80015':
        print('\n'+str_+'\n')
        L = (0.80015*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80020':
        print('\n'+str_+'\n')
        L = (0.80020*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80025':
        print('\n'+str_+'\n')
        L = (0.80025*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80030':
        print('\n'+str_+'\n')
        L = (0.80030*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80035':
        print('\n'+str_+'\n')
        L = (0.80035*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80040':
        print('\n'+str_+'\n')
        L = (0.80040*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80045':
        print('\n'+str_+'\n')
        L = (0.80040*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    elif str_ == 'Reasonable cavity R_99 L_80050':
        print('\n'+str_+'\n')
        L = (0.80040*9.717)*1e-6
        print(' L =',round(L * 1e6,4),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m
    else:
        print('There is no set of parameter with the name: '+str_)
        return
    

    FSR = c / 2 / L # Hz
    
    nu_0 = c / lambda_0 # Hz
    
    Finesse = pi * np.sqrt(R) / 2 / (1 - R)
    
    delta_nu = FSR / Finesse # Hz
    
    Q_nu = nu_0 / delta_nu
    
    FSR_lambda = FSR * c / nu_0**2 # m
    
    solid_angle= 8*np.pi*(b/L)**2
    
    gamma_total = np.array(gammaTotal(L,R,lambda_*1e-9,solid_angle))
    inhibition = 1-min(gamma_total)
    enhancement = max(gamma_total)
    
    print(' inhibition: ',round(100*(1-min(gamma_total)),2),'%\n',
          'enhancement: ',round(max(gamma_total),2),'%')
    print(' Cavity length over mirror diameter:',L/(2*b))
    print(' Central wavelength:', round(lambda_0*1e9,2),'nm')
    print(' Free Spectral Range:',round(FSR_lambda*1e9,2),'nm')
    
    if gamma_plot == True:
        fig, ax = plt.subplots(1, 1, dpi = 300)
        ax.plot(lambda_, gamma_total)
        plt.axhline(y=1, color='r', linestyle='dashed')
        
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel(r'$\gamma$ / $\gamma_{sp}$')
        ax.set_title('Ratio of spont-emi rate into the cavity to spont-emi free space rate into the same solid angle.', fontsize = 5)
        ax.grid()

    #################### Spectra plots #################################################
    convo = []
    for i in range(len(ECS_E_c)):
        convo.append(ECS_E_c[i] * gamma_total)

    if spectra_plot == True:
        #%
        mean_regular = np.zeros(len(ECS_E_c))
        mean_cavity = np.zeros(len(ECS_E_c))
        
        ## Esto hace los plots para todas las temperaturas
        
        fig, ax = plt.subplots(1, 2, figsize=(8,4), dpi = 300)
        for i in range(len(ECS_E_c)):
            if i == 0:
                ax[0].plot(lambda_, gamma_total,color = 'red', lw = 2, alpha = .5,
                            label = r'$F_P^{res}$ = ')
                ax[0].plot(lambda_, ACS_E_c[i],
                        label = str(T[i])+' K')
                ax[1].plot(lambda_, convo[i],
                        label = str(T[i])+' K')
            else:
                ax[0].plot(lambda_, ACS_E_c[i],
                        label = str(T[i])+' K')
                ax[1].plot(lambda_, convo[i],
                        label = str(T[i])+' K')
            mean_regular[i] = mean_(lambda_, ACS_E_c[i])
            mean_cavity[i] = mean_(lambda_, convo[i])
        
        ax[0].set_ylim(-1,6)
        ax[0].set_xlabel('wavelength (nm)')
        ax[0].set_ylabel(r'ACS ($10^{-20}$ cm$^{2}$)', fontsize = font_size)
        ax[0].set_title('Regular spectra E//c', fontsize = font_size)
        ax[0].grid()
        ax[0].legend()
        ax[1].set_ylim(-.1,50)
        ax[1].set_xlabel('wavelength (nm)')
        ax[1].set_ylabel('Amplitude (AU)')
        ax[1].set_title('Spectra in the cavity direction', fontsize = font_size)
        ax[1].grid()
        ax[1].legend(loc = 'upper right')
        plt.tight_layout()

        fig, ax = plt.subplots(1, 2, figsize=(10,4), dpi = 200)

        ax[0].plot(lambda_, ECS_E_c[0],
                label = 'Free space', color = 'C0')        
        ax[0].plot(lambda_, gamma_total,
                label = 'gamma total', color = 'C2')

        
        # ax[0].set_ylim(-2.5,52)
        # ax[0].set_xlim(971,973), ax[0].set_ylim(-.1,3)
        
        ax[0].set_xlabel('wavelength (nm)')
        ax[0].set_ylabel('Amplitude (AU)')
        ax[0].set_title('T = 78 K', fontsize = font_size)
        ax[0].grid(alpha = 0.5)
        ax[0].legend(loc = 'upper left')
        
        # Create an inset axis
        # inset_axes([x0, y0, width, height])
        ax_inset = fig.add_axes([0.25, 0.5, 0.15, 0.2])
        
        # Sample data for the inset
        x_inset = lambda_
        
        # Plot data in the inset
        ax_inset.plot(x_inset, ECS_E_c[0], color='C0', lw = 1)
        ax_inset.plot(x_inset, gamma_total/100, color='C2', lw = 1)
        ax_inset.set_xlim(960, 1000)
        ax_inset.set_ylim(-.25,12)
        ax_inset.set_xlabel('')
        ax_inset.set_ylabel('')
        ax_inset.grid(alpha = 0.5)
        
        ax[1].plot(lambda_, convo[0], color = 'C3',
                label = 'Cavity modified')
        
        ax[1].set_ylim(-5,95)
        ax[1].set_xlabel('wavelength (nm)')
        ax[1].set_ylabel('Amplitude (AU)')
        ax[1].set_title('T = 78 K', fontsize = font_size)
        ax[1].grid(alpha = 0.5)
        ax[1].legend(loc = 'upper left')
        # Create an inset axis
        # inset_axes([x0, y0, width, height])
        ax_inset = fig.add_axes([0.75, 0.5, 0.15, 0.2])
        
        # Sample data for the inset
        x_inset = lambda_
        
        # Plot data in the inset
        ax_inset.plot(x_inset, convo[0], color='C3', lw = 1)
        ax_inset.plot(x_inset, ECS_E_c[0], color='C0', lw = 1)
        ax_inset.set_xlim(960, 1000)
        ax_inset.set_ylim(-.25,12)
        ax_inset.set_xlabel('')
        ax_inset.set_ylabel('')
        ax_inset.grid(alpha = 0.5)
        
        # Define starting points and directions for the arrows
        start_points = [(960, -2), (1000, -2), (982, 21.5)]
        directions = [(40, 0), (22, 23), (40, 0)]
        
        # Loop to create multiple arrows
        for (x, y), (dx, dy) in zip(start_points, directions):
            ax[0].arrow(x, y, dx, dy, color='C1', zorder = 2,
                      head_width=0, head_length=0, lw = 3)

        plt.tight_layout()

    ##################### transition linewidths plots ################################################
        
    if transition_linewidth == True:
        lambda_gamma_fit = np.linspace(lambda_[0],lambda_[-1],100000)
        gamma_total_fit = np.array(gammaTotal(L,R,lambda_gamma_fit*1e-9,solid_angle))
        
        mask_cav = (lambda_gamma_fit > 971) & (lambda_gamma_fit < 972.4)
        
        x_subset = lambda_gamma_fit[mask_cav]
        y_subset_gamma = gamma_total_fit[mask_cav]
        # Initial guess for parameters: amplitude, center, width, offset
        initial_guess = [1, 971.7, 0.1, 0.1]
        if R > 0.996:
            initial_guess = [450, 971.7, 0.1, 0.1]

        # Fit the Lorentzian model to the subset
        params_gamma, covariance_gamma = curve_fit(lorentzian,
                                                   x_subset,
                                                   y_subset_gamma,
                                                   p0=initial_guess,
                                                   bounds=(-1, [2000, 1000, 5, 1]))
        # Extract fitted parameters
        amplitude_fit_gamma, center_fit_gamma, width_fit_gamma, offset_fit_gamma = params_gamma
        
        mask = (lambda_ > 971) & (lambda_ < 975)
        x_subset = lambda_[mask]
        # Generate fitted data
        x_fit = lambda_[mask]
        x_fit = np.linspace(lambda_[mask][0], lambda_[mask][-1],500)
        y_fit_gamma = lorentzian(x_fit, *params_gamma)
        
        # Initial guess for parameters: amplitude, center, width, offset
        initial_guess = [1, 971.7, 0.1, 0.1]
        
        ECS_centers = np.zeros(len(ECS_E_c))
        ECS_widths = np.zeros(len(ECS_E_c))
        fig, ax = plt.subplots(1, 1, figsize=(8,4), dpi = 200)
        ax.axvline(x=center_fit_gamma, color='C7', linestyle='--', linewidth=1,
           label='Cavity central wavelength')
        ax.grid()
        ax.scatter(lambda_gamma_fit[mask_cav], 2 * y_subset_gamma/max(y_fit_gamma), color='blue', s = 5)
        ax.plot(x_fit, 2 * y_fit_gamma/max(y_fit_gamma), color='red', lw = 1, label = f"Cavity line | $\lambda_0$ = {center_fit_gamma:.2f} nm | $\gamma$ = {width_fit_gamma:.2f} nm" )
                    # label = f"Amplitude: {amplitude_fit_gamma:.2f} \n Center: {center_fit_gamma:.2f} nm \n Width: {width_fit_gamma:.2f} nm")
        ax.set_title(r'Lorentzian Fit to the 971.7 emission line')
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('ECS')
        ax.set_ylim(0,2.1)
        

        for i in range(len(T)-2):

            y_subset_ECS = ECS_E_c[i][mask]
            # Fit the Lorentzian model to the subset
            params_ECS, covariance_ECS = curve_fit(lorentzian,
                                                        x_subset,
                                                        y_subset_ECS,
                                                        p0=initial_guess,
                                                        bounds=(-1, [2, 1000, 10, 1]))
            # Extract fitted parameters
            amplitude_fit_ECS, center_fit_ECS, width_fit_ECS, offset_fit_ECS = params_ECS
            # amplitude_fit_ECS, center_fit_ECS, width_fit_ECS, gamma_fit_ECS = params_ECS
            
            # Generate fitted data
    
            y_fit_ECS = lorentzian(x_fit, *params_ECS)
            # y_fit_ECS = voigt(x_fit, *params_ECS)
            ECS_centers[i] = center_fit_ECS
            ECS_widths[i] = width_fit_ECS
            
        
            # Plot the original data and the Lorentzian fit

            # ax.scatter(x_subset, y_subset_ECS, label='Data', color='blue', s = 5)
            ax.plot(x_fit, y_fit_ECS, lw = 1,
                        label = f"T = {T[i]:.2f}, $\lambda_0$ = {center_fit_ECS:.1f} nm | $\gamma$ = {width_fit_ECS:.2f} nm")
                        # label = f"Amplitude: {amplitude_fit_ECS:.2f} \n Center: {center_fit_ECS:.2f} nm \n Width gaussian: {width_fit_ECS:.2f} nm \n Width lorentzian: {gamma_fit_ECS:.2f} nm")
        ax.legend(fontsize = 8)

        fig, ax = plt.subplots(1, 1, figsize=(5,4), dpi = 200)
        ax.grid(alpha = 0.2)
        ax.scatter(T[:len(T)-2], ECS_centers[:len(T)-2], color='C3', s = 10, label = 'Central wavelength', zorder = 3)
        ax.tick_params(axis='y')
        ax.errorbar(T[:len(T)-2], ECS_centers[:len(T)-2], yerr=ECS_widths[:len(T)-2]/2, alpha = 0.8, lw = 1,
                    fmt='none', ecolor='C3', capsize=5, label='Transition linewidth', zorder = 3)
        ax.set_xlim(68,260)
        ax.set_ylim(971,973.75)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('wavelength (nm)')
        
        # Banda de la cavidad
        h_0 = center_fit_gamma - width_fit_gamma/2
        h_1 = center_fit_gamma
        h_2 = center_fit_gamma + width_fit_gamma/2
        ax.axhline(y=h_1, color='C2', linestyle='-', linewidth=0.5,
                   label ='Cavity resonance')
        v_min, v_max = -0.1, 1.1
        T_fill = np.array([68, 260])
        grad1 = ax.imshow(np.linspace(0, 0.5, 256).reshape(-1, 1),
                          cmap='Greens',
                          vmin= v_min, vmax = v_max, aspect='auto',
                          extent=[T_fill.min(), T_fill.max(), h_0, h_1],
                          origin='lower')
        poly_pos = ax.fill_between(T_fill, h_0, h_1, alpha=0.1)
        grad1.set_clip_path(poly_pos.get_paths()[0], transform=ax.transData)
        poly_pos.remove()
        
        grad2 = ax.imshow(np.linspace(0, 0.5, 256).reshape(-1, 1),
                          cmap='Greens',
                          vmin= v_min, vmax = v_max, aspect='auto',
                          extent=[T_fill.min(), T_fill.max(), h_1, h_2],
                          origin='upper')
        poly_pos = ax.fill_between(T_fill, h_1, h_2, alpha=0.1)
        grad2.set_clip_path(poly_pos.get_paths()[0], transform=ax.transData)
        poly_pos.remove()
        
        ax.legend()
    
    #################### Lines plots #################################################
    index_ = np.array([199,453,606,665,
                      694,882,944,1139,
                      1154,1199,1460,1739])
    
    lambda_lines = lambda_[index_]
    ACS_E_c_lines = []
    convo_lines = []
    for i in range(len(ACS_E_c)):
        ACS_E_c_lines.append( ACS_E_c[i][index_] )
        convo_lines.append( convo[i][index_] )
    if lines_plot == True:
        fig, ax = plt.subplots(1, 2, figsize=(8.5,4), dpi = 300)
        for i in range(len(ACS_E_c_lines)):
            ax[0].scatter(lambda_lines, ACS_E_c_lines[i],
                       s = 20, label = 'T ='+str(T[i]))
            ax[1].scatter(lambda_lines, convo_lines[i],
                       s = 20, label = 'T ='+str(T[i]))
        ax[0].set_title('Free space emission')
        ax[1].set_title('Cavity')
        ax[0].grid(), ax[1].grid()
        ax[0].legend(), ax[1].legend()
    
    ################### Save data ##################################################
    if save_data == True:
        np.savez(str_+'_ecs',
                  lambda_lines = lambda_lines,
                  convo_lines_T = convo_lines,
                  param_ = [L,b,R,width_fit_gamma,center_fit_gamma],
                  param_lines = [ECS_centers,ECS_widths],
                  inhibition = inhibition,
                  enhancement = enhancement)
    
    return
#%%
cavities_ = ['Reasonable cavity R_95',
              'Reasonable cavity R_96','Reasonable cavity R_97','Reasonable cavity R_98',
              'Reasonable cavity R_99',
              'Reasonable cavity R_995','Reasonable cavity R_999',
              'Cavity b_0', 'Cavity b_1', 'Cavity b_2', 'Cavity b_3', 'Cavity b_4']

L_N = [79965,79970,79975,79980,79985,79990,79995,
       80000,
       80005,80010,80015,80020,80025,80030,80035,80040,80045,80050]

for i in range(len(L_N)):
    cavities_.append('Reasonable cavity R_96 L_'+str(L_N[i]))

for i in range(len(L_N)):
    cavities_.append('Reasonable cavity R_99 L_'+str(L_N[i]))

for cavity in cavities_:
    cavity_spont_emi_rate(cavity, gamma_plot = False,
                                  spectra_plot = False,
                                  transition_linewidth = True,
                                  lines_plot = True,
                                  save_data = True) 
#%%
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
c = 3e8

cavities_ = ['Reasonable cavity R_95',
             'Reasonable cavity R_96','Reasonable cavity R_97','Reasonable cavity R_98',
             'Reasonable cavity R_99',
             'Reasonable cavity R_995','Reasonable cavity R_999',
             'Cavity b_0', 'Cavity b_1', 'Cavity b_2', 'Cavity b_3', 'Cavity b_4']

L_N = [79965,79970,79975,79980,79985,79990,79995,
       80000,
       80005,80010,80015,80020,80025,80030,80035,80040,80045,80050]

for i in range(len(L_N)):
    cavities_.append('Reasonable cavity R_96 L_'+str(L_N[i]))

for i in range(len(L_N)):
    cavities_.append('Reasonable cavity R_99 L_'+str(L_N[i]))

for i in range(len(cavities_)):
    cavities_[i] = cavities_[i]+'_ecs.npz'
#%%
convo_lines_T_n = []
# param_ = [L,b,R,width_fit_gamma]
cavity_parameters_lw = []
cavity_parameters_b = []
cavity_parameters_L_R_96 = []
cavity_parameters_L_R_99 = []
cavity_inhibition_lw = np.zeros(len(cavities_[:7]))
cavity_enhancement_lw = np.zeros(len(cavities_[:7]))
cavity_inhibition_b = np.zeros(len(cavities_[7:12]))
cavity_enhancement_b = np.zeros(len(cavities_[7:12]))
cavity_inhibition_L_R_96 = np.zeros(len(cavities_[12:30]))
cavity_enhancement_L_R_96 = np.zeros(len(cavities_[12:30]))
cavity_inhibition_L_R_99 = np.zeros(len(cavities_[30:]))
cavity_enhancement_L_R_99 = np.zeros(len(cavities_[30:]))
j, k = 0, 0
for i in range(len(cavities_)):
    if i == 0:
        with np.load(cavities_[i], allow_pickle = True) as data:
            lambda_lines = data['lambda_lines']
            convo_lines_T_n.append( data['convo_lines_T'] )
            cavity_parameters_lw.append( data['param_'] )
            cavity_inhibition_lw[i] = data['inhibition']
            cavity_enhancement_lw[i] = data['enhancement']
    else:
        if i <= 6 :
            with np.load(cavities_[i], allow_pickle = True) as data:
                convo_lines_T_n.append( data['convo_lines_T'] )
                cavity_parameters_lw.append( data['param_'] )
                cavity_inhibition_lw[i] = data['inhibition']
                cavity_enhancement_lw[i] = data['enhancement']
        if 6 < i <= 11:
            with np.load(cavities_[i], allow_pickle = True) as data:
                convo_lines_T_n.append( data['convo_lines_T'] )
                cavity_parameters_b.append( data['param_'] )
                cavity_inhibition_b[j] = data['inhibition']
                cavity_enhancement_b[j] = data['enhancement']
                j = j + 1
        if 11 < i <= 29:
            with np.load(cavities_[i], allow_pickle = True) as data:
                convo_lines_T_n.append( data['convo_lines_T'] )
                cavity_parameters_L_R_96.append( data['param_'] )
                cavity_inhibition_L_R_96[k] = data['inhibition']
                cavity_enhancement_L_R_96[k] = data['enhancement']
        if 29 < i:
            with np.load(cavities_[i], allow_pickle = True) as data:
                convo_lines_T_n.append( data['convo_lines_T'] )
                cavity_parameters_L_R_99.append( data['param_'] )
                cavity_inhibition_L_R_99[k] = data['inhibition']
                cavity_enhancement_L_R_99[k] = data['enhancement']
                k = k + 1

cavity_parameters_lw = np.array(cavity_parameters_lw)
cavity_parameters_b = np.array(cavity_parameters_b)
cavity_parameters_L_R_96 = np.array(cavity_parameters_L_R_96)
cavity_parameters_L_R_99 = np.array(cavity_parameters_L_R_99)

ECS_centers = np.load(cavities_[0], allow_pickle = True)['param_lines'][0]
ECS_widths = np.load(cavities_[0], allow_pickle = True)['param_lines'][1]

T = [78, 100, 125, 150, 175, 200, 225, 250, 275, 300]
T_str = [str(T[i]) + " K" for i in range(len(T))]
font_size = 14

"""
Bueno, una explicación de cómo está ordenado esto.

convo_lines_T_n tiene toda la información

convo_lines_T_n[n]:
    Emision:
    si n = 0: la fabry perot con una cavidad irrealmente mala.
    si n = 2: la fabry perot con una cavidad razonable.
    si n = 4: la fabry perot con una cavidad irrealmente buena.
    Absorcion:
    si n = 1: la fabry perot con una cavidad irrealmente mala.
    si n = 3: la fabry perot con una cavidad razonable.
    si n = 5: la fabry perot con una cavidad irrealmente buena.

convo_lines_T_n[n][i]:
    Los i-es representan cada una de las temperaturas.
        Cada vector de estos tiene la amplitud de las 12 resonancias para su
    respectiva temperatura T[i]

Después, para addresear individualmente las amplitudes de la emisión para
cada longitud de onda en cada caso habría que buscar

convo_lines_T_n[n][i][j]

donde el j está ordenado por longitud de onda.

El orden quedó esta vez, al ordenarlo por longitud de onda

06 05 16 26 04 15 25 14 36 24 35 34

que en el vector del espectro entero se busca con los índices así

index_ = [199,453,606,665,
           694,882,944,1139,
           1154,1199,1460,1739]
"""
#%%
title_ = cavities_

for n in range(len(cavities_)):
    fig, ax = plt.subplots(1, 1, figsize=(6.5,5), dpi = 200)
    for i in range(len(T)):
        ax.scatter(lambda_lines, convo_lines_T_n[n][i],
                   s = 20, label = 'T ='+str(T[i]))
        ax.set_ylabel('Emission Cross Section')
    ax.set_title(title_[n])
    ax.grid()
    ax.legend()
#%%
# param_ = [L,b,R,width_fit_gamma]
def ini_and_enh_behaviour(str_):
    c_e_i = ['k', 'C3']
    if str_ == 'Mirror reflectance':
        fig, ax = plt.subplots(1, 1, figsize=(4,4), dpi = 200)
        line1 = ax.scatter(cavity_parameters_lw[:,2],
                           cavity_enhancement_lw,
                   s = 80, marker = r"$\clubsuit$", color = c_e_i[0],
                   label = 'Max enhancement')
        ax.set_xlabel('Mirror reflectance')
        ax.set_ylabel('Maximum Enhancement')
        ax.set_ylim(-20, 500)
        
        # Create the second y-axis
        ax2 = ax.twinx()
        
        # Plot the second y-axis data
        line2 = ax2.scatter(cavity_parameters_lw[:,2],
                            cavity_inhibition_lw,
                    s = 30, marker = 'D', color =  c_e_i[1],
                    label = 'Max Inhibition')
        ax2.set_ylabel('Maximum Inhibition', color= c_e_i[1])
        ax2.tick_params(axis='y', labelcolor= c_e_i[1])
        ax2.set_ylim(0.42, 0.49)
        # Combine legends from both axes
        lines = [line1, line2]
        labels = [line.get_label() for line in lines]
        ax.legend(lines, labels, loc='upper left')
    elif str_ == 'Cavity linewidth':

        fig, ax = plt.subplots(1, 1, figsize=(4,4), dpi = 200)
        line1 = ax.scatter(cavity_parameters_lw[:,3],
                           cavity_enhancement_lw,
                   s = 80, marker = r"$\clubsuit$", color =  c_e_i[0],
                   label = 'Max enhancement')
        
        ax.set_xlabel('Cavity linewidth (nm)')
        ax.set_ylabel('Maximum Enhancement')
        ax.set_ylim(-20, 500)
        
        # Create the second y-axis
        ax2 = ax.twinx()
        
        # Plot the second y-axis data
        line2 = ax2.scatter(cavity_parameters_lw[:,3],
                            cavity_inhibition_lw,
                    s = 30, marker = 'D', color =  c_e_i[1],
                    label = 'Max Inhibition')
        ax2.set_ylabel('Maximum Inhibition', color= c_e_i[1])
        ax2.tick_params(axis='y', labelcolor= c_e_i[1])
        ax2.set_ylim(0.42, 0.49)
        # ax2.set_xscale('log')
        # Combine legends from both axes
        lines = [line1, line2]
        labels = [line.get_label() for line in lines]
        ax.legend(lines, labels, loc='upper right')

    elif str_ == 'Mirror diameter':
        fig, ax = plt.subplots(1, 1, figsize=(4,4), dpi = 200)
        
        ax.scatter(cavity_parameters_b[:,1]*1e6,
                   cavity_enhancement_b,
                   s = 80, marker = r"$\clubsuit$", color =  c_e_i[0],
                   label = 'Max enhancement')
        
        ax.scatter(cavity_parameters_b[:,1]*1e6,
                   cavity_inhibition_b,
                   s = 30, marker = 'D', color =  c_e_i[1],
                   label = 'Max Inhibition')
        
        ax.set_xlabel(r'Mirror diameter ($\mu m$)')
        ax.set_ylabel('Enhancement and inhibition')
        
        ax.set_yscale('log')
        ax.legend()
    elif str_ == 'Cavity length':
        fig, ax = plt.subplots(1, 1, figsize=(4,4), dpi = 200)
        
        ax.scatter(cavity_parameters_L_R_96[:,0]*1e6,
                    cavity_enhancement_L_R_96,
                    s = 80, marker = r"$\clubsuit$", color =  c_e_i[0],
                    label = 'Max enhancement')
        ax.scatter(cavity_parameters_L_R_96[:,0]*1e6,
                    cavity_inhibition_L_R_96,
                    s = 30, marker = 'D', color =  c_e_i[1],
                    label = 'Max Inhibition')
        
        # ax.scatter(cavity_parameters_L_R_99[:,0]*1e6,
        #             cavity_enhancement_L_R_99,
        #             s = 80, marker = r"$\circ$", color =  'C0',
        #             label = 'Max enhancement')
        # ax.scatter(cavity_parameters_L_R_99[:,0]*1e6,
        #             cavity_inhibition_L_R_99,
        #             s = 30, marker = 's', color =  'C1',
        #             label = 'Max Inhibition')
        
        ax.set_xlabel(r'Cavity length ($\mu m$)')
        ax.set_ylabel('Enhancement and inhibition')

        ax.legend()
    else:
        print('There is no set of parameter with the name: '+str_)
    
    ax.grid(alpha = 0.2)
    return

#%
# ini_and_enh_behaviour('Mirror reflectance')
# ini_and_enh_behaviour('Mirror diameter')
# ini_and_enh_behaviour('Cavity linewidth')
ini_and_enh_behaviour('Cavity length')

#%%
fig, ax = plt.subplots(1, 1, figsize=(5,4), dpi = 200)
ax.grid(alpha = 0.2)
ax.scatter(T[:len(T)-2], ECS_centers[:len(T)-2], color='k', s = 10, label = r'$2 \leftrightarrow 5$ emission', zorder = 3)
ax.tick_params(axis='y')
ax.errorbar(T[:len(T)-2], ECS_centers[:len(T)-2], yerr=ECS_widths[:len(T)-2]/2, alpha = 0.8, lw = 1,
            fmt='none', ecolor='k', capsize=5, label='emission linewidth', zorder = 3)
ax.set_xlim(68,260)
ax.set_ylim(971,973.75)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('wavelength (nm)')


cmaps_ = ['Greens', 'Blues', 'Reds']
colors_ = ['C2', 'C0', 'C3']
labels_ = ['Smaller cavity resonance','Cavity on resonance',
           'Larger cavity resonance']
j = 0
for i in [0,7,17]:
    width_fit_gamma = cavity_parameters_L[:,3][i]
    center_fit_gamma = cavity_parameters_L[:,4][i]
    # Banda de la cavidad
    h_0 = center_fit_gamma - width_fit_gamma/2
    h_1 = center_fit_gamma
    h_2 = center_fit_gamma + width_fit_gamma/2
    ax.axhline(y=h_1, color=colors_[j], linestyle='-', linewidth=0.5,
               label =labels_[j])
    v_min, v_max = -0.1, 1.1
    T_fill = np.array([68, 260])
    grad1 = ax.imshow(np.linspace(0, 0.5, 256).reshape(-1, 1),
                      cmap=cmaps_[j],
                      vmin= v_min, vmax = v_max, aspect='auto',
                      extent=[T_fill.min(), T_fill.max(), h_0, h_1],
                      origin='lower')
    poly_pos = ax.fill_between(T_fill, h_0, h_1, alpha=0.1)
    grad1.set_clip_path(poly_pos.get_paths()[0], transform=ax.transData)
    poly_pos.remove()
    
    grad2 = ax.imshow(np.linspace(0, 0.5, 256).reshape(-1, 1),
                      cmap=cmaps_[j],
                      vmin= v_min, vmax = v_max, aspect='auto',
                      extent=[T_fill.min(), T_fill.max(), h_1, h_2],
                      origin='upper')
    poly_pos = ax.fill_between(T_fill, h_1, h_2, alpha=0.1)
    grad2.set_clip_path(poly_pos.get_paths()[0], transform=ax.transData)
    poly_pos.remove()
    j = j + 1

ax.legend()