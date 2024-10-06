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
                          lines_plot,
                          save_data):
    # Esto asume siempre que es una cavidad confocal, es decir dos espejos curvos.
    # Por lo tanto si el diametro de apertura de la lente es mayor que su largo,
    # en este caso querria decir que los espejos se atravesarian y no tiene sentido
    
    # L: Largo de la cavidad
    # b: radio de la lente (aperture diameter)
    # R: Reflectividad de los espejos

    if str_ == 'Jason M Smith 2012':
        L = 4e-6
        b = 0.5 * 2*1e-6
        R = 0.994
        lambda_0 = 630e-9 # m
    elif str_ == 'Hunger 2016':
        L = 10e-6
        b = 0.5 * 11e-6
        R = 0.99
        lambda_0 = 556e-9
    elif str_ == 'Feld 1987':
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

    elif str_ == 'Intermediate cavity for 971 nm':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * (L / 12.5) * 10
        print(' b =',round(2 * b * 1e6,2),r'$\mu$ m')
        R = 0.99
        print(' R =',round(R,4))
        lambda_0 = 971e-9 # m

    elif str_ == 'Spherical cavity for 971 nm':
        print('\n'+str_+'\n')
        L = (0.8*9.717)*1e-6
        print(' L =',round(L * 1e6,2),r'$\mu$ m')
        b = 0.5 * L
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
    gamma_total += (1 - np.min(gamma_total))
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

    convo = []
    for i in range(len(ECS_E_c)):
        convo.append(ECS_E_c[i] * gamma_total)
    
    if spectra_plot == True:
        #%
        mean_regular = np.zeros(len(ECS_E_c))
        mean_cavity = np.zeros(len(ECS_E_c))
        
        fig, ax = plt.subplots(2, 1, figsize=(8.5,6), dpi = 300)
        for i in range(len(ECS_E_c)):
            if i == 0:
                ax[0].plot(lambda_, gamma_total,color = 'red', lw = 2, alpha = .5,
                           label = r'$F_P^{res}$ = ')
                ax[0].plot(lambda_, ECS_E_c[i],
                        label = str(T[i])+' K')
                ax[1].plot(lambda_, convo[i],
                        label = str(T[i])+' K')
            else:
                ax[0].plot(lambda_, ECS_E_c[i],
                        label = str(T[i])+' K')
                ax[1].plot(lambda_, convo[i],
                        label = str(T[i])+' K')
            mean_regular[i] = mean_(lambda_, ECS_E_c[i])
            mean_cavity[i] = mean_(lambda_, convo[i])
        
        
        ax[0].set_ylim(-.1,20)
        ax[0].set_xlabel('wavelength (nm)')
        ax[0].set_ylabel(r'ACS ($10^{-20}$ cm$^{2}$)', fontsize = font_size)
        ax[0].set_title('Regular spectra E//c', fontsize = font_size)
        ax[0].grid()
        ax[0].legend()
        ax[1].set_ylim(-.1,20)
        ax[1].set_xlabel('wavelength (nm)')
        ax[1].set_ylabel('Amplitude (AU)')
        ax[1].set_title('Spectra in the cavity direction', fontsize = font_size)
        ax[1].grid()
        ax[1].legend(loc = 'upper right')
        plt.tight_layout()


    index_ = np.array([199,453,606,665,
                      694,882,944,1139,
                      1154,1199,1460,1739])
    
    lambda_lines = lambda_[index_]
    ECS_E_c_lines = []
    convo_lines = []
    for i in range(len(ECS_E_c)):
        ECS_E_c_lines.append( ECS_E_c[i][index_] )
        convo_lines.append( convo[i][index_] )
    
    
    if lines_plot == True:
        fig, ax = plt.subplots(1, 2, figsize=(8.5,4), dpi = 300)
        for i in range(len(ECS_E_c_lines)):
            ax[0].scatter(lambda_lines, ECS_E_c_lines[i],
                       s = 20, label = 'T ='+str(T[i]))
            ax[1].scatter(lambda_lines, convo_lines[i],
                       s = 20, label = 'T ='+str(T[i]))
        ax[0].set_title('Free space emission')
        ax[1].set_title('Cavity')
        ax[0].grid(), ax[1].grid()
        ax[0].legend(), ax[1].legend()
    
    if save_data == True:
        np.savez(str_+"mod",
                  lambda_lines = lambda_lines,
                  convo_lines_T = convo_lines)
    
    return

cavity_spont_emi_rate('Feld 1987 for 971 nm', gamma_plot = False,
                                              spectra_plot = False,
                                              lines_plot = False,
                                              save_data = True)
cavity_spont_emi_rate('Intermediate cavity for 971 nm', gamma_plot = False,
                                              spectra_plot = False,
                                              lines_plot = False,
                                              save_data = True)
cavity_spont_emi_rate('Spherical cavity for 971 nm', gamma_plot = False,
                                              spectra_plot = False,
                                              lines_plot= False,
                                              save_data = True)
#%%

import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
c = 3e8

names = ['Feld 1987 for 971 nm.npz',
         'Feld 1987 for 971 nm_abs.npz',
         'Intermediate cavity for 971 nm.npz',
         'Intermediate cavity for 971 nm_abs.npz',
         'Spherical cavity for 971 nm.npz',
         'Spherical cavity for 971 nm_abs.npz']

convo_lines_T_n = []

for i in range(len(names)):
    if i == 0:
        with np.load(names[i], allow_pickle = True) as data:
            lambda_lines = data['lambda_lines']
            convo_lines_T_n.append( data['convo_lines_T'] )
    else:
        with np.load(names[i], allow_pickle = True) as data:
            convo_lines_T_n.append( data['convo_lines_T'] )

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

title_ = ['Cavidad mala',
          'Cavidad mala',
          'Cavidad razonable',
          'Cavidad razonable',
          'Cavidad demasiado buena',
          'Cavidad demasiado buena']

for n in range(len(names)):
    fig, ax = plt.subplots(1, 1, figsize=(6.5,5), dpi = 300)
    for i in range(len(T)):
        ax.scatter(lambda_lines, convo_lines_T_n[n][i],
                   s = 20, label = 'T ='+str(T[i]))
        
    if n % 2 == 0:
        ax.set_ylabel('Emission Cross Section')
    else:
        ax.set_ylabel('Absorption Cross Section')
    ax.set_title(title_[n])
    ax.grid()
    ax.legend()
