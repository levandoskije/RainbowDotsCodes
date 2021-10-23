#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 16:38:26 2021

@author: levandoskije
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tkinter import *


def conc_enzimas(E1, i=1):
    
    E2 = E3 = E4 = E5 = E6 = E7 = E8 = E9 = E10 = E1
    n = 2
    if i == 2:
        E2 = n*E1
    elif i == 3:
        E3 = n*E1
    elif i == 4:
        E4 = n*E1
    elif i == 5:
        E5 = n*E1
    elif i == 6:
        E6 = n*E1
    elif i == 7:
        E7 = n*E1
    elif i == 8:
        E8 = n*E1
    elif i == 9:
        E9 = n*E1
    elif i == 10:
        E10 = n*E1
    elif i == 0:
        
        E10 = n*E1
        E3 =  n*E1
        E8 =  n*E1
        
        
        
    
    return E2, E3, E4, E5, E6, E7, E8, E9, E10
    


def ODEs(y, t, I, i=1):
    """
    
    """
    S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, P = y
    # ----------------------Parameters---------------------------
        
    S0 = 1
    
    k1 = 10 #s-1
    k2 = 200 #s-1
    k_2 = 3000 #s-1 
    k3 = 6 #s-1
    k4 = 10 #s-1
    k5 = 20 #s-1
    k6 = 4 #s-1
    k7 = 1 #s-1
    
    k8 = 8.27 #s-1
    k9 = 1.4 #s-1
    k10 = 0.079 #s-1
    k_on = 0.0043 
    
    Km1 = 0.65 #mM
    Ki1 = 0.014 #mM
    Km2 = 1 #mM
    Km_2 = 0.01 #mM
    Km3a = 0.015 #mM
    Km3b = 0.003 #mM
    Ki3 = 0.01 #mM
    Km4 = 0.02 #mM
    Km5 = 0.1 #mM
    Km6 = 0.9 #mM
    Km7 = 0.2 #mM
    
    Km8 = 0.0225 #mM
    Ki8 = 0.299 #mM
    Km9 = 0.39 #mM
    Ki9 = 0.573 #mM
    Km10 = 0.00017 #mM
    
    A = 10 
    e = 0.01
    
    E1 = 1
    E2, E3, E4, E5, E6, E7, E8, E9, E10 =  conc_enzimas(E1, i)
    
    kl2 = e*k2*E2/(k1*E1) #0.2
    kl_2 = e*k_2*E2/(k1*E1) #3
    kl3 = k3*E3/(k1*E1) #0.6
    kl4 = k4*E4/(k1*E1) #1
    kl5 = k5*E5/(k1*E1) #2
    kl6 = k6*E6/(k1*E1) #0.4
    kl7 = k7*E7/(k1*E1) #0.1
    
    
    kl8 = k8*E8/(k1*E1) #0.827
    kl9 = k9*E9/(k1*E1) #0.14
    kl10 = k10*E10/(e*k1*E1) #0.79
    kl_on =  k_on/e
    
    Klm1 = Km1/S0 # 0.65
    Kli1 = Ki1/(e*S0) #1.4
    Klm2 = Km2/S0 #1
    Klm_2 = Km_2/(e*S0) #1
    Klm3a = Km3a/(e*S0) #1.5
    Klm3b = Km3b/(e*S0) #0.3
    Kli3 = Ki3/(e*S0) #1
    Klm4 =  Km4/(e*S0) #2
    Klm5 =  Km5/(S0) # 0.1
    Klm6 =  Km6/(S0) #0.9
    Klm7 =  Km7/(S0) #0.2
    
    Klm8 =  Km8/(e*S0) #2.25
    Kli8 = Ki8/(S0) #0.299
    Klm9 =  Km9/(S0) #0.39
    Kli9 = Ki9/(S0) #0.573
    Klm10 =  Km10/(e*e*S0) #1.7
    
    Al = A*S0/(k1*E1) #1
    
    
    
    #----------------------------------------
    dS1dt = (I*kl_on - 1)*(e*Kli1*S1)/(e*Kli1*(S1 + Klm1) + S1*S2)
    
    dS2dt = (e*Kli1*S1)/(e*Kli1*(S1 + Klm1) + S1*S2) \
        - kl2*S2/(e*(S2 + Klm2)) + kl_2*S3/(e*(S3 + e*Klm_2)) - Al*S2
    
    dS3dt = kl2*S2/(e*(S2 + Klm2)) - kl_2*S3/(e*(S3 + e*Klm_2)) \
        - (kl3*Kli3*S2*S3)/(Kli3*S2*S3 + Klm3a*S3*(S3 + e*Kli3) + e*Kli3*Klm3b*S2)
    
    dS4dt = (kl3*Kli3*S2*S3)/(Kli3*S2*S3 + Klm3a*S3*(S3 * e*Kli3) + e*Kli3*Klm3b*S2) \
        - kl4*S4/(S4 + e*Klm4)
    
    dS5dt = kl4*S4/(S4 + e*Klm4) - kl5*S5/(S5 + Klm5)
    
    dS6dt = kl5*S5/(S5 + Klm5) - kl6*S6/(S6 + Klm6)
    
    dS7dt = kl6*S6/(S6 + Klm6) - kl7*S7/(S7 + Klm7)
    
    #IPP
    dS8dt = (kl7*S7)/(S7 + e*Klm7) \
        - (e*kl8*Kli8*S8)/(e*Kli8*(S8 + Klm8) + S8*S9) 
    #DMAP
    dS9dt = (e*kl8*Kli8*S8)/(e*Kli8*(S8 + Klm8) + S8*S9) \
         	-(e*Kli9*S9)/(e*Kli9*(S9 + Klm9))
                       
    #GPP
    dS10dt = (e*Kli9*S9)/(e*Kli9*(S9 + Klm9))\
    		- kl10*S10/(S10 + Klm10)
    
    dPdt = kl10*S10/(S10 + Klm10)
    # --------------------------------------
   
    dydt = [dS1dt, dS2dt, dS3dt, dS4dt, dS5dt, dS6dt, dS7dt, dS8dt,  dS9dt, dS10dt, dPdt]
    
    return dydt
#------------------------------------------------

def optimized_system():
    
    print('Optimized system!')
    
    s1 = 'Pyruvate'
    s2 = 'Acetyl-CoA'
    s3 = 'Acetoacetyl-CoA'
    s4 = 'HMG-CoA'
    s5 = 'Mevalonate'
    s6 = 'Mevalonate phosphate'
    s7 = 'Mevalonate diphosphate'
    s8 = 'IPP'
    s9 = 'DMAPP'
    s10 = 'GPP'
    p = '1,8-Cineole'
    
    
    y0 = np.array([1,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0])+1e-14
    t = np.linspace(0, 100, 1000) # s
    
    lista = [0,]#[0, 0.25, 0.5, 0.75, 1]
    
    lista = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    plt.figure(dpi=150, figsize=(7, 5))
    I = 0
    lista = [1, 0]
    for i in lista:
    
        sol = odeint(ODEs, y0, t, args=(I, i))
        S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, P= sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], \
            sol[:, 4], sol[:, 5], sol[:, 6], sol[:, 7], sol[:, 8], sol[:, 9], sol[:, 10]
        
        if i == 1:
            Pref = P
            s8ref = S8
            plt.plot(t, P, label=f'{p:20}         Ref.')
            plt.plot(t, S8, label=f'{s8:20}       Ref.')
        elif i == 0:
            Popt = P
            s8opt = S8
            plt.plot(t, P, '-.', label=f'{p:20}         Optimized')
            plt.plot(t, S8,  '--', label=f'{s8:20}         Optimized')
            
    plt.title(f"IPP and 1,8 - Cineole\n[cs] = [idi] = [hmgcs] = 2*[pdh]")
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(bbox_to_anchor=(1.1, 1.05),
          ncol=1, fancybox=True, shadow=True)
    plt.xlim(0, max(t))
    plt.ylim(0,)
    plt.savefig('Production_optimized.png', type='png')
    plt.show()

#--------------------------------------------------------

def def_enzymes():
    
    print('Analysis of Enzymes influence!')
    
    s1 = 'Pyruvate'
    s2 = 'Acetyl-CoA'
    s3 = 'Acetoacetyl-CoA'
    s4 = 'HMG-CoA'
    s5 = 'Mevalonate'
    s6 = 'Mevalonate phosphate'
    s7 = 'Mevalonate diphosphate'
    s8 = 'IPP'
    s9 = 'DMAPP'
    s10 = 'GPP'
    p = '1,8-Cineole'
    
    
    y0 = np.array([1,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0])+1e-14
    t = np.linspace(0, 100, 1000) 
    
    
    lista = [2, 3, 4, 5, 6, 7, 8, 9, 10, 1]
    
    
    plt.figure(dpi=150, figsize=(7, 5))
    I = 0
    for i in lista:
    
        sol = odeint(ODEs, y0, t, args=(I, i))
        S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, P= sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], \
            sol[:, 4], sol[:, 5], sol[:, 6], sol[:, 7], sol[:, 8], sol[:, 9], sol[:, 10]
        
        if i == 1:
            plt.plot(t, P, label=f'{p:20} Ref')
        else:
            plt.plot(t, P, label=f'{p:20} E{i} = 2*E1')
    
    plt.title(f"Production of 1,8 - Cineole")
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(bbox_to_anchor=(1.1, 1.05),
          ncol=1, fancybox=True, shadow=True)
    plt.xlim(0, max(t))
    plt.ylim(0,)
    plt.savefig('Production_cineole.png', type='png')
    plt.show()
    
    
    
    
    plt.figure(dpi=150, figsize=(7, 5))
    I = 0
    for i in lista:
    
        sol = odeint(ODEs, y0, t, args=(I, i))
        S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, P= sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], \
            sol[:, 4], sol[:, 5], sol[:, 6], sol[:, 7], sol[:, 8], sol[:, 9], sol[:, 10]
        
        if i == 1:
            plt.plot(t, S8, label=f'{s8:5} Ref')
        else:
            plt.plot(t, S8, label=f'{s8:5} E{i} = 2*E1')
    
    plt.title(f"Accumulation of IPP")
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(bbox_to_anchor=(1.1, 1.05),
          ncol=1, fancybox=True, shadow=True)
    plt.xlim(0, max(t))
    plt.ylim(0,)
    plt.savefig('Accumulation_ipp.png', type='png')
    plt.show()


#------------------------------------------------------------
def light_boost():
    
    print('Light influence!')
    s1 = 'Pyruvate'
    s2 = 'Acetyl-CoA'
    s3 = 'Acetoacetyl-CoA'
    s4 = 'HMG-CoA'
    s5 = 'Mevalonate'
    s6 = 'Mevalonate phosphate'
    s7 = 'Mevalonate diphosphate'
    s8 = 'IPP'
    s9 = 'DMAPP'
    s10 = 'GPP'
    p = '1,8-Cineole'
    
    power = 60
    
    
    y0 = np.array([1,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0])+1e-14
    t = np.linspace(0, 100, 1000) 
    
    
    lista = [0, 0.5, 1, 1.5, 2, 2.25, 2.5]
    lista = np.linspace(0, 1, 3)
    
    plt.figure(dpi=150, figsize=(7, 5))
    i = 0
    for I in lista:
    
        sol = odeint(ODEs, y0, t, args=(I, i))
        S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, P= sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], \
            sol[:, 4], sol[:, 5], sol[:, 6], sol[:, 7], sol[:, 8], sol[:, 9], sol[:, 10]
            
        plt.plot(t, S1, '-', label=f'{s1:20} I:{I*power}' + r'$\frac{mW}{cm^{-2}}$')
            
        plt.plot(t, P, '-', label=f'{p:20} I:{I*power}' + r'$\frac{mW}{cm^{-2}}$')
            
    plt.title(f"Influence of light in production of 1,8 - Cineole")
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.legend(bbox_to_anchor=(1.1, 1.05),
          ncol=1, fancybox=True, shadow=True)
    
    plt.xlim(0, max(t))
    plt.ylim(0,)
    plt.savefig('Influence_light.png', type='png')
    plt.show()
    
    
def command_line():
    
        
    choice = 1
    
    while choice in [1, 2, 3]:
        print('-.-#'*10)
        print('Analysis of Enzymes influence: 1\nOptimized system: 2\nLight influence: 3\nExit: any key')
        
        choice = float(input('Insert:\n'))
        
        print('Calculating:')
        if choice == 1:
            def_enzymes()
        elif choice == 2:
            optimized_system()
        elif choice == 3:
            light_boost()        
        print('Done!\n\n\n\n\n')
        
        
def gui():
    master = Tk()
    master.geometry('300x130')
    
    master.configure(bg='light grey')
    
    Button(master, text="Analysis of Enzymes influence", command=def_enzymes).place(x=20, y=20)
    Button(master, text="Optimized system", command=optimized_system).place(x=20, y=50)
    Button(master, text="Light influence", command=light_boost).place(x=20, y=80)
    Button(master, text="Exit", command=master.destroy).place(x=20, y=120)
    
    mainloop()

if __name__ == "__main__":
    
    # command_line()
    gui()
    
    
    
    
    
    
    
    
    
    
    
    
    