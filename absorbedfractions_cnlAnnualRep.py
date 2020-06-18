#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 12:27:58 2019

@author: margarita
Calculating absorber fractions based on:
    A practical method for assessment of dose conversion coefficients for aquatic biota (Ulanovski&Proehl)
CSDA ranges and mean free paths can be obtained from NIST
Phi0 are tabulated in the cited paper

This is currently set up to calculate absorbed fractions from the following nuclides:
    - Tritium (H3)
(will be added as needed)
"""

import math
import numpy as np

######################### INPUTS #############################
#cm
Biota = ['Pike','Pumpkinseed','Egg']

axis1=[60, 11, 0.1]
axis2=[7, 3, 0.1]
axis3=[9, 5, 0.1]

BiotaDensity = np.array([0.001,0.001,0.001,0.001])#kg/cm3
BiotaVolume = []  #cm3
BiotaMass=[]#kg
for a in range(len(axis1)):
    v=(4/3)*(math.pi)*(axis1[a]/2)*(axis2[a]/2)*(axis3[a]/2)
    #print(v)
    m=v*BiotaDensity[a]
    BiotaVolume.append(v)
    BiotaMass.append(m)
    
    

    
##############################################################
#set variables needed later
photon_mfp = 0.
CSDArange = 0.
#values of the fitting parameters to calculate s [ electron,photon]
a =[0.783,0.677]
b =[0.549,1.830]
x0 =[0.506,4.97]
c =[0.235,0.079]
d =[0.495,0.247]
x1 =[0.940,9.90]
C =[0.0147,0.071]

nuclide = ['H3','Sr90', 'Y90']
#Phi0 for each nuclide and each biota (pike, bullhead, clam, mollusc)
Phi0H3 = [1, 1, 1] #0.00568 MeV
Phi0Sr90 = [1, 9.90E-01, 7.28E-01] # 0.0546 MeV
Phi0Y90 = [9.73E-01, 8.81E-01, 9.90E-02] # 2.28 MeV interpolated from the given energies

for i in range(len(Biota)):
    print(str(Biota[i]))
    #biota specific
    axes=[axis1[i]/2,axis2[i]/2,axis3[i]/2] #cm
        
    mass = BiotaMass[i] #kg
    print('Volume of '+ str(Biota[i]) + ': '+str(BiotaVolume[i])+'cm^3')
    print('Mass of '+ str(Biota[i]) + ': '+str(BiotaMass[i])+'kg')
        
    #radius of a sphere of the same mass as the organism’s ellipsoid
    R0=math.pow((axes[0]*axes[1]*axes[2]),(1/3))
    #print(axes[0],axes[1], axes[2])
    #print(R0)
    ################## GEOM SPECIFIC STUFF #################################
        
    #ellipsoid axes expressed through R0
    norm_axes= []
    for j in range(len(axes)):
        norm_axes.append(axes[j]/R0)
        
    #test statement for normalization
    if abs(1-(norm_axes[0]*norm_axes[1]*norm_axes[2]))>1E-15:
        print("Normalization went wrong!")
    else:
        print('Normalized correctly, proceed')
        
        #ksi and chi are the lengths of the two shorter semiaxes, expressed in terms of the length of the longest semi-axis.
    ksi = norm_axes[1]/norm_axes[0]
    chi = norm_axes[2]/norm_axes[0]
        
     #surface of the sphere
    S0=4*math.pi*R0**2
        
    #surface of ellipsoid (approximation), approximation parameter
    p=1.6075 
    #non-sphericity parameter eta = S0/S
    ###testing the eta calculation#####
    S=4*math.pi* (((axes[0]*axes[1])**p+(axes[0]*axes[2])**p+(axes[2]*axes[1])**p)/(3))**(1/p)
    eta_alt=S0/S
    ###################################
    eta = (1/(math.pow((ksi*chi),(1/3))))*(math.pow((3/(1+math.pow(ksi,(-p))+math.pow(chi,(-p)))),(1/p))) 
    print('ratio of the surface areas of the ellipsoid and the equivalent sphere: '+str(eta))
        
    for n in range(len(nuclide)):
        if nuclide[n] == 'H3':
            radtype = 'electron'
            energy = 0.00568 #MeV
            CSDArange = 0.000103 #[g/cm^2] CSDA or cm in water 
            Phi0=Phi0H3[i]
            
        if nuclide[n] == 'Sr90':
            radtype = 'electron'
            energy = 0.1958 #MeV
            CSDArange = 0.0455 #[g/cm^2] CSDA or cm in water
            Phi0=Phi0Sr90[i]
            
        if nuclide[n] == 'Y90':
            radtype = 'electron'
            energy = 0.9936 #MeV
            CSDArange =0.408 #[g/cm^2] CSDA or cm in water
            Phi0=Phi0Y90[i]
        
        
        if radtype=='electron':
            r=0
            #CSDA range approximated for the energy
            energydep_range = CSDArange
        elif radtype=='photon':
            r=1
            energydep_range = photon_mfp
        #r0=[99353.684693568,1]
        
        def sigmoid(x):
          return 1 / (1 + math.exp(-x))
        
    
        
        ################## RAD TYPE SPECIFIC STUFF############################
        #combination of the effects of the body mass and the source particle energy.
        #This parameter is defined as the radius of the equal-mass sphere scaled
        r0=R0/energydep_range
        #print(r0)
        #parameter s for the re-scaling factor
        s= C[r]+(a[r]/(1+math.pow((r0/x0[r]),b[r])))+(c[r]/(d[r]+math.pow(sigmoid(r0/x1[r]),2)))

        #Rescaling factor
        RescalingFactor=math.pow(1-math.pow(abs(1-eta),(1/s)),s)
        print('Rescaling Factor for '+str(nuclide[n])+': ' + str(RescalingFactor))
        AbsorbedFraction = Phi0*RescalingFactor
        print('Absorbed Fraction for '+str(nuclide[n])+' (Phi): '+str(AbsorbedFraction)+'\n')
    
