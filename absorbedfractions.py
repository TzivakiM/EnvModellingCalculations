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

######################### INPUTS #############################
#cm
axis1=30
axis2=6
axis3=6
    
BiotaMass = 0.1 #kg
    
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

nuclide = ['H3','tritium']

for n in range(len(nuclide)):
    
    if nuclide[n] == 'H3':
        radtype = 'electron'
        energy = 0.00568 #MeV
        CSDArange = 0.00010326598064 #[g/cm^2] CSDA range in water
        Phi0=1
        
    #placeholder for other nuclides
    if nuclide[n] == 'tritium':
        radtype = 'electron'
        energy = 0.00568 #MeV
        CSDArange = 0.00010326598064 #[g/cm^2] CSDA range in water
        Phi0=1
    
    if radtype=='electron':
        r=0
        #CSDA range approximated for the energy
        energydep_range = CSDArange
    elif radtype=='photon':
        r=1
        energydep_range = photon_mfp
    r0=[99353.684693568,1]
    
    def sigmoid(x):
      return 1 / (1 + math.exp(-x))
    
    
    #biota specific
    axes=[axis1,axis1,axis3] #cm
    
    mass = BiotaMass #kg
    
    #radius of a sphere of the same mass as the organism’s ellipsoid
    R0=math.pow((axes[0]*axes[1]*axes[2]),(1/3))
    
    ################## RAD TYPE SPECIFIC STUFF############################
    #combination of the effects of the body mass and the source particle energy.
    #This parameter is defined as the radius of the equal-mass sphere scaled
    r0=R0/energydep_range
    #parameter s for the re-scaling factor
    s= C[r]+(a[r]/(1+math.pow((r0/x0[r]),b[r])))+(c[r]/(d[r]+math.pow(sigmoid(r0/x1[r]),2)))
    
    ################## GEOM SPECIFIC STUFF #################################
    #ellipsoid axes exressed through R0
    norm_axes= []
    for i in range(len(axes)):
        norm_axes.append(axes[i]/R0)
    
    #test statement for normalization
    if abs(1-(norm_axes[0]*norm_axes[1]*norm_axes[2]))>1E-15:
        print("Normalization went wrong!")
    else:
        print('Normalized correctly, proceed')
    
    #ksi and chi are the lengths of the two shorter semiaxes, expressed in terms of the length of the longest semi-axis.
    ksi = norm_axes[1]/norm_axes[0]
    chi = norm_axes[2]/norm_axes[0]
    
    #surface of the sphere
    S0=4*math.pi*R0
    
    #surface of ellipsoid (approximation), approximation parameter
    p=1.6075 
    #non-sphericity parameter eta = S0/S
    eta = (1/(math.pow((ksi*chi),(1/3))))*(math.pow((3/(1+math.pow(ksi,(-p))+math.pow(chi,(-p)))),(1/p))) 
    print('ratio of the surface areas of the ellipsoid and the equivalent sphere: '+str(eta))
    
    #Rescaling factor
    RescalingFactor=math.pow(1-math.pow(abs(1-eta),(1/s)),s)
    print('Rescaling Factor for '+str(nuclide[n])+': ' + str(RescalingFactor))
    AbsorbedFraction = Phi0*RescalingFactor
    print('Absorbed Fraction for '+str(nuclide[n])+' (Phi): '+str(AbsorbedFraction))