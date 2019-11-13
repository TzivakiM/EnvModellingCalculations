#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 13:07:31 2019

@author: margarita
"""
import math
from decimal import Decimal as dec
#Calculatioon for large fish (pike) and tritium contamination CNL site study
BiotaList = ['Pike', 'Bullhead', 'Clam']

nuclideList = ['H3', 'Tritium']

def BiotaConcentration(Cw, BAF):
    BiotaConc_calc = BAF*Cw #Bq/kg
    return BiotaConc_calc

def BiotaConH3(Cw, WC, WEQ, Rf):
    HTOConcfw = WC *Cw #Bq/kg
    OBTConcfw = (1-WC)*WEQ*Rf*Cw #Bq/kg
    Biotaconc_calc_H3 = HTOConcfw+OBTConcfw #Bq/kg
    return HTOConcfw,OBTConcfw,Biotaconc_calc_H3

def SedimentConcentration(Cw, theta, Kd, rhoS, rhoW):
    SedimentConc_calc=(theta*Cw+(1-theta)*Cw*Kd*rhoS)/(theta*rhoW+(1-theta)*rhoS) #Bq/kg
    return SedimentConc_calc

def DoseCalculation(Cw,Cs,Cb,OF):
    Dint = 5.76E-4*E_nu*Y_nu*Phi*Cb #microGy/h
    Dw=5.76E-4*E_nu*Y_nu*(1-Phi)*Cw #microGy/h
    #print(Cs,OF,Phi)
    Dsed=2.88E-4*E_nu*Y_nu*(1-Phi)*Cs*OF #microGy/h
    Dext = Dw+Dsed
    Dtotal = Dext+Dint
    '''
    print(Dint)
    print(Dsed)
    print(Dw)
    print(Dtotal)
    print(Dext)
    '''
    return Dint, Dw, Dsed, Dext, Dtotal

#nuclide specific stuff
DistrCoeffList = [0.,0.]
E_nuList = [0.00568,0.00568] #MeV
Y_nuList = [1.,1.]#yield

#exposure specific stuff
WaterConcList = [4449.,4449.,2903.]#Bq/L
SedimentConcList = [3047.,3047.,259.] #Bq/kg
BiotaConcList = [4029,3946,4535.] # Bq/kg fresh weight
porosityList = [0.7,0.7,0.35]
densityofsolidsList = [1.,1., 2.5]
densityofwaterList=[1.,1,1]
#Pike specific
OccupancyFactorList = [0.3,0.5,1]
BioAccFactorList = [1.,1.,1.] #L/kg
WatercontentList = [0.78,0.78,0.78]
WaterEquivalentFactorList = [0.65,0.65,0.65]#this is also the recommended value
PartitionFactorList = [0.84,0.77,0.65]
PhiList = [0.99992496032721,0.999982983572409,0.999999737628655] #absorbed fraction

for i in range(len(nuclideList)):
    print("#######################################################")
    print("Calculations for " +str(nuclideList[i])+'\n')
    DistrCoeff = DistrCoeffList[i]
    E_nu = E_nuList[i] #MeV
    Y_nu = Y_nuList[i]#yield
    for n in range(len(BiotaList)):
        print(str(BiotaList[n])+'\n')
        Biota=BiotaList[n]
        #exposure specific stuff
        WaterConc = WaterConcList[n]#Bq/L
        SedimentConc = SedimentConcList[n] #Bq/kg
        BiotaConc = BiotaConcList[n] # Bq/kg fresh weight
        porosity = porosityList[n]
        densityofsolids = densityofsolidsList[n]
        densityofwater=densityofwaterList[n]
        OccupancyFactor = OccupancyFactorList[n]
        BioAccFactor = BioAccFactorList[n] #L/kg
        Watercontent = WatercontentList[n]
        WaterEquivalentFactor = WaterEquivalentFactorList[n]#this is also the recommended value
        PartitionFactor = PartitionFactorList[n]
        Phi = PhiList[n]
         
        InternalDose, WaterDose, SedimentDose, ExternalDose, TotalDose = DoseCalculation(WaterConc,SedimentConc,BiotaConc,OccupancyFactor)
    
        print("Doses calculated from field values:")
        print("Concentration in " +str(Biota)+" : "+str('%.2E' % dec(BiotaConc))+ " [Bq/kg]")
        print("Concentration in Sediment "+str('%.2E' % dec(SedimentConc))+ " [Bq/kg]")
        print("Internal Dose for " +str(Biota)+" : "+ str('%.2E' % dec(InternalDose))+ " [microGy/h]")
    
        print("Dose from Water for "  +str(Biota)+" : "+ str('%.2E' % dec(WaterDose))+ " [microGy/h]")
    
        print("Dose from Sediment for " +str(Biota)+" : " + str('%.2E' % dec(SedimentDose))+ " [microGy/h]")
    
        print("External Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(ExternalDose))+ " [microGy/h]")
    
        print("Total Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(TotalDose))+ " [microGy/h]")
        ##########################################################################################
        print('\n')
        print("Calculated Values from water concentration using regular partitioning for biota:")
        Calc_Biotaconc = BiotaConcentration(WaterConc,BioAccFactor)
        print("Concentration in " +str(Biota)+" : "+str('%.2E' % dec(Calc_Biotaconc))+ " [Bq/kg]")
        Calc_SedConc = SedimentConcentration(WaterConc, porosity,DistrCoeff,densityofsolids,densityofwater)
        print("Concentration in Sediment "+str('%.2E' % dec(Calc_SedConc))+ " [Bq/kg]")
    
        InternalDose_calc, WaterDose_calc, SedimentDose_calc, ExternalDose_calc, TotalDose_calc= DoseCalculation(WaterConc,Calc_SedConc,Calc_Biotaconc,OccupancyFactor)
    
        print("Internal Dose for " +str(Biota)+" : "+ str('%.2E' % dec(InternalDose_calc))+ " [microGy/h]")
    
        print("Dose from Water for " +str(Biota)+" : " + str('%.2E' % dec(WaterDose_calc))+ " [microGy/h]")
    
        print("Dose from Sediment for "  +str(Biota)+" : "+ str('%.2E' % dec(SedimentDose_calc))+ " [microGy/h]")
    
        print("External Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(ExternalDose_calc))+ " [microGy/h]")
    
        print("Total Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(TotalDose_calc))+ " [microGy/h]")
        ##########################################################################################
        print('\n')
        print("Calculated Values from water concentration using tritium specific partitioning for biota:")
        Calc_BiotaconcentrationH3 = BiotaConH3(WaterConc,Watercontent,WaterEquivalentFactor,PartitionFactor)
        print("Concentration in "+str(Biota)+" : "+str('%.2E' % dec(Calc_BiotaconcentrationH3[2]))+ " [Bq/kg]")
        print("HTO Concentration in "+str(Biota)+" : "+str('%.2E' % dec(Calc_BiotaconcentrationH3[0]))+ " [Bq/kg]")
        print("OBT Concentration in "+str(Biota)+" : "+str('%.2E' % dec(Calc_BiotaconcentrationH3[1]))+ " [Bq/kg]")
        Calc_SedConcH3 = SedimentConcentration(WaterConc, porosity,DistrCoeff,densityofsolids,densityofwater)
        print("Concentration in Sediment "+str('%.2E' % dec(Calc_SedConc))+ " [Bq/kg]")
    
        H3InternalDose_calc, H3WaterDose_calc, H3SedimentDose_calc, H3ExternalDose_calc, H3TotalDose_calc= DoseCalculation(WaterConc,Calc_SedConc,Calc_BiotaconcentrationH3[2],OccupancyFactor)
    
        print("Internal Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(H3InternalDose_calc))+ " [microGy/h]")
    
        print("Dose from Water for "  +str(Biota)+" : " + str('%.2E' % dec(H3WaterDose_calc))+ " [microGy/h]")
    
        print("Dose from Sediment for " +str(Biota)+" : "+ str('%.2E' % dec(H3SedimentDose_calc))+ " [microGy/h]")
    
        print("External Dose for "   +str(Biota)+" : " + str('%.2E' % dec(H3ExternalDose_calc))+ " [microGy/h]")
    
        print("Total Dose for "   +str(Biota)+" : "+ str('%.2E' % dec(H3TotalDose_calc))+ " [microGy/h]")
        print('\n')
        print('\n')
