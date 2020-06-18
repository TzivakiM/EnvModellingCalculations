#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 13:07:31 2019

@author: margarita
"""
#Perch Lake for Pumpkinseed and Eggs and Duke Stream for Pike
import math
from decimal import Decimal as dec
#Calculatioon for large fish (pike) and tritium contamination CNL site study
BiotaList = ['Pike', 'Pumpkinseed', 'Egg']

nuclideList = ['H3', 'Sr90', 'Y90']

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

#nuclide specific stuff (tritium, strontium, ytrium)
DistrCoeffList = [0.,621.,0.]
E_nuList = [0.00568,0.1958, 0.9936] #MeV
Y_nuList = [1.,1., 1.]#yield

#exposure specific stuff (pike bullhead,clam,mollusc)
WaterConcListH3 = [13E3,2.7E3,2.7E3]#Bq/L #2903
WaterConcListSr90 = [0.15,55,55]#Bq/L 
WaterConcListY90 = [0.15,55,55]#Bq/L 
#SedimentConcListH3 = [3047.,3047.,259.,0] #Bq/kg
#SedimentConcListSr90 = [0,0,0,3225] #Bq/kg
#SedimentConcListY90 = [0,0,0,3225] #Bq/kg
#BiotaConcListH3 = [4029,3946,4535.,0] # Bq/kg fresh weight
#BiotaConcListSr90 = [0,0,0,37190] # Bq/kg fresh weight
#BiotaConcListY90 = [0,0,0,37190] # Bq/kg fresh weight

porosityList = [0.6,0.6,0.6]
densityofsolidsList = [1.,1.,1.]
densityofwaterList=[1.,1.,1.]
#Biota specific (pike, bullhead, clam,mollusc)
OccupancyFactorList = [0.3,0.5,1.]
BioAccFactorListH3 = [1.,1.,1.] #L/kg
#fish, benthic, bivalve, gastropods from TRS479 (=4.9E+2 for srnot used),used site value from ecometrix
BioAccFactorListSr90 = [945,430,430] # L/kg
#fish, fish, bivalve, bivale from TRS479
BioAccFactorListY90 = [0,0,0] #L/kg


WatercontentList = [0.6,0.6,0.6]
WaterEquivalentFactorList = [0.65,0.65,0.65]#this is also the recommended value
PartitionFactorList = [0.66,0.66,0.66] #default value
PhiListH3 = [0.9995276286535315,0.9999910197235308,1] #absorbed fraction
PhiListSr90 = [0.9986867911677773,0.9898924378236275,0.728] #absorbed fraction
PhiListY90 = [0.9674166919430548,0.8795563884560679,0.099] #absorbed fraction

########################################################
#set up dictionaries for all factors
#BiotaConcDict={'H3': BiotaConcListH3,
#                  'Sr90': BiotaConcListSr90,
#                  'Y90': BiotaConcListY90
#                    }
WaterConcDict={'H3': WaterConcListH3,
                  'Sr90': WaterConcListSr90,
                  'Y90': WaterConcListY90
                    }
#SedimentConcDict={'H3': SedimentConcListH3,
#                  'Sr90': SedimentConcListSr90,
#                  'Y90': SedimentConcListY90
#                    }

BioAccFactorDict={'H3': BioAccFactorListH3,
                  'Sr90': BioAccFactorListSr90,
                  'Y90': BioAccFactorListY90
                    }
PhiDict = {'H3': PhiListH3,
           'Sr90': PhiListSr90,
           'Y90': PhiListY90
          }
#quick and dirty useful variable later. This is horrible and should be fixed at some point
Sr90SedConc=[]
Sr90BiotaConc=[]
#############################################################

for i in range(len(nuclideList)):
    print("#######################################################")
    print("Calculations for " +str(nuclideList[i])+'\n')
    DistrCoeff = DistrCoeffList[i]
    E_nu = E_nuList[i] #MeV
    Y_nu = Y_nuList[i]#yield
    #print(str(nuclideList[i]))
    for n in range(len(BiotaList)):

            print(str(BiotaList[n])+'\n')
            Biota=BiotaList[n]     
    
            #exposure specific stuff
            WaterConc = WaterConcDict[str(nuclideList[i])][n]#Bq/L
            #SedimentConc = SedimentConcDict[str(nuclideList[i])][n] #Bq/kg
            #BiotaConc = BiotaConcDict[str(nuclideList[i])][n] # Bq/kg fresh weight
            porosity = porosityList[n]
            densityofsolids = densityofsolidsList[n]
            densityofwater=densityofwaterList[n]
            OccupancyFactor = OccupancyFactorList[n]
            BioAccFactor = BioAccFactorDict[str(nuclideList[i])][n] #L/kg
            Watercontent = WatercontentList[n]
            WaterEquivalentFactor = WaterEquivalentFactorList[n]#this is also the recommended value
            PartitionFactor = PartitionFactorList[n]
            Phi = PhiDict[str(nuclideList[i])][n]
            
            ##########################################################################################
            print('\n')
            print("Calculated Values from water concentration using regular partitioning for biota:")
            if str(nuclideList[i])=='H3':
                Calc_Biotaconc = BiotaConcentration(WaterConc,BioAccFactor)
            elif str(nuclideList[i])=='Sr90':
                Calc_Biotaconc = BiotaConcentration(WaterConc,BioAccFactor)
                #print(Calc_Biotaconc)
                Sr90BiotaConc.append(Calc_Biotaconc)#this ONLY works if Sr comes before Y-90!!
                #print(Sr90BiotaConc)
            elif str(nuclideList[i])=='Y90':
                Calc_Biotaconc=Sr90BiotaConc[n]
                #print(Calc_Biotaconc)
            else:
                print('I do not know this nuclide')
        
            print("Concentration in " +str(Biota)+" : "+str('%.2E' % dec(Calc_Biotaconc))+ " [Bq/kg]")
            
            if str(nuclideList[i])=='H3':
                Calc_SedConc = SedimentConcentration(WaterConc, porosity,DistrCoeff,densityofsolids,densityofwater)
            elif str(nuclideList[i])=='Sr90':
                Calc_SedConc = SedimentConcentration(WaterConc, porosity,DistrCoeff,densityofsolids,densityofwater)
                Sr90SedConc.append(Calc_SedConc)
                #print(Sr90SedConc)
            elif str(nuclideList[i])=='Y90':
                Calc_SedConc=Sr90SedConc[n]#this ONLY works if Sr comes before Y-90!!
                print(Calc_SedConc)
            else:
                print('I do not know this nuclide')
            
            print("Concentration in Sediment "+str('%.2E' % dec(Calc_SedConc))+ " [Bq/kg]")
        
            InternalDose_calc, WaterDose_calc, SedimentDose_calc, ExternalDose_calc, TotalDose_calc= DoseCalculation(WaterConc,Calc_SedConc,Calc_Biotaconc,OccupancyFactor)
        
            print("Internal Dose for " +str(Biota)+" : "+ str('%.2E' % dec(InternalDose_calc))+ " [microGy/h]")
        
            print("Dose from Water for " +str(Biota)+" : " + str('%.2E' % dec(WaterDose_calc))+ " [microGy/h]")
        
            print("Dose from Sediment for "  +str(Biota)+" : "+ str('%.2E' % dec(SedimentDose_calc))+ " [microGy/h]")
        
            print("External Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(ExternalDose_calc))+ " [microGy/h]")
        
            print("Total Dose for "  +str(Biota)+" : "+ str('%.2E' % dec(TotalDose_calc))+ " [microGy/h]")
            
            ##########################################################################################
            print('\n')
            if str(nuclideList[i])=='H3':
                
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
