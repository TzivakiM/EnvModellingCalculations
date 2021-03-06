#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 13:07:31 2019

@author: margarita
"""
import math
from decimal import Decimal as dec
#Calculatioon for large fish (pike) and tritium contamination CNL site study
BiotaList = ['Pike', 'Bullhead', 'Clam','Mollusc']

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
DistrCoeffList = [0.,621.,0.] #L/kg
E_nuList = [0.00568,0.1958, 0.9336] #MeV
Y_nuList = [1.,1., 1.]#yield

#exposure specific stuff (pike bullhead,clam,mollusc)
WaterConcListH3 = [4449.,4449.,2903.,0]#Bq/L #2903
WaterConcListSr90 = [0,0,0,5.]#Bq/L 
WaterConcListY90 = [0,0,0,5.]#Bq/L 
SedimentConcListH3 = [3047.,3047.,259.,0] #Bq/kg
SedimentConcListSr90 = [0,0,0,3225] #Bq/kg
SedimentConcListY90 = [0,0,0,3225] #Bq/kg
BiotaConcListH3 = [4029,3946,4535.,0] # Bq/kg fresh weight
BiotaConcListSr90 = [0,0,0,37190] # Bq/kg fresh weight
BiotaConcListY90 = [0,0,0,37190] # Bq/kg fresh weight

porosityList = [0.7,0.7,0.35,0.7]
densityofsolidsList = [1.,1.,2.5,1.]
densityofwaterList=[1.,1.,1.,1.]
#Biota specific (pike, bullhead, clam,mollusc)
OccupancyFactorList = [0.3,0.5,1.,1.]
BioAccFactorListH3 = [1.,1.,1.,1.] #L/kg
#fish, benthic, bivalve, gastropods from TRS479 (=4.9E+2 for srnot used),used site value from ecometrix
BioAccFactorListSr90 = [ 8.9E+2,1.2E+3,3.8E+2,4092] # L/kg Gastropod value
#fish, fish, bivalve, bivale from TRS479
BioAccFactorListY90 = [3.1E-1,3.1E-1,2.3E+3,2.3E+3] #L/kg


WatercontentList = [0.78,0.78,0.78,0.78]
WaterEquivalentFactorList = [0.65,0.65,0.65,0.65]#this is also the recommended value
PartitionFactorList = [0.84,0.77,0.66,0.66]#trs 479
PhiListH3 = [0.9999193953193054,0.9999808257562249,0.9999996938089017,0.9999999535986788] #absorbed fraction
PhiListSr90 = [0.990627695629863,0.9898008055794054,0.9909932638891716,0.9709938794319874] #absorbed fraction
PhiListY90 = [0.9413861587808662,0.8787452862454077,0.9437947086606648,0.7506622235028259] #absorbed fraction

########################################################
#set up dictionaries for all factors
BiotaConcDict={'H3': BiotaConcListH3,
                  'Sr90': BiotaConcListSr90,
                  'Y90': BiotaConcListY90
                    }
WaterConcDict={'H3': WaterConcListH3,
                  'Sr90': WaterConcListSr90,
                  'Y90': WaterConcListY90
                    }
SedimentConcDict={'H3': SedimentConcListH3,
                  'Sr90': SedimentConcListSr90,
                  'Y90': SedimentConcListY90
                    }

BioAccFactorDict={'H3': BioAccFactorListH3,
                  'Sr90': BioAccFactorListSr90,
                  'Y90': BioAccFactorListY90
                    }
PhiDict = {'H3': PhiListH3,
           'Sr90': PhiListSr90,
           'Y90': PhiListY90
          }
#quick and dirty useful variable later. This is horrible and should be fixed at some point
Sr90SedConc=0.
Sr90BiotaConc=0.
#############################################################

for i in range(len(nuclideList)):
    print("#######################################################")
    print("Calculations for " +str(nuclideList[i])+'\n')
    DistrCoeff = DistrCoeffList[i]
    E_nu = E_nuList[i] #MeV
    Y_nu = Y_nuList[i]#yield
    #print(str(nuclideList[i]))
    for n in range(len(BiotaList)):
        if nuclideList[i]=='H3' and BiotaList[n]== 'Mollusc':
            continue
        elif nuclideList[i]=='Sr90' and (BiotaList[n]== 'Pike' or BiotaList[n]=='Bullhead' or BiotaList[n]=='Clam'):
            continue
        elif nuclideList[i]=='Y90' and (BiotaList[n]== 'Pike' or BiotaList[n]=='Bullhead' or BiotaList[n]=='Clam'):
            continue
        else:
            print(str(BiotaList[n])+'\n')
            Biota=BiotaList[n]     
    
            #exposure specific stuff
            WaterConc = WaterConcDict[str(nuclideList[i])][n]#Bq/L
            SedimentConc = SedimentConcDict[str(nuclideList[i])][n] #Bq/kg
            BiotaConc = BiotaConcDict[str(nuclideList[i])][n] # Bq/kg fresh weight
            porosity = porosityList[n]
            densityofsolids = densityofsolidsList[n]
            densityofwater=densityofwaterList[n]
            OccupancyFactor = OccupancyFactorList[n]
            BioAccFactor = BioAccFactorDict[str(nuclideList[i])][n] #L/kg
            Watercontent = WatercontentList[n]
            WaterEquivalentFactor = WaterEquivalentFactorList[n]#this is also the recommended value
            PartitionFactor = PartitionFactorList[n]
            Phi = PhiDict[str(nuclideList[i])][n]
            
            '''
            if str(nuclideList[i])=='Y90':#set Y-90 to equilibrium concentrations With Sr90
                SedimentConc= SedimentConcDict[str(nuclideList[(i-1)])][n]
                BiotaConc = BiotaConcDict[str(nuclideList[i-1])][n]
                print('YTRIUM'+str(SedimentConc)+str(BiotaConc))
            '''    


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
            if str(nuclideList[i])=='H3':
                Calc_Biotaconc = BiotaConcentration(WaterConc,BioAccFactor)
            elif str(nuclideList[i])=='Sr90':
                Calc_Biotaconc = BiotaConcentration(WaterConc,BioAccFactor)
                Sr90BiotaConc= Calc_Biotaconc#this ONLY works if Sr comes before Y-90!!
            elif str(nuclideList[i])=='Y90':
                Calc_Biotaconc=Sr90BiotaConc
                #print(Calc_SedConc)
            else:
                print('I do not know this nuclide')
        
            print("Concentration in " +str(Biota)+" : "+str('%.2E' % dec(Calc_Biotaconc))+ " [Bq/kg]")
            
            if str(nuclideList[i])=='H3':
                Calc_SedConc = SedimentConcentration(WaterConc, porosity,DistrCoeff,densityofsolids,densityofwater)
            elif str(nuclideList[i])=='Sr90':
                Calc_SedConc = SedimentConcentration(WaterConc, porosity,DistrCoeff,densityofsolids,densityofwater)
                Sr90SedConc= Calc_SedConc
                #print(Sr90SedConc)
            elif str(nuclideList[i])=='Y90':
                Calc_SedConc=Sr90SedConc#this ONLY works if Sr comes before Y-90!!
                #print(Calc_SedConc)
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
