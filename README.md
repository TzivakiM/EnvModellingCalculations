# EnvModellingCalculations 
######################## INFO ########################################
Files included are work performed during the completion of a PhD thesis.
The correspond to results presented in Chapter 3.

First part (3.1):
absorbedfractions_cnlAnnualRep.py
AquBiotaDose_cnlAnnualRep.py

Calculates the predicted dose to aquatic biota based on concentration measurements from the CNL Annual Dose Report as
described in Section 3.1
These modules include: extended sediment-water partitioning


Second Part (3.2):
absorbedfractions.py
AquBiotaDose.py

Calculates the predicted dose to aquatic biota based on concentration measurements from the CNL Site as
described in Section 3.2, based on environmental measurements
These modules include: extended sediment-water partitioning, HTO+OBT partitioning model for H3

#################### INSTRUCTIONS ###################################

Instructions on use of the files.
1. It is recommended to use the files: absorbedfractions.py, AquBiotaDose.py which contain the newest models of this work
2. Inputs needed to run absorbedfractions.py:
    - Phi_0 (absorbed fraction for organism of the same mass) can be looked up in Ulanovsky and Proehl, 2006
    - stopping power for the energy required can be downloaded from NIST
    - geometry of ellipsoid representing the organism 
3. Run absorbedfractions.py
4. Use results of Phi (absorbed fraction for each nuclide and organisms as input for AquBiotaDose.py
5. All parameters in this file can be changed based on site specific measurements
6. Run AquBiotaDose.py
7. Output:
    - Calculated Sediment concentration based on water-sediment partitioning
    - if Tritium is used, concentration of Tritium based on OBT+HTO partitioning and BAF partitioning
    - Internal, water and sediment dose based on measured values and the 2 partitioning methods

#################### SOURCES ###################################

The stopping power of electron used for both files calculating absorbed fractions is from NIST:
https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
absorbed fractions Phi_0 are from: A. Ulanovsky, G.Proehl. 
A practical method for assessment of dose conversion coefficientsfor aquatic biota;
Radiat Environ Biophys (2006) 45:203-214
