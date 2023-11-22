# Simulation Notes
'''   '''

# Simulation Modes
ECM_Production = True
Agar_Drying = True
Decreasing_Nutrient = False
Random_Cell_Position = True

# Model dimensions
Xdim = 600
Ydim = 120
Zdim = 1

# Initial values for diffusion fields
Initial_Nutrient = 150.0
Minimum_Nutrient = 50
Initial_Oxygen = 50.0

# Agar thickness in pixels
Agar_Thickness = 20
Agar_Cell_Width = 50

# ECM parameters
InitEcm_targetVolume = 10
ECM_lambdaVolume = 10.0
ECM_secretion_time = 95

# Initial yeast cell placements
Start_Point = 187
End_Point = 412
Cell_Number = 32
S1_to_S2_Ratio = 1.0 # between 0.0 and 1.0
'''****************************************************************************************************************************************'''
# Shared parameters by both strains

# Parameters influencing the rate of cell targetVolume increase
Conversion_Rate = 0.0004
Growth_Rate = 0.3

# Ratio of active/passive cell metabolism
AP_Metabolism_Ratio = 0.5


# Time spent under threshold before cell type switchting
Threshold_Time_Limit = 20

# Initial target attributes
InitActiveCell_targetVolume = 20

InitCell_lambdaVolume = 80.0
Cell_Target_Pressure = 40

# Cell division volume dependence
CellDiv_Volume_Ratio = 1.5
'''*****************************************************************************************************************************************'''
# Strain Parameters

# Strain 1

# How much dose the current pressure of the cell influence the targetVolume increase
S1_Pressure_Dependence = 1.4

# How high the contact enregy goes with agar with time
S1_Max_Agar_Contact = 18

# Cell type switch thresholds

# Thresholds for Active/Passive cell status swtiching
S1_AP_Nutrient_Threshold = 9.0

# Threshold for Active/Fermenting cell status switching
S1_AF_Oxygen_Threshold = 8.0

# Threshold for Fermenting/Passive cell states switching
S1_FP_Nutrient_Threshold = 6.0

# Threshold for Passive/Dead cell status switchting
S1_PD_Oxygen_Threshold = 4.0
S1_PD_Nutrient_Threshold = 4.0
'''*****************************************************************************************************************************************'''