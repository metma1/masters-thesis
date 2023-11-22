from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np
from os.path import exists
import os
from math import ceil
from pathlib import Path
import shutil
import time

import parameters as p
from gif_maker import make_gif


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        for x in range(0, p.Xdim):
            for y in range(0, p.Agar_Thickness):
                AgarCell = self.new_cell(self.AGAR)
                self.cell_field[x, y, 0] = AgarCell
        
        # Gather all agar cell in a single cluster       
        Agcell0 = next(iter(self.cell_list_by_type(self.AGAR)))
        cluster_id = Agcell0.clusterId
        for cell in self.cell_list_by_type(self.AGAR):
            self.reassign_cluster_id(cell, cluster_id)
        
        # Placing initial yeast cells
        # cell_middle_arr = np.empty(shape=[p.Cell_Number, 1], dtype=int)
        cell_middle_list = np.empty(shape=[p.Cell_Number, 1], dtype=int)
        
        # Define the starting zone length
        start_zone_width = p.End_Point - p.Start_Point
        
        # Space between cells for even placement
        stagger = round(start_zone_width/p.Cell_Number)
        
        # Generate x coordinate for cell middles
        cell_middle_arr = np.array([p.Start_Point + n*stagger for n in range(p.Cell_Number)])
        cell_middle_list = cell_middle_arr.tolist()
        
        # Random distribution of cell types
        if p.Random_Cell_Position: 
            # Generate list with appropriate number of cells from each strain    
            n_S1 = round(p.S1_to_S2_Ratio*p.Cell_Number)
            cell_strain_arr = np.array([self.S1_ACTIVE]*n_S1 + [self.S2_ACTIVE]*(p.Cell_Number-n_S1))
            # Shuffle the order of the list
            np.random.shuffle(cell_strain_arr)
        
        # Manually assigned cell type per placement
        else:
            cell_strain_arr = np.array([self.S2_ACTIVE,
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE,
                                        self.S2_ACTIVE,
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE,
                                        self.S2_ACTIVE,
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE,
                                        self.S2_ACTIVE,
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE, 
                                        self.S2_ACTIVE, 
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE,
                                        self.S2_ACTIVE,
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE,
                                        self.S2_ACTIVE,
                                        self.S1_ACTIVE, self.S1_ACTIVE, self.S1_ACTIVE,
                                        self.S2_ACTIVE,
                                        self.S1_ACTIVE])
            
        # Assign shuffled cell types for cell field placements
        for n in range(p.Cell_Number):
            current_cell = self.new_cell(int(cell_strain_arr[n]))
            current_cell_middle = cell_middle_list[n]
            
            self.cell_field[current_cell_middle:current_cell_middle + 1, p.Agar_Thickness:p.Agar_Thickness + 1, 0] = current_cell


        #Initial cell attribute values
        for cell in self.cell_list:
            cell.dict['CellAge']=0
            cell.dict['Pressure']=0.0
            cell.dict['OxygenContent'] = 0.0
            cell.dict['NutrientContent'] = 0.0
            cell.dict['CellDivision'] = 0.0
            cell.dict['ECMProduction'] = 0.0
            cell.dict['TimeSinceLastDivision'] = 0
            
            # Strain 1
            if cell.type == self.S1_ACTIVE:
                cell.dict['OxygenContent'] = 20.0
                cell.dict['NutrientContent'] = 20.0
                cell.targetVolume = 1 # p.InitActiveCell_targetVolume 
                cell.lambdaVolume = p.InitCell_lambdaVolume
                cell.dict['TimeUnderThreshold'] = 0
                cell.dict['CellDivision'] = 0
                cell.dict['ECMProduction'] = 0
                cell.dict['NumberOfDivisions'] = 0
            
            # Strain 2
            if cell.type == self.S2_ACTIVE:
                cell.dict['OxygenContent'] = 20.0
                cell.dict['NutrientContent'] = 20.0
                cell.targetVolume = 1
                cell.lambdaVolume = p.InitCell_lambdaVolume
                cell.dict['TimeUnderThreshold'] = 0
                cell.dict['CellDivision'] = 0
                cell.dict['ECMProduction'] = 0
                cell.dict['NumberOfDivisions'] = 0
            
            if cell.type == self.AGAR:
                cell.lambdaVolume = 300.0
                cell.targetVolume = 1


class AgarParameterChangesSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def step(self, mcs):
        
        # Decrease constant nutrient amount at borders
        if mcs%5000 == 0 and p.Decreasing_Nutrient:
            Nutrient_X_Min = self.get_xml_element('Nutri_X_Min')
            if float(Nutrient_X_Min.Value) > p.Minimum_Nutrient:
                Nutrient_X_Min.Value = float(Nutrient_X_Min.Value) - 25
            
            Nutrient_X_Max = self.get_xml_element('Nutri_X_Max')
            if float(Nutrient_X_Max.Value) > p.Minimum_Nutrient:
                Nutrient_X_Max.Value = float(Nutrient_X_Max.Value) - 25
            
            Nutrient_Y_Min = self.get_xml_element('Nutri_Y_Min')
            if float(Nutrient_Y_Min.Value) > p.Minimum_Nutrient:
                Nutrient_Y_Min.Value = float(Nutrient_Y_Min.Value) - 25
        
        # Increase the contact energy between yeast cells and agar
        if mcs%95 == 0 and p.Agar_Drying:
            agar_s1a_energy = self.get_xml_element('agar_s1a_energy')
            if float(agar_s1a_energy.cdata) < p.S1_Max_Agar_Contact:
                agar_s1a_energy.cdata = float(agar_s1a_energy.cdata) + 1
  
            agar_s1p_energy = self.get_xml_element('agar_s1p_energy')
            if float(agar_s1p_energy.cdata) < p.S1_Max_Agar_Contact:
                agar_s1p_energy.cdata = float(agar_s1p_energy.cdata) + 1
            
            agar_s1d_energy = self.get_xml_element('agar_s1d_energy')
            if float(agar_s1d_energy.cdata) < p.S1_Max_Agar_Contact:
                agar_s1d_energy.cdata = float(agar_s1d_energy.cdata) + 1
            
            # agar_s1f_energy = self.get_xml_element('agar_s1f_energy')
            # if float(agar_s1f_energy.cdata) < 20:
                # agar_s1f_energy.cdata = float(agar_s1f_energy.cdata) + 1
                
            agar_s2a_energy = self.get_xml_element('agar_s2a_energy')
            if float(agar_s2a_energy.cdata) < p.S2_Max_Agar_Contact:
                agar_s2a_energy.cdata = float(agar_s2a_energy.cdata) + 1
            
            agar_s2p_energy = self.get_xml_element('agar_s2p_energy')
            if float(agar_s2p_energy.cdata) < p.S2_Max_Agar_Contact:
                agar_s2p_energy.cdata = float(agar_s2p_energy.cdata) + 1

            agar_s2d_energy = self.get_xml_element('agar_s2d_energy')
            if float(agar_s2d_energy.cdata) < p.S2_Max_Agar_Contact:
                agar_s2d_energy.cdata = float(agar_s2d_energy.cdata) + 1

            # agar_s2f_energy = self.get_xml_element('agar_s2f_energy')
            # if float(agar_s2f_energy.cdata) < 20:
                # agar_s2f_energy.cdata = float(agar_s2f_energy.cdata) + 1


class UptakeSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SecretionBasePy.__init__(self, frequency)

    def step(self, mcs):
        nutrient_field = self.get_field_secretor('Nutrient')      
        oxygen_field = self.get_field_secretor('Oxygen')

        # Oxygen, nutrient uptake for passive, active, fermenting cell types and nutrient release for dead cells
        
        # Strain 1
        for cell in self.cell_list_by_type(self.S1_ACTIVE):
            N = nutrient_field.uptakeInsideCellTotalCount(cell, 0.8, 0.8)
            O = oxygen_field.uptakeInsideCellTotalCount(cell,  0.8, 0.8)
            cell.dict['NutrientContent']+= np.abs(N.tot_amount)/cell.volume                    
            cell.dict['OxygenContent']+= np.abs(O.tot_amount)/cell.volume
              
        for cell in self.cell_list_by_type(self.S1_PASSIVE):
            N = nutrient_field.uptakeInsideCellTotalCount(cell, 0.5, 0.8)
            O = oxygen_field.uptakeInsideCellTotalCount(cell, 0.5, 0.8)
            cell.dict['NutrientContent']+= np.abs(N.tot_amount)/cell.volume                    
            cell.dict['OxygenContent']+= np.abs(O.tot_amount)/cell.volume
            
        for cell in self.cell_list_by_type(self.S1_FERMENT):
            N = nutrient_field.uptakeInsideCellTotalCount(cell, 0.5, 0.8)
            cell.dict['NutrientContent']+= np.abs(N.tot_amount)/cell.volume                  
            
        for cell in self.cell_list_by_type(self.S1_DEAD):
            nutrient_field.secreteOutsideCellAtBoundary(cell, cell.dict['NutrientContent'])
            cell.dict['NutrientContent'] = 0.0
            
        #Strain 2
        for cell in self.cell_list_by_type(self.S2_ACTIVE):
            N = nutrient_field.uptakeInsideCellTotalCount(cell, 0.8, 0.8)
            O = oxygen_field.uptakeInsideCellTotalCount(cell,  0.8, 0.8)
            cell.dict['NutrientContent']+= np.abs(N.tot_amount)/cell.volume                    
            cell.dict['OxygenContent']+= np.abs(O.tot_amount)/cell.volume
              
        for cell in self.cell_list_by_type(self.S2_PASSIVE):
            N = nutrient_field.uptakeInsideCellTotalCount(cell, 0.5, 0.8)
            O = oxygen_field.uptakeInsideCellTotalCount(cell, 0.5, 0.8)
            cell.dict['NutrientContent']+= np.abs(N.tot_amount)/cell.volume                    
            cell.dict['OxygenContent']+= np.abs(O.tot_amount)/cell.volume

        for cell in self.cell_list_by_type(self.S2_FERMENT):
            N = nutrient_field.uptakeInsideCellTotalCount(cell, 0.5, 0.8)
            cell.dict['NutrientContent']+= np.abs(N.tot_amount)/cell.volume
            
        for cell in self.cell_list_by_type(self.S2_DEAD):
            nutrient_field.secreteOutsideCellAtBoundary(cell, cell.dict['NutrientContent'])
            cell.dict['NutrientContent'] = 0.0
            
            
        # Increase Dead cells cell age, decrease cell targetVolume, secrete nutrient into surrounding
        for cell in self.cell_list_by_type(self.S1_DEAD, self.S2_DEAD):
            cell.dict['CellAge'] += 1
            if cell.dict['CellAge']%50 == 0 and cell.targetVolume > 10:
                cell.targetVolume -= 1
                nutrient_field.secreteOutsideCellAtBoundary(cell, 3)


class TypeSwitchSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SecretionBasePy.__init__(self, frequency)

    def step(self, mcs):
        # Threshold detection
        for cell in self.cell_list_by_type(self.S1_ACTIVE, self.S1_PASSIVE, self.S1_DEAD, self.S1_FERMENT, self.S2_ACTIVE, self.S2_PASSIVE, self.S2_DEAD, self.S2_FERMENT):
            # Strain 1
            if cell.type == self.S1_ACTIVE and p.S1_AP_Nutrient_Threshold > cell.dict['NutrientContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S1_PASSIVE
                    cell.dict['TimeUnderThreshold'] = 0
            
            if cell.type == self.S1_PASSIVE and p.S1_AP_Nutrient_Threshold <= cell.dict['NutrientContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S1_ACTIVE
                    cell.dict['TimeUnderThreshold'] = 0
                    
            if cell.type == self.S1_ACTIVE and p.S1_AF_Oxygen_Threshold  > cell.dict['OxygenContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S1_FERMENT
                    cell.dict['TimeUnderThreshold'] = 0
                    
            if cell.type == self.S1_FERMENT and p.S1_AF_Oxygen_Threshold <= cell.dict['OxygenContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S1_ACTIVE
                    cell.dict['TimeUnderThreshold'] = 0
                    
            if cell.type == self.S1_FERMENT and  p.S1_FP_Nutrient_Threshold > cell.dict['NutrientContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S1_PASSIVE
                    cell.dict['TimeUnderThreshold'] = 0

            if cell.type == self.S1_PASSIVE and (p.S1_PD_Oxygen_Threshold > cell.dict['OxygenContent'] or p.S1_PD_Nutrient_Threshold > cell.dict['NutrientContent']):
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S1_DEAD
                    cell.dict['TimeUnderThreshold'] = 0
                    
            # Strain 2        
            if cell.type == self.S2_ACTIVE and p.S2_AP_Nutrient_Threshold > cell.dict['NutrientContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S2_PASSIVE
                    cell.dict['TimeUnderThreshold'] = 0
            
            if cell.type == self.S2_PASSIVE and p.S2_AP_Nutrient_Threshold <= cell.dict['NutrientContent']: 
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S2_ACTIVE
                    cell.dict['TimeUnderThreshold'] = 0
                    
            if cell.type == self.S2_ACTIVE and p.S2_AF_Oxygen_Threshold  > cell.dict['OxygenContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S2_FERMENT
                    cell.dict['TimeUnderThreshold'] = 0
                    
            if cell.type == self.S2_FERMENT and p.S2_AF_Oxygen_Threshold <= cell.dict['OxygenContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S2_ACTIVE
                    cell.dict['TimeUnderThreshold'] = 0
                    
            if cell.type == self.S2_FERMENT and p.S2_FP_Nutrient_Threshold > cell.dict['NutrientContent']:
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S2_PASSIVE
                    cell.dict['TimeUnderThreshold'] = 0
                
            if cell.type == self.S2_PASSIVE and (p.S2_PD_Oxygen_Threshold > cell.dict['OxygenContent'] or p.S2_PD_Nutrient_Threshold > cell.dict['NutrientContent']):
                cell.dict['TimeUnderThreshold'] += 1
                if cell.dict['TimeUnderThreshold'] >= p.Threshold_Time_Limit:
                    cell.type = self.S2_DEAD
                    cell.dict['TimeUnderThreshold'] = 0


class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
        self.track_cell_level_scalar_attribute(field_name="CellNutrient",attribute_name= 'NutrientContent')
        self.track_cell_level_scalar_attribute(field_name="CellOxygen",attribute_name= 'OxygenContent')
        self.track_cell_level_scalar_attribute(field_name ='Cell_Pressure', attribute_name = 'Pressure')

    def step(self, mcs):
  
        # Current cell pressure of the cell
        for cell in self.cell_list_by_type(self.S1_ACTIVE, self.S1_PASSIVE, self.S2_ACTIVE, self.S2_PASSIVE):
            
            cell.dict['Pressure'] = 2*cell.lambdaVolume*(cell.targetVolume-cell.volume)
        
        #Metabolism of nutrient and oxygen inside of cell
        
        #Strain 1
        for cell in self.cell_list_by_type(self.S1_ACTIVE):
            Metabolism_Rate = max(0,p.Conversion_Rate*cell.dict['NutrientContent']*cell.dict['OxygenContent'])
            
            cell.dict['NutrientContent'] -= max(0,Metabolism_Rate*cell.dict['NutrientContent'])
            cell.dict['OxygenContent'] -= max(0,Metabolism_Rate/2*cell.dict['OxygenContent'])
            cell.targetVolume += p.Growth_Rate*max(0,p.Cell_Target_Pressure - p.S1_Pressure_Dependence*cell.dict['Pressure'])*float(Metabolism_Rate)
            cell.dict['CellAge']+=1
            
        for cell in self.cell_list_by_type(self.S1_PASSIVE):
            Metabolism_Rate = max(0,p.AP_Metabolism_Ratio*p.Conversion_Rate*cell.dict['NutrientContent']*cell.dict['OxygenContent'])
            
            cell.dict['NutrientContent'] -= max(0,Metabolism_Rate*cell.dict['NutrientContent'])
            cell.dict['OxygenContent'] -= max(0,Metabolism_Rate/2*cell.dict['OxygenContent']) 
            cell.dict['CellAge']+=1
            
        for cell in self.cell_list_by_type(self.S1_FERMENT):
            Metabolism_Rate = max(0,p.AF_Metabolism_Ratio*p.Conversion_Rate*cell.dict['NutrientContent']*cell.dict['OxygenContent'])
            
            cell.dict['NutrientContent'] -= max(0,Metabolism_Rate*cell.dict['NutrientContent'])
            cell.dict['OxygenContent'] -= max(0,Metabolism_Rate/2*cell.dict['OxygenContent'])
            cell.targetVolume += p.Growth_Rate*max(0,p.Cell_Target_Pressure - p.S1_Pressure_Dependence*cell.dict['Pressure'])*float(Metabolism_Rate)
            cell.dict['CellAge']+=1
            

        #Strain 2
        for cell in self.cell_list_by_type(self.S2_ACTIVE):
            Metabolism_Rate = max(0,p.Conversion_Rate*cell.dict['NutrientContent']*cell.dict['OxygenContent'])
            
            cell.dict['NutrientContent'] -= max(0,Metabolism_Rate*cell.dict['NutrientContent'])
            cell.dict['OxygenContent'] -= max(0,Metabolism_Rate/2*cell.dict['OxygenContent'])
            cell.targetVolume += p.Growth_Rate*p.Growth_Rate_Modifier*max(0,p.Cell_Target_Pressure - p.S2_Pressure_Dependence/2*cell.dict['Pressure'])*float(Metabolism_Rate)
            cell.dict['CellAge']+=1
            
            
        for cell in self.cell_list_by_type(self.S2_PASSIVE):
            Metabolism_Rate = max(0,p.AP_Metabolism_Ratio*p.Conversion_Rate*cell.dict['NutrientContent']*cell.dict['OxygenContent'])
            
            cell.dict['NutrientContent'] -= max(0,Metabolism_Rate*cell.dict['NutrientContent'])
            cell.dict['OxygenContent'] -= max(0,Metabolism_Rate/2*cell.dict['OxygenContent']) 
            cell.dict['CellAge']+=1
        
        for cell in self.cell_list_by_type(self.S2_FERMENT):
            Metabolism_Rate = max(0,p.AF_Metabolism_Ratio*p.Conversion_Rate*cell.dict['NutrientContent']*cell.dict['OxygenContent'])
            
            cell.dict['NutrientContent'] -= max(0,Metabolism_Rate*cell.dict['NutrientContent'])
            cell.dict['OxygenContent'] -= max(0,Metabolism_Rate/2*cell.dict['OxygenContent']) 
            cell.targetVolume += p.Growth_Rate*p.Growth_Rate_Modifier*max(0,p.Cell_Target_Pressure - p.S2_Pressure_Dependence/2*cell.dict['Pressure'])*float(Metabolism_Rate)
            cell.dict['CellAge']+=1

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):
        # Strain 1
        S1_Division_At_Volume = p.CellDiv_Volume_Ratio*p.InitActiveCell_targetVolume
        S2_Division_At_Volume = p.CellDiv_Volume_Ratio*p.InitActiveCell_targetVolume
        
        Cells_at_division_state = []
        Cells_To_Divide=[]
        # Strain 1
        for cell in self.cell_list_by_type(self.S1_ACTIVE):
            if cell.dict['NutrientContent'] >= p.S1_AP_Nutrient_Threshold and cell.volume >= S1_Division_At_Volume:
                Cells_To_Divide.append(cell)
                cell.dict['CellDivision'] = 1
                cell.dict['NumberOfDivisions'] += 1
                
        for cell in self.cell_list_by_type(self.S1_FERMENT):
            if cell.volume >= S1_Division_At_Volume:
                Cells_To_Divide.append(cell)
                cell.dict['CellDivision'] = 1
                cell.dict['NumberOfDivisions'] += 1
  
        for cell in self.cell_list_by_type(self.S1_FERMENT):    
            if p.ECM_Production and cell.dict['CellAge']%p.ECM_secretion_time == 0:
                Cells_To_Divide.append(cell)
                cell.dict['ECMProduction'] = 1
        
        # Strain 2  
        for cell in self.cell_list_by_type(self.S2_ACTIVE):
            if cell.dict['NutrientContent'] >= p.S2_AP_Nutrient_Threshold and cell.volume > S2_Division_At_Volume:
                Cells_To_Divide.append(cell)
                cell.dict['CellDivision'] = 1
                cell.dict['NumberOfDivisions'] += 1
                
        for cell in self.cell_list_by_type(self.S2_FERMENT):
            if cell.volume > S2_Division_At_Volume:
                Cells_To_Divide.append(cell)
                cell.dict['CellDivision'] = 1
                cell.dict['NumberOfDivisions'] += 1
            
        for cell in self.cell_list_by_type(self.S2_FERMENT):    
            if p.ECM_Production and cell.dict['CellAge']%p.ECM_secretion_time == 0:
                Cells_To_Divide.append(cell)
                cell.dict['ECMProduction'] = 1
                
        for cell in Cells_To_Divide:
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # Strain 1
        if self.parent_cell.type == self.S1_ACTIVE:               
            
            # Mitosis event
            if self.parent_cell.dict['CellDivision'] == 1:
                # Putting the child cell into the same cluster
                self.reassign_cluster_id(self.child_cell, self.parent_cell.clusterId)
                
                self.clone_parent_2_child()
                
                self.child_cell.targetVolume = int((p.CellDiv_Volume_Ratio - 1) * p.InitActiveCell_targetVolume)
                self.parent_cell.targetVolume = p.InitActiveCell_targetVolume
                
                volumeRatio = self.child_cell.targetVolume/self.parent_cell.targetVolume
                
                self.child_cell.dict['NutrientContent']= volumeRatio*self.parent_cell.dict['NutrientContent']
                self.child_cell.dict['OxygenContent']= volumeRatio*self.parent_cell.dict['OxygenContent']
                
                self.parent_cell.dict['NutrientContent']-= max(0,self.child_cell.dict['NutrientContent'])
                self.parent_cell.dict['OxygenContent']-= max(0,self.child_cell.dict['OxygenContent'])
                
                self.parent_cell.dict['CellDivision'] = 0
            
                self.child_cell.dict['NumberOfDivisions'] = 0
                self.child_cell.dict['CellAge'] = 0
                self.child_cell.dict['NumberOfDivisions'] = 0
                
        # Strain 2
        if self.parent_cell.type == self.S2_ACTIVE:              
            
            # Mitosis event
            if self.parent_cell.dict['CellDivision'] == 1:
                # Putting the child cell into the same cluster
                self.reassign_cluster_id(self.child_cell, self.parent_cell.clusterId)
                
                self.clone_parent_2_child()
       
                self.child_cell.targetVolume = int((p.CellDiv_Volume_Ratio - 1) * p.InitActiveCell_targetVolume)
                self.parent_cell.targetVolume = p.InitActiveCell_targetVolume
                
                volumeRatio = self.child_cell.targetVolume/self.parent_cell.targetVolume
                
                self.child_cell.dict['NutrientContent']= volumeRatio*self.parent_cell.dict['NutrientContent']
                self.child_cell.dict['OxygenContent']= volumeRatio*self.parent_cell.dict['OxygenContent']
                
                self.parent_cell.dict['NutrientContent']-= max(0,self.child_cell.dict['NutrientContent'])
                self.parent_cell.dict['OxygenContent']-= max(0,self.child_cell.dict['OxygenContent'])
                
                self.parent_cell.dict['CellDivision'] = 0
                
                self.child_cell.dict['NumberOfDivisions'] = 0
                self.child_cell.dict['CellAge'] = 0
                self.child_cell.dict['NumberOfDivisions'] = 0

            
        # ECM secretion event    
        if self.parent_cell.dict['ECMProduction'] == 1:
            
            # Putting the child cell into the same cluster
            self.reassign_cluster_id(self.child_cell, self.parent_cell.clusterId)
            
            self.child_cell.type = self.ECM
            self.child_cell.dict['CellAge'] = 0
            self.child_cell.targetVolume = p.InitEcm_targetVolume
            self.child_cell.lambdaVolume = p.ECM_lambdaVolume
            self.parent_cell.dict['ECMProduction'] = 0


class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.S1_DEAD, self.S2_DEAD):
            cell.lambdaVolume = 100


class DataSaveSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        os.mkdir(self.output_dir + "\SavedData")
        os.mkdir(self.output_dir + "\SavedData\\NeighborData")
        os.mkdir(self.output_dir + "\SavedData\\PixelData")
        os.mkdir(self.output_dir + "\SavedData\\PlotData")
        # os.mkdir(self.output_dir + "\SavedData\\EdgeData")
        os.mkdir(self.output_dir + "\SavedData\\CellFieldData")
        
        
        # Plot for cell numbers in different cell states
        self.S1_cell_state_numbers = self.add_new_plot_window(title='Volume of Strain1 cells in each cell state',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Total cell volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        
        self.S1_cell_state_numbers.add_plot("S1 Active Cell", style='Lines', color='blue', size=3)
        self.S1_cell_state_numbers.add_plot("S1 Passive Cell", style='Lines', color='orange', size=3)
        self.S1_cell_state_numbers.add_plot("S1 Dead Cell", style='Lines', color='yellow', size=3)
        self.S1_cell_state_numbers.add_plot("S1 Fermenting Cell", style='Lines', color='green', size=3)
        
        self.S2_cell_state_numbers = self.add_new_plot_window(title='Volume of Strain2 cells in each cell state',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Total cell volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        
        self.S2_cell_state_numbers.add_plot("S2 Active Cell", style='Lines', color='blue', size=3)
        self.S2_cell_state_numbers.add_plot("S2 Passive Cell", style='Lines', color='orange', size=3)
        self.S2_cell_state_numbers.add_plot("S2 Dead Cell", style='Lines', color='yellow', size=3)
        self.S2_cell_state_numbers.add_plot("S2 Fermenting Cell", style='Lines', color='green', size=3)
        
        # Plot for colony dimensions
        self.colony_dims = self.add_new_plot_window(title='Colony dimensions',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Pixels', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        
        self.colony_dims.add_plot("Width", style='Lines', color='red', size=3)
        self.colony_dims.add_plot("Height", style='Lines', color='green', size=3)

    def step(self, mcs):
        if mcs%95 == 0:
            
            # Add data points to the cell pixel numbers in each cell state plot
            S1_Active_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S1_ACTIVE):
                S1_Active_Pixel_Num += cell.volume
            self.S1_cell_state_numbers.add_data_point("S1 Active Cell", mcs, S1_Active_Pixel_Num)
            
            S1_Passive_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S1_PASSIVE):
                S1_Passive_Pixel_Num += cell.volume
            self.S1_cell_state_numbers.add_data_point("S1 Passive Cell", mcs, S1_Passive_Pixel_Num)
            
            S1_Dead_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S1_DEAD):
                S1_Dead_Pixel_Num += cell.volume
            self.S1_cell_state_numbers.add_data_point("S1 Dead Cell", mcs, S1_Dead_Pixel_Num)
            
            S1_Ferment_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S1_FERMENT):
                S1_Ferment_Pixel_Num += cell.volume
            self.S1_cell_state_numbers.add_data_point("S1 Fermenting Cell", mcs, S1_Ferment_Pixel_Num)

            # Add data points to the cell pixel numbers in each cell state plot
            S2_Active_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S2_ACTIVE):
                S2_Active_Pixel_Num += cell.volume
            self.S2_cell_state_numbers.add_data_point("S2 Active Cell", mcs, S2_Active_Pixel_Num)
            
            S2_Passive_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S2_PASSIVE):
                S2_Passive_Pixel_Num += cell.volume
            self.S2_cell_state_numbers.add_data_point("S2 Passive Cell", mcs, S2_Passive_Pixel_Num)
            
            S2_Dead_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S2_DEAD):
                S2_Dead_Pixel_Num += cell.volume
            self.S2_cell_state_numbers.add_data_point("S2 Dead Cell", mcs, S2_Dead_Pixel_Num)
            
            S2_Ferment_Pixel_Num = 0
            for cell in self.cell_list_by_type(self.S2_FERMENT):
                S2_Ferment_Pixel_Num += cell.volume
            self.S2_cell_state_numbers.add_data_point("S2 Fermenting Cell", mcs, S2_Ferment_Pixel_Num)           
            
            # Find the current highest point of the colony
            colony_height_min = p.Ydim
            colony_height_max = 0
            for cell in self.cell_list:
                for pixel_data in self.get_cell_pixel_list(cell):
                    if pixel_data.pixel.y > colony_height_max:
                        colony_height_max = pixel_data.pixel.y
                    if pixel_data.pixel.y < colony_height_min:
                        colony_height_min = pixel_data.pixel.y
            
            # Find the current farthest points of the colony widthwise
            colony_width_min = p.Xdim
            colony_width_max = 0
            for cell in self.cell_list:
                if cell.type != self.AGAR:
                    for pixel_data in self.get_cell_pixel_list(cell):
                        if pixel_data.pixel.x > colony_width_max:
                            colony_width_max = pixel_data.pixel.x
                        if pixel_data.pixel.x < colony_width_min:
                            colony_width_min = pixel_data.pixel.x
            
            # Calculate maximum colony width and height
            colony_width = colony_width_max - colony_width_min
            colony_height = colony_height_max - p.Agar_Thickness
            # Add data point to the colony dimension plot
            self.colony_dims.add_data_point("Width", mcs, colony_width)
            self.colony_dims.add_data_point("Height", mcs, colony_height)

            
            # Save every pixels cell strain as a number 1 or -1, MEDIUM AGAR and ECM is 0

            # Save the current row of the cell field in a csv
            savedir = "SavedData\\PixelData\\"
            txt_path = Path(self.output_dir).joinpath(savedir + str(mcs) + ".txt")        
            with open(txt_path, "a") as pixelfile:
                pixelfile.write(' '.join(str(e) for e in [y for y in range(p.Ydim)]) + "\n")

            Column_data = []
            # Go trough every pixel
            for x,y,z in self.every_pixel():
                # If we are at the end of the column
                if (y == (p.Ydim-1)):
                    cell = self.cell_field[x,y,z]
                    # If it is a cell
                    if cell:  # this is true if the pixel is not MEDIUM
                        # If it is a cell from S1
                        if ((cell.type == 3) | (cell.type == 4) | (cell.type == 5) | (cell.type == 6)):
                            Column_data.append(-1)
                        # If it is a cell from S2
                        elif ((cell.type == 7) | (cell.type == 8) | (cell.type == 9) | (cell.type == 10)):
                            Column_data.append(1)
                        # If it is Agar or ECM
                        elif ((cell.type == 1) | (cell.type == 2)):
                            Column_data.append(0)
                    # If it is medium
                    else:
                        Column_data.append(0)
                        # Save the current row of the cell field in a csv
                        with open(txt_path, "a") as pixelfile:
                            pixelfile.write(' '.join(str(e) for e in Column_data) + "\n")
                        # Empty the Row_data array
                        Column_data = []
                else:
                    # If we are not at the end of the column
                    cell = self.cell_field[x,y,z]
                    # If it is a cell
                    if cell:  # this is true if the pixel is not MEDIUM
                        # If it is a cell from S1
                        if ((cell.type == 3) | (cell.type == 4) | (cell.type == 5) | (cell.type == 6)):
                            Column_data.append(-1)
                        # If it is a cell from S2
                        elif ((cell.type == 7) | (cell.type == 8) | (cell.type == 9) | (cell.type == 10)):
                            Column_data.append(1)
                        # If it is Agar or ECM
                        elif ((cell.type == 1) | (cell.type == 2)):
                            Column_data.append(0)
                    # If it is medium
                    else:
                        Column_data.append(0)

            # X
            # Contains each cell id along the current x axis sample
            sample_Xdimensions = []
            
            # Sample the cells at the current x coordinate
            for x in range(p.Xdim+1):
                
                # Contains the cell ids along the current axis
                sample_set_x = set()
                
                # Iterate over the current axis in the y direction
                for y in range(p.Ydim+1):
                    current_cell = self.cell_field[x, y, 0]
                    # Only if the current pixel is part of a cell (not MEDIUM or AGAR)
                    if current_cell:
                        if (current_cell.type != 1):
                            sample_set_x.add(current_cell.id)
                        
                # Append the current columns start and end points
                sample_Xdimensions.append(sample_set_x)
                
            # Iterate over the selected cells in X dim and count each neighbor by type
            sample_Xdim_print = 0
            for Width_cell_Set in sample_Xdimensions:
                for cell_id in Width_cell_Set:
                    cell = self.fetch_cell_by_id(cell_id)
                    Cell_neigbor_by_type = {"MCS": mcs, "SampleXdim": sample_Xdim_print, "CellID": cell.id, "CellType" : cell.type, "XCOM" : round(cell.xCOM,1), "YCOM" : round(cell.yCOM, 1), "S1_ACTIVE" : 0, "S1_PASSIVE" : 0, "S1_FERMENT" : 0, "S1_DEAD" : 0, "S2_ACTIVE" : 0, "S2_PASSIVE" : 0, "S2_FERMENT" : 0, "S2_DEAD" : 0, "ECM" : 0, "AGAR" : 0, "MEDIUM" : 0}
                    # Count how many neighbors from each cell type a cell has
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                        if neighbor:
                            if neighbor.type == self.S1_ACTIVE:
                                Cell_neigbor_by_type["S1_ACTIVE"] += 1
                            elif neighbor.type == self.S1_PASSIVE:
                                Cell_neigbor_by_type["S1_PASSIVE"] += 1
                            elif neighbor.type == self.S1_FERMENT:
                                Cell_neigbor_by_type["S1_FERMENT"] += 1
                            elif neighbor.type == self.S1_DEAD:
                                Cell_neigbor_by_type["S1_DEAD"] += 1
                            elif neighbor.type == self.S2_ACTIVE:
                                Cell_neigbor_by_type["S2_ACTIVE"] += 1
                            elif neighbor.type == self.S2_PASSIVE:
                                Cell_neigbor_by_type["S2_PASSIVE"] += 1
                            elif neighbor.type == self.S2_FERMENT:
                                Cell_neigbor_by_type["S1_FERMENT"] += 1
                            elif neighbor.type == self.S2_DEAD:
                                Cell_neigbor_by_type["S2_DEAD"] += 1
                            elif neighbor.type == self.ECM:
                                Cell_neigbor_by_type["ECM"] += 1
                            elif neighbor.type == self.AGAR:
                                Cell_neigbor_by_type["AGAR"] += 1
                        else: # Medium cell
                            Cell_neigbor_by_type["MEDIUM"] += 1
                    savedir = "SavedData\\NeighborData\\"
                    txt_path = Path(self.output_dir).joinpath(savedir + "CellNeighbors.txt")        
                    with open(txt_path, "a") as neighborfile:
                        neighborfile.write(str(Cell_neigbor_by_type) + "\n")
                sample_Xdim_print += 1 
    
    def finish(self): 
        savedir = "SavedData\\PlotData\\"
        
        # Save cell state number plot data as png and txt file
        txt_path = Path(self.output_dir).joinpath(savedir + "S1_CellStateNumbers.txt")
        self.S1_cell_state_numbers.save_plot_as_data(txt_path, CSV_FORMAT)
        self.S1_cell_state_numbers.save_plot_as_png(savedir + "S1_CellStateNumbers.png", 1000, 1000)
        
        txt_path = Path(self.output_dir).joinpath(savedir + "S2_CellStateNumbers.txt")
        self.S2_cell_state_numbers.save_plot_as_data(txt_path, CSV_FORMAT)
        self.S2_cell_state_numbers.save_plot_as_png(savedir + "S2_CellStateNumbers.png", 1000, 1000)
        
        # Save colony dimension plot data as png and txt file
        txt_path = Path(self.output_dir).joinpath(savedir + "ColonyDimensions.txt")
        self.colony_dims.save_plot_as_data(txt_path, CSV_FORMAT)
        self.colony_dims.save_plot_as_png(savedir + "ColonyDimensions.png", 1000, 1000)
        
        # Moving all cell field images to same directory
        source_path1 = Path(self.output_dir).joinpath("Cell_Field_CellField_2D_XY_0")
        source_path2 = Path(self.output_dir).joinpath("Cell_Pressure_ScalarFieldCellLevel_2D_XY_0")
        source_path3 = Path(self.output_dir).joinpath("CellNutrient_ScalarFieldCellLevel_2D_XY_0")
        source_path4 = Path(self.output_dir).joinpath("CellOxygen_ScalarFieldCellLevel_2D_XY_0")
        source_path5 = Path(self.output_dir).joinpath("Nutrient_ConField_2D_XY_0")
        source_path6 = Path(self.output_dir).joinpath("Oxygen_ConField_2D_XY_0")
        
        destination_path = Path(self.output_dir).joinpath("SavedData\\CellFieldData")
        
        shutil.move(source_path1, destination_path)
        shutil.move(source_path2, destination_path)
        shutil.move(source_path3, destination_path)
        shutil.move(source_path4, destination_path)
        shutil.move(source_path5, destination_path)
        shutil.move(source_path6, destination_path)
        
        time.sleep(1)
        
        make_gif(str(destination_path), "Cell_Field_CellField_2D_XY_0")
        make_gif(str(destination_path), "Cell_Pressure_ScalarFieldCellLevel_2D_XY_0")
        make_gif(str(destination_path), "CellNutrient_ScalarFieldCellLevel_2D_XY_0")
        make_gif(str(destination_path), "CellOxygen_ScalarFieldCellLevel_2D_XY_0")
        make_gif(str(destination_path), "Nutrient_ConField_2D_XY_0")
        make_gif(str(destination_path), "Oxygen_ConField_2D_XY_0")
        
        return

