from cc3d import CompuCellSetup



from TwoStrainModel_2DColonyDevelopmentSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from TwoStrainModel_2DColonyDevelopmentSteppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from TwoStrainModel_2DColonyDevelopmentSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from TwoStrainModel_2DColonyDevelopmentSteppables import DeathSteppable

CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=1))




from TwoStrainModel_2DColonyDevelopmentSteppables import AgarParameterChangesSteppable

CompuCellSetup.register_steppable(steppable=AgarParameterChangesSteppable(frequency=1))


 
 
from TwoStrainModel_2DColonyDevelopmentSteppables import UptakeSteppable

CompuCellSetup.register_steppable(steppable=UptakeSteppable(frequency=1))


   

from TwoStrainModel_2DColonyDevelopmentSteppables import TypeSwitchSteppable

CompuCellSetup.register_steppable(steppable=TypeSwitchSteppable(frequency=1))




from TwoStrainModel_2DColonyDevelopmentSteppables import DataSaveSteppable

CompuCellSetup.register_steppable(steppable=DataSaveSteppable(frequency=1))

from InitialConcentrations import ConcentrationInit


CompuCellSetup.run()
