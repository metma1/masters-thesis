# Repository structure

## CompuCell3D Models:
Contains the three different model projects, with all necessary files for running a simulation on a computer which has CompuCell3D installed.
CompuCell3D download page [link](https://compucell3d.org/SrcBin#MostRecent)

The model can be opened from the Tweedit++ editor by opening the .cc3d file with it.

## Image Processing Python Scripts and Results:
Contains all python scripts which were used to process all section images and simulation results and characteristic metric calculations and generates plots of them:

- SectionDataPlot_v2.py: Processes the colony sections.
- SimulationDataPlot_v2.py: Processes the simulation results.
- Section-Simulation_Comparison.ipynb: Plots the comparison metrics.
- Section-Simulation Comparison: Contains .csv files used by 'Section-Simulation_Comparison.ipynb' to plot the result comparisons.

## Section_Data_Gen.jim:
ImageJ Macro script, which can be run on a computer with ImageJ installed.
ImageJ download page: [link](https://imagej.net/ij/download.html)

Iterates over a whole strain folder and generates the necessary data from two .tif section images for further processing and groups them in a folder.


## Other recources
Google Drive folder with the strain folders containing example section images and model simulation results: [link](https://drive.google.com/drive/folders/1s5FahRhtLXPPjGtBNqd1eofk7OXduhBy?usp=sharing)
