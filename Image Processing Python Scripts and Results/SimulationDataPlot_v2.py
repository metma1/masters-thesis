import pandas as pd
import numpy as np
import tkinter.filedialog
import os
import seaborn as sns
import matplotlib.pyplot as plt
from natsort import natsorted

pixel_to_micron = 23.20267324
mcs_to_day = 190
sec_pixel_in_sim_pixel = 14.39014713


sns.color_palette("bright")

# Re-plots the recorded colony dimension changes in time
def SimulationDimensionPlot():
    tkinter.Tk().withdraw()

    file = tkinter.filedialog.askopenfile(mode='r', filetypes=[('Text Files', '*.txt')], title="Select the ColonyDimensions file")

    if file:
        filepath = os.path.abspath(file.name)
        output_directory = filepath.split("ColonyDimensions.txt")[0]
    
    if 'S1' in output_directory:
        strain_name = 'S1 strain'
    elif 'S2' in output_directory:
        strain_name = 'S2 strain'
    else:
        strain_name = 'Mixed strains'

    simulation_dimensions_data = pd.read_csv(filepath, skipinitialspace = True)
    height_start = simulation_dimensions_data[simulation_dimensions_data['Width_x'] == 'Height_x'].index.values
    
    simulation_dimensions = pd.DataFrame()
    simulation_dimensions['X'] = simulation_dimensions_data['Width_x'][0:height_start[0]]
    simulation_dimensions['Width'] = simulation_dimensions_data['Width_y'][0:height_start[0]]
    array = simulation_dimensions_data['Width_y'][height_start[0]+1:].to_numpy()
    simulation_dimensions['Height'] = array 

    simulation_dimensions['Width'] = simulation_dimensions['Width'].apply(float)
    simulation_dimensions['Width'] = simulation_dimensions['Width']*pixel_to_micron
    simulation_dimensions['Height'] = simulation_dimensions['Height'].apply(float)
    simulation_dimensions['Height'] = simulation_dimensions['Height']*pixel_to_micron
    simulation_dimensions['X'] = simulation_dimensions['X'].apply(float)/mcs_to_day + 1

    plot = sns.lineplot(x='X', y='value', hue='variable', legend='full', data=pd.melt(simulation_dimensions, ['X']))
    plot.set(xlabel='Days', ylabel='Microns', title='Dimensions of ' + strain_name + ' colony simulation')
    plot.legend(title='')
    plt.grid()

    fig = plot.get_figure()
    fig.savefig(output_directory + "SimulationDimensions.png")





# Calculates the colony dimensions per row and column
def ColonyDimensions(df: pd.DataFrame):
    # Rotate the array into the correct position

    # Column
    # Find the last non-zero value in each column
    Y_colony_end = df.ne(0).idxmax(axis=0)
    # Find the first non-zero value in each column
    ud_flipped_df = df[::-1]
    ud_flipped_df.reset_index(drop=True)
    Y_colony_start = ud_flipped_df.ne(0).idxmax(axis=0)


    # Row
    # Find the first non-zero value in each row
    X_colony_start = df.ne(0).idxmax(axis=1)
    # Find the last non-zero value in each row
    lr_flipped_df = df[df.columns[::-1]]
    X_colony_end = lr_flipped_df.ne(0).idxmax(axis=1)
    

    # Calculate the heights of the colony
    Y_colony_end = Y_colony_end.astype('int')
    Y_colony_start += 1
    Y_colony_start = Y_colony_start.astype('int')
    Y_length_of_colony = Y_colony_start - Y_colony_end
    Y_length_of_colony[Y_length_of_colony == df.shape[0]] = 0
    height = Y_length_of_colony.max()

    # Calculate the widths of the colony
    X_colony_end = X_colony_end.astype('int')
    X_colony_end += 1
    X_colony_start = X_colony_start.astype('int')
    X_length_of_colony = X_colony_end - X_colony_start
    X_length_of_colony[X_length_of_colony == df.shape[1]] = 0
    width = X_length_of_colony.max()

    return [X_length_of_colony, Y_length_of_colony, width, height]





# Calculates and plots all descriptive metrics
def BatchSimulationPlot():
    tkinter.Tk().withdraw()
    # Get the directory of the pixel data files
    directory_path = tkinter.filedialog.askdirectory(title="Select the input pixel data directory")
    # Get the file list in the directory
    files = natsorted(os.listdir(directory_path))
    
    # Check if a Plots directory already exists
    if not os.path.exists(directory_path + "/Plots"):
        # If not, create one
        os.makedirs(directory_path + "/Plots")

    # Check if a NumericalFeatures directory already exists
    if not os.path.exists(directory_path + "/NumericalFeatures"):
        # If not, create one
        os.makedirs(directory_path + "/NumericalFeatures")
    
    # Create dataframe to contain the mean intermixing values per axes
    Mean_Intermixing = pd.DataFrame(columns=["MCS", "Column-wise", "Row-wise"])

    # Create dataframe to contain the strain ratios
    Strain_Ratio = pd.DataFrame(columns=["MCS", "S1", "S2"])

    # Iterate over all the files
    for file in files:
        if os.path.isfile(directory_path + "/" + file):
            print("Processing ", file)
            # Convert it to a Dataframe
            SegmentedSection = pd.read_csv(directory_path + "/" + file, delimiter=' ')

            # Rotate the dataframe to get the correct physical orientation
            df_array = SegmentedSection.to_numpy()
            df_array_rotated = np.rot90(df_array, axes=(1,0))
            SegmentedSection = pd.DataFrame(df_array_rotated, columns=SegmentedSection.index)

            # Get the widths and heights along the axes
            colony_widths, colony_heights, max_width, max_height = ColonyDimensions(SegmentedSection)

            # Strain ratios

            # Counting the number of the two strains per columns
            Strain1Num_perColumn = np.count_nonzero(SegmentedSection == -1, axis=0)
            Strain2Num_perColumn = np.count_nonzero(SegmentedSection == 1, axis=0)

            # Writing them into a dataframe
            StrainNums_perColumn = pd.DataFrame(columns=["Index", "Channel", "Value"])
            for index in range(len(Strain1Num_perColumn)):
                if (Strain1Num_perColumn[index] != 0 or Strain2Num_perColumn[index] != 0):
                    StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Green", Strain1Num_perColumn[index]/max_height]
                    StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perColumn[index]/max_height]
                else:
                    StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Green", Strain1Num_perColumn[index]]
                    StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perColumn[index]]
            
            # Reverse the order of the Value column (because the (0,0) point is on the bottom left of the simulation field) so the plot is oriented correctly     
            StrainNums_perColumn["Value"] = StrainNums_perColumn["Value"].values[::-1]
            StrainNums_perColumn["Channel"] = StrainNums_perColumn["Channel"].values[::-1]

            # Create a plot of the color changes along the x axis
            px = 1/plt.rcParams['figure.dpi']
            fig, ax = plt.subplots(figsize=(2000*px, 500*px))
            fig.tight_layout()
            plot = sns.lineplot(ax=ax, data=StrainNums_perColumn, x="Index", y="Value", hue="Channel", palette=["m", "g"])
            plot.set(xlabel='X (micron)', ylabel='Ratio', title='Strain ratio change along the x axis per channel')
            ax.set_xlim(left=0, right=len(Strain1Num_perColumn)*pixel_to_micron)

            # Save the plot as a png
            fig = plot.get_figure()
            mcs = file.split('.')
            fig.savefig(directory_path + "/Plots/" + mcs[0] + "MCS_StrainNumsPerColumns.png", bbox_inches='tight')
            plt.close()

            # Counting the number of the two strains per row
            Strain1Num_perRow = np.count_nonzero(SegmentedSection == -1, axis=1)
            Strain2Num_perRow = np.count_nonzero(SegmentedSection == 1, axis=1)

            # Writing them into a dataframe
            StrainNums_perRow = pd.DataFrame(columns=["Index", "Channel", "Value"])
            for index in range(len(Strain1Num_perRow)):
                if (Strain1Num_perRow[index] != 0 or Strain2Num_perRow[index] != 0):
                    StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Green", Strain1Num_perRow[index]/max_width]
                    StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perRow[index]/max_width]
                else:
                    StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Green", Strain1Num_perRow[index]]
                    StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perRow[index]]

            # Reverse the order of the Value column (because the (0,0) point is on the bottom left of the simulation filed) so the plot is oriented correctly     
            StrainNums_perRow["Value"] = StrainNums_perRow["Value"].values[::-1]
            StrainNums_perRow["Channel"] = StrainNums_perRow["Channel"].values[::-1]

            # Create a plot of the color changes along the y axis
            fig, ax = plt.subplots(figsize=(1000*px, 500*px))
            fig.tight_layout()
            plot = sns.lineplot(ax=ax, data=StrainNums_perRow, x="Index", y="Value", hue="Channel", palette=["m", "g"])
            plot.set(xlabel='Y (micron)', ylabel='Ratio', title='Strain ratio change along the y axis per channel')
            ax.set_xlim(left=0, right=len(Strain1Num_perRow)*pixel_to_micron)

            # Save it as a png
            fig = plot.get_figure()
            fig.savefig(directory_path + "/Plots/" + mcs[0] +  "MCS_StrainNumsPerRows.png", bbox_inches='tight')
            plt.close()

            # Calculate the mean strain ratios for the whole colony
            S1_Ratio = StrainNums_perRow[StrainNums_perRow["Channel"] == "Green"]["Value"].sum()
            S2_Ratio = StrainNums_perRow[StrainNums_perRow["Channel"] == "Magenta"]["Value"].sum()
            Strain_Ratio.loc[len(Mean_Intermixing.index)] = [mcs[0], S1_Ratio/(S1_Ratio+S2_Ratio), S2_Ratio/(S1_Ratio+S2_Ratio)]

            # Intermixing

            # Find all the indicies per column where the value is not 0
            Nonzero_Indices = SegmentedSection.apply(np.flatnonzero, axis=0)
            ColorChangeNum_perColumn = []
            for column in range(len(Nonzero_Indices)):
                ColorChange = 0
                for row in Nonzero_Indices[column][:-1]:
                    if SegmentedSection.iloc[row, column]*(-1) == SegmentedSection.iloc[row + 1, column]:
                        ColorChange += 1
                ColorChangeNum_perColumn.append(ColorChange)

            # Writing the number of changes per column into a dataframe
            ColorChangeNum_perColumn_df = pd.DataFrame(columns=["Index", "Value", "Thickness"])
            for index in range(len(ColorChangeNum_perColumn)):
                ColorChangeNum_perColumn_df.loc[len(ColorChangeNum_perColumn_df.index)] = [index*pixel_to_micron, ColorChangeNum_perColumn[index], colony_heights[index]] #
            
            # Subtract one from each element in the Thickness column if the value is not zero
            NonZeroValues = ColorChangeNum_perColumn_df["Thickness"] != 0
            MaxColorChange = ColorChangeNum_perColumn_df["Thickness"] - NonZeroValues
            
            # Weigh the color change numbers with the maximum possible color changes in the column
            ColorChangeNum_perColumn_df["Intermixing"] = ColorChangeNum_perColumn_df["Value"].div(MaxColorChange)
            ColorChangeNum_perColumn_df["Intermixing"].fillna(0, inplace=True)
            # Scaling the Intermixing values according to how many section pixels would be in a simulation pixel
            ColorChangeNum_perColumn_df["Intermixing"] /= sec_pixel_in_sim_pixel

            # Reverse the order of the Intermixing column (because the (0,0) point is on the bottom left of the simulation filed) so the plot is oriented correctly
            ColorChangeNum_perColumn_df["Intermixing"] = ColorChangeNum_perColumn_df["Intermixing"].values[::-1]

            # Create a plot from it
            fig, ax = plt.subplots(figsize=(2000*px, 500*px))
            fig.tight_layout()
            plot = sns.lineplot(ax=ax, data=ColorChangeNum_perColumn_df, x="Index", y="Intermixing")
            plot.set(xlabel='X (micron)', ylabel='Intermixing', title='Intermixing along the x axis')
            ax.set_xlim(left=0, right=len(Nonzero_Indices)*pixel_to_micron)

            # Save it as a png file
            fig = plot.get_figure()
            fig.savefig(directory_path + "/Plots/" + mcs[0] + "MCS_IntermixingPerColumns.png", bbox_inches='tight')
            plt.close()

            # Find all the indicies per row where the value is not 0
            Nonzero_Indices = SegmentedSection.apply(np.flatnonzero, axis=1)
            ColorChangeNum_perRow = []
            for row in range(len(Nonzero_Indices)):
                ColorChange = 0
                for column in Nonzero_Indices[row][:-1]:
                    if SegmentedSection.iloc[row, column]*(-1) == SegmentedSection.iloc[row, column + 1]:
                        ColorChange += 1
                ColorChangeNum_perRow.append(ColorChange)
            
            # Writing the number of changes per row into a dataframe
            ColorChangeNum_perRow_df = pd.DataFrame(columns=["Index", "Value", "Thickness"])
            for index in range(len(ColorChangeNum_perRow)):
                ColorChangeNum_perRow_df.loc[len(ColorChangeNum_perRow_df.index)] = [index*pixel_to_micron, ColorChangeNum_perRow[index], colony_widths[index]] #
                    

            # Subtract one from each element in the Thickness column if the value is not zero
            NonZeroValues = ColorChangeNum_perRow_df["Thickness"] != 0
            MaxColorChange = ColorChangeNum_perRow_df["Thickness"] - NonZeroValues

            # Weigh the color change numbers with the maximum possible color changes in the row
            ColorChangeNum_perRow_df["Intermixing"] = ColorChangeNum_perRow_df["Value"].div(MaxColorChange)
            ColorChangeNum_perRow_df["Intermixing"].fillna(0, inplace=True)
            # Scaling the Intermixing values according to how many section pixels would be in a simulation pixel
            ColorChangeNum_perRow_df["Intermixing"] /= sec_pixel_in_sim_pixel

            # Reverse the order of the Value column (because the (0,0) point is on the bottom left of the simulation filed) so the plot is oriented correctly
            ColorChangeNum_perRow_df["Intermixing"] = ColorChangeNum_perRow_df["Intermixing"].values[::-1]
            
            # Create a plot from it
            fig, ax = plt.subplots(figsize=(1000*px, 500*px))
            fig.tight_layout()
            plot = sns.lineplot(ax=ax, data=ColorChangeNum_perRow_df, x="Index", y="Intermixing")
            plot.set(xlabel='Y (micron)', ylabel='Intermixing', title='Intermixing along the y axis')
            ax.set_xlim(left=0, right=len(Nonzero_Indices)*pixel_to_micron)

            # Save it as a png
            fig = plot.get_figure()
            fig.savefig(directory_path + "/Plots/" + mcs[0] + "MCS_IntermixingPerRows.png", bbox_inches='tight')
            plt.close()

            # Calculate the mean intermixing along the two axes
            Columns_Mean_Intermixing = ColorChangeNum_perColumn_df["Intermixing"].mean()
            Rows_Mean_Intermixing = ColorChangeNum_perRow_df["Intermixing"].mean()
            Mean_Intermixing.loc[len(Mean_Intermixing.index)] = [mcs[0], Columns_Mean_Intermixing, Rows_Mean_Intermixing]

    # Save the strain ratios of the whole colony as a csv
    Strain_Ratio.to_csv(directory_path + "/NumericalFeatures/StrainRatio.csv", sep="\t", index=False)

    # Save the mean intermixing values as a csv
    Mean_Intermixing.to_csv(directory_path + "/NumericalFeatures/IntermixingIndex.csv", sep="\t", index=False)


# SimulationDimensionPlot() # Re-plotting the colony dimensions plot so that it has the same plot style as all the others
BatchSimulationPlot()