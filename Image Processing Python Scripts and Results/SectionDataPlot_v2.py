import pandas as pd
import numpy as np
import tkinter
import tkinter.filedialog
import os
import seaborn as sns
import matplotlib.pyplot as plt

pixel_to_micron = 1.6124
sns.color_palette("bright")

# Finds the maximum width and height of the colony
def ColonyDimensions(df = pd.DataFrame):
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
    # Get the maximum height
    height = Y_length_of_colony.max()

    # Calculate the widths of the colony
    X_colony_end = X_colony_end.astype('int')
    X_colony_end += 1
    X_colony_start = X_colony_start.astype('int')
    X_length_of_colony = X_colony_end - X_colony_start
    X_length_of_colony[X_length_of_colony == df.shape[1]] = 0
    # Get the maximum width
    width = X_length_of_colony.max()

    return [X_length_of_colony, Y_length_of_colony, width, height]





# Counting and saving the number of pixels associated to the given strains per rows and columns, calculating and plotting the strain ratios
def SegmentationStrainRatio(SegmentedSection, colony_widths, colony_heights, output_path):

    # Counting the number of the two strains per columns
    Strain1Num_perColumn = np.count_nonzero(SegmentedSection == -1, axis=0)
    Strain2Num_perColumn = np.count_nonzero(SegmentedSection == 1, axis=0)

    # Writing them into a dataframe
    StrainNums_perColumn = pd.DataFrame(columns=["Index", "Channel", "Value"])
    for index in range(len(Strain1Num_perColumn)):
        # If there are cells present in the axes
        if (Strain1Num_perColumn[index] != 0 or Strain2Num_perColumn[index] != 0):
            StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Green", Strain1Num_perColumn[index]/colony_heights.max()]
            StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perColumn[index]/colony_heights.max()]
        # If there are no cells present
        else:
            StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Green", Strain1Num_perColumn[index]]
            StrainNums_perColumn.loc[len(StrainNums_perColumn.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perColumn[index]]

    # Create a plot of the color changes along the x axis
    px = 1/plt.rcParams['figure.dpi']
    fig, ax = plt.subplots(figsize=(2000*px, 800*px))
    # fig, ax = plt.subplots(figsize=(20, 8))
    plot = sns.lineplot(ax=ax, data=StrainNums_perColumn, x="Index", y="Value", hue="Channel", palette=["g", "m"])
    plot.set(xlabel='X (micron)', ylabel='Ratio', title='Strain ratio change along the x axis per channel')
    ax.set_xlim(left=0, right=len(Strain1Num_perColumn)*pixel_to_micron)

    # Check if a Plots directory already exists
    if not os.path.exists(output_path + "/Plots"):
        # If not, create one
        os.makedirs(output_path + "/Plots")

    # Save the plot as a png
    fig = plot.get_figure()
    fig.savefig(output_path + "/Plots/StrainNumsPerColumns.png")
    plt.close()

    # Counting the number of the two strains per row
    Strain1Num_perRow = np.count_nonzero(SegmentedSection == -1, axis=1)
    Strain2Num_perRow = np.count_nonzero(SegmentedSection == 1, axis=1)

    # Writing them into a dataframe
    StrainNums_perRow = pd.DataFrame(columns=["Index", "Channel", "Value"])
    for index in range(len(Strain1Num_perRow)):
        # If there are cells present in the axes
        if (Strain1Num_perRow[index] != 0 or Strain2Num_perRow[index] != 0):
            StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Green", Strain1Num_perRow[index]/colony_widths.max()]
            StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perRow[index]/colony_widths.max()]
        # If there are no cells present
        else:
            StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Green", Strain1Num_perRow[index]]
            StrainNums_perRow.loc[len(StrainNums_perRow.index)] = [index*pixel_to_micron, "Magenta", Strain2Num_perRow[index]]
  
    # Create a plot of the color changes along the y axis
    fig, ax = plt.subplots(figsize=(1000*px, 800*px))
    # fig, ax = plt.subplots(figsize=(20, 8))
    plot = sns.lineplot(ax=ax, data=StrainNums_perRow, x="Index", y="Value", hue="Channel", palette=["g", "m"])
    plot.set(xlabel='Y (micron)', ylabel='Ratio', title='Strain ratio change along the y axis per channel')
    ax.set_xlim(left=0, right=len(Strain1Num_perRow)*pixel_to_micron)

    # Save it as a png
    fig = plot.get_figure()
    fig.savefig(output_path + "/Plots/StrainNumsPerRows.png")
    plt.close()

    # Create dataframe to contain the strain ratios
    Strain_Ratio = pd.DataFrame(columns=["Green", "Magenta"])

    # Check if a NumericalFeatures directory already exists
    if not os.path.exists(output_path + "/NumericalFeatures"):
        # If not, create one
        os.makedirs(output_path + "/NumericalFeatures")

    # Calculate the strain ratio for the whole colony
    S1_Ratio = StrainNums_perRow[StrainNums_perRow["Channel"] == "Green"]["Value"].sum()
    S2_Ratio = StrainNums_perRow[StrainNums_perRow["Channel"] == "Magenta"]["Value"].sum()
    Strain_Ratio.loc[len(Strain_Ratio.index)] = [S1_Ratio/(S1_Ratio+S2_Ratio), S2_Ratio/(S1_Ratio+S2_Ratio)]
    Strain_Ratio.to_csv(output_path + "/NumericalFeatures/StrainRatio.csv", sep="\t", index=False)




# Counting and saving the number of color changes per rows and columns, calculating and plotting the intermixing indexes
def SegmentationIntermixing(SegmentedSection: pd.DataFrame, colony_widths, colony_heights, output_path):

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
        ColorChangeNum_perColumn_df.loc[len(ColorChangeNum_perColumn_df.index)] = [index*pixel_to_micron, ColorChangeNum_perColumn[index], colony_heights[index]]
    
    # Subtract one from each element in the Thickness column if the value is not zero
    NonZeroValues = ColorChangeNum_perColumn_df["Thickness"] != 0
    MaxColorChange = ColorChangeNum_perColumn_df["Thickness"] - NonZeroValues
    
    # Weigh the color change numbers with the maximum possible color changes in the column
    ColorChangeNum_perColumn_df["Intermixing"] = ColorChangeNum_perColumn_df["Value"].div(MaxColorChange)
    ColorChangeNum_perColumn_df["Intermixing"].fillna(0, inplace=True)

    # Create a plot from it
    px = 1/plt.rcParams['figure.dpi']
    fig, ax = plt.subplots(figsize=(2000*px, 800*px))
    # fig, ax = plt.subplots(figsize=(20, 8))
    plot = sns.lineplot(ax=ax, data=ColorChangeNum_perColumn_df, x="Index", y="Intermixing")
    plot.set(xlabel='X (micron)', ylabel='Intermixing', title='Intermixing along the x axis')
    ax.set_xlim(left=0, right=len(Nonzero_Indices)*pixel_to_micron)

    # Save it as a png file
    fig = plot.get_figure()
    fig.savefig(output_path + "/Plots/IntermixingPerColumns.png")
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
        ColorChangeNum_perRow_df.loc[len(ColorChangeNum_perRow_df.index)] = [index*pixel_to_micron, ColorChangeNum_perRow[index], colony_widths[index]]
            

    # Subtract one from each element in the Thickness column if the value is not zero
    NonZeroValues = ColorChangeNum_perRow_df["Thickness"] != 0
    MaxColorChange = ColorChangeNum_perRow_df["Thickness"] - NonZeroValues

    # Weigh the color change numbers with the maximum possible color changes in the row
    ColorChangeNum_perRow_df["Intermixing"] = ColorChangeNum_perRow_df["Value"].div(MaxColorChange)
    ColorChangeNum_perRow_df["Intermixing"].fillna(0, inplace=True)

    # Create a plot from it
    fig, ax = plt.subplots(figsize=(1000*px, 800*px))
    # fig, ax = plt.subplots(figsize=(20, 8))
    plot = sns.lineplot(ax=ax, data=ColorChangeNum_perRow_df, x="Index", y="Intermixing")
    plot.set(xlabel='Y (micron)', ylabel='Intermixing', title='Intermixing along the y axis')
    ax.set_xlim(left=0, right=len(Nonzero_Indices)*pixel_to_micron)

    # Save it as a png
    fig = plot.get_figure()
    fig.savefig(output_path + "/Plots/IntermixingPerRows.png")
    plt.close()

    # Create dataframe to contain the mean intermixing values per axes
    Mean_Intermixing = pd.DataFrame(columns=["Column-wise", "Row-wise"])

    # Check if a NumericalFeatures directory already exists
    if not os.path.exists(output_path + "/NumericalFeatures"):
        # If not, create one
        os.makedirs(output_path + "/NumericalFeatures")

    # Calculate the mean intermixing along the two axes
    Columns_Mean_Intermixing = ColorChangeNum_perColumn_df["Intermixing"].mean()
    Rows_Mean_Intermixing = ColorChangeNum_perRow_df["Intermixing"].mean()
    Mean_Intermixing.loc[len(Mean_Intermixing.index)] = [Columns_Mean_Intermixing, Rows_Mean_Intermixing]
    Mean_Intermixing.to_csv(output_path + "/NumericalFeatures/IntermixingIndex.csv", sep="\t", index=False)





# Segments the section pixel values according to the two channels pixel intensity at the given pixel
def SectionSegmentation(channel1: pd.DataFrame, channel2: pd.DataFrame, colony_widths, colony_heights, output_path, file_name):
    # SegmentedSectionDf_array = np.empty(channel1.shape)

    # Normalizing the intensity values    
    channel1_array = channel1.to_numpy()
    ch1_max = channel1_array.max()
    channel1_array = np.divide(channel1_array, ch1_max)
 
    # Normalizing the intensity values
    channel2_array = channel2.to_numpy()
    ch2_max = channel2_array.max()
    channel2_array = np.divide(channel2_array, ch2_max)

    # Check if the difference between the two channels is big enough
    channels_diff = abs(channel1_array - channel2_array)
    channel_diff_std = np.std(channels_diff)
    cleaned_channels = (channels_diff > 0.1*channel_diff_std).astype(int)

    # Replace the ambiguous values with 0-s
    channel1_array = np.multiply(channel1_array,cleaned_channels)
    channel2_array = np.multiply(channel2_array,cleaned_channels)

    # Check where the values in channel2 are higher than channel1
    channel2_higher = (channel1_array < channel2_array).astype(int)

    # Check where the values in channel1 are higher than channel2
    channel1_higher = (channel2_array < channel1_array).astype(int)
    channel1_higher = channel1_higher*-1

    # Add the two array together
    SegmentedSectionDf_array = channel1_higher + channel2_higher

    # Convert to dataframe
    SegmentedSectionDf = pd.DataFrame(SegmentedSectionDf_array, columns=channel1.columns)

    # Save the segmented section data as a csv file
    SegmentedSectionDf.to_csv(output_path + "/" + file_name + "_Segmentation.csv", index=False)

    SegmentationStrainRatio(SegmentedSectionDf_array, colony_widths, colony_heights, output_path)
    SegmentationIntermixing(SegmentedSectionDf, colony_widths, colony_heights, output_path)

    return SegmentedSectionDf





# Creates dataframes from the text image files and applies the necessary functions on them
def ColonySectionEval(mask, channel1, channel2, output_path, file_name):
    # Apply ColonyMask onto the the two channel dataframes
    channel1_array = channel1.to_numpy()
    mask_array = mask.to_numpy()
    channel1_array = np.multiply(channel1_array, mask_array)

    masked_channel1 = pd.DataFrame(channel1_array, columns=channel1.columns)

    channel2_array = channel2.to_numpy()
    channel2_array = np.multiply(channel2_array, mask_array)

    masked_channel2= pd.DataFrame(channel2_array, columns=channel1.columns)

    # Get the colony dimensions for the current section image data
    colony_widths, colony_heights, max_colony_width, max_colony_height = ColonyDimensions(mask)
   
    # Get the segmentation data for the current section image data
    SegmentedSectionDf = SectionSegmentation(masked_channel1, masked_channel2, colony_widths, colony_heights , output_path, file_name)
    # Save it as a csv file
    # SegmentedSectionDf.to_csv(output_path + "/SegmentedSection.csv", index=False)

    return (max_colony_width, max_colony_height)




# Iterates over a whole strain directory
def BatchColonySectionPlot():
    tkinter.Tk().withdraw()

    section_dimensions = pd.DataFrame(columns=['X', 'Width', 'Height'])
    # Average colony width and height at day 1
    section_dimensions.loc[len(section_dimensions.index)] = ['1', 3237.783105, 27.05845042]

    directory_path = tkinter.filedialog.askdirectory(title="Select the input strain directory")
    split_directory_path = directory_path.split('/')
    strain_name = split_directory_path[-1]

    subdirectory_list = os.listdir(directory_path)
    # Iterating over all the subdirectories per strain
    for subdirectory in subdirectory_list:
        subdirectory_path = directory_path + "/" + subdirectory
        # If it is a subdirectory
        if os.path.isdir(subdirectory_path):
            # Iterate over all the directories in the current subdirectory
            for dir in os.listdir(subdirectory_path):
                if os.path.isdir(subdirectory_path + "/" + dir):

                    ColonyMask = pd.DataFrame()
                    GFPChannelDf = pd.DataFrame()
                    mCherryChannelDf = pd.DataFrame()

                    print("Processing: ", dir)
                    file_paths = os.listdir(subdirectory_path + "/" + dir)

                    # Load in all the necessary data frames 
                    for current_file in file_paths:
                        if current_file.find(".txt") != -1:
                            channel = int(current_file[current_file.find("=")+1])
                            if current_file.find("Mask") != -1:
                                ColonyMask = pd.read_csv(subdirectory_path + "/" + dir + "/" + current_file, sep = "\t", header = None)
                                ColonyMask.replace(255, 1, inplace=True)       
                            if channel != 2:
                                if current_file.find("Intensity") != -1:
                                    GFPChannelDf = pd.read_csv(subdirectory_path + "/" + dir + "/" + current_file, sep = "\t", header = None)
                            else:
                                if current_file.find("Intensity") != -1:
                                    mCherryChannelDf = pd.read_csv(subdirectory_path + "/" + dir + "/" + current_file, sep = "\t", header = None)

                    # Calling the ColonySectionEval function on the current section data
                    colony_width, colony_height = ColonySectionEval(ColonyMask, GFPChannelDf, mCherryChannelDf, subdirectory_path + "/" + dir, dir)

                    # Add the current sections dimension into a data frame
                    section_dimensions.loc[len(section_dimensions.index)] = [subdirectory[-1], colony_width, colony_height]

    # Save the colony dimensions as a csv
    section_dimensions.to_csv(directory_path + "/" + strain_name + "SectionDimensions.csv", index=False)





def SectionDimensionPlot():
    tkinter.Tk().withdraw()

    file_path = tkinter.filedialog.askopenfilename(title="Select the input colony dimensions csv")

    split_file_path = file_path.split('/')
    strain_name = split_file_path[-2]
    section_dimensions = pd.read_csv(file_path)
     # Calculating the mean dimensions per colony age
    mean_section_dimensions = section_dimensions.groupby('X', as_index=False)[['Width','Height']].mean()
    
    # Converting the dimensions from pixels to microns
    mean_section_dimensions['Width'] = mean_section_dimensions['Width']*pixel_to_micron
    mean_section_dimensions['Height'] = mean_section_dimensions['Height']*pixel_to_micron
    
    # PLotting the dimensions
    plot = sns.lineplot(x='X', y='value', hue='variable', legend='full', data=pd.melt(mean_section_dimensions, id_vars=['X'], value_vars=['Width', 'Height']))
    plot.set(xlabel='Days', ylabel='Microns', title='Dimensions of ' + strain_name + ' colony sections')
    plot.legend(title='')
    plt.grid()

    # Saving the plot as a PNG
    directory_path = file_path[:file_path.rfind("/")]

    fig = plot.get_figure()
    fig.savefig(directory_path + "/" + strain_name + "SectionDimensions.png")
    plt.close()

BatchColonySectionPlot()
SectionDimensionPlot() # Can only be run after BatchColonySectionPlot()