#@ File(label = "Input directory", style = "directory") input
#@ File(label = "Output directory", style = "directory") output
#@ String(label = "File suffix (.tif/.nd2)", value = ".tif") suffix

// Replace "\" characters with "/"
input = replace(input, "\\", "/");
output = replace(output, "\\", "/");
dir_list = getFileList(input)

for(i = 0; i < dir_list.length; i++) {
    input_folder = input + "/" + dir_list[i];
	output_folder = output + "/" + dir_list[i];
	
	// Scan folders/subfolders/files to find files with correct suffix
	file_list = getFileList(input_folder);

	for (j= 0; j < file_list.length; j++) {
			if(endsWith(file_list[j], suffix))
				processFile(input_folder, output_folder, file_list[j]);
	}
}

function processFile(input_folder, output_folder, file_name)  {
	file_path = input_folder + file_name;
	split_filename = split(file_name, ".");
	output_subfolder = split_filename[0];
	channel_idx = split(split_filename[1], "-");
	open(file_path);
	
	// Create output folder for current image
	if (!File.exists(output_folder + output_subfolder)){
		File.makeDirectory(output_folder + output_subfolder);
		}
	
	width = getWidth();
	height = getHeight();
	
	/*run("Enhance Contrast", "saturated=0.35");
	run("Apply LUT");*/
  	saveAs("Text Image", output_folder + output_subfolder + "/" + split_filename[0] + "_" + channel_idx[1] + "_Intensity.txt");
  	saveAs("PNG", output_folder + output_subfolder + "/" + split_filename[0] + "_" + channel_idx[1] + ".png");
  	
  	channel_idx = indexOf(file_name, "C=", 0) + 2;
  	channel = parseInt(substring(file_name, channel_idx, channel_idx+1));
	if (channel % 2 != 0){
		img1 = getTitle();
		other_channel = split(file_name, "=");
		open(input_folder + "/" + other_channel[0] + "=2" + substring(other_channel[1], 1, lengthOf(other_channel[1])));
		img2 = getTitle();
		run("Merge Channels...", "c2=[" + img1 + "] c6=[" + img2 + "] create keep ignore");
		saveAs("PNG", output_folder + output_subfolder + "/" + split_filename[0] + "_Composite.png");
		imageCalculator("Add create", img1, img2);
		run("Gaussian Blur...", "sigma=30 scaled");
		setAutoThreshold("Triangle dark");

		run("Convert to Mask");
		run("Options...", "iterations=1 count=1 do=[Fill Holes]");
		run("Options...", "iterations=20 count=1 do=Erode");
		run("Options...", "iterations=20 count=1 do=Dilate");
		saveAs("Text Image", output_folder + output_subfolder + "/" + split_filename[0] + "_Mask.txt");
		saveAs("PNG", output_folder + output_subfolder + "/" + split_filename[0] + "_Mask.png");
		/*run("Options...", "iterations=2 count=1 do=Outline");
		saveAs("Text Image", output_folder + output_subfolder + "/" + split_filename[0] + "_ColonyOutline.txt");*/
	}
	close("*");
}