This folder contains codes and data used to study the morphodynamics 
of bacterial flocks during the colonisation of plant roots.

RProfile_.java 
is an ImageJ plugin used to extract pixel intensity data along a root.
Usage:
* Use segmented line to outline the filament of interest
* Use compile and run function from imageJ and executre RProfile_.java
* Save the result of the analysis as .txt file

Parameter to adjust in RProfile_.java:
* int w_profile = 200		// width of the profile in pixels
* int l_profile = 10		// length of segments in pixels (after resampling)	


* Copy the result in DATA folder for further analysis:
	All filament traced for on given sample are recorded within one 
	folder within data
* Record the metadata of all the filament traced in 
* Use metadata.csv to indicate which timepoint the filament is recorded


RProfile_process_timelapse.py
Series of python scripts to extract length and diameter of filaments
* Copy the result in DATA folder for further analysis:
	All filament traced for on given sample are recorded within one 
	folder within data
* Use metadata.csv in DATA folder to associate timepoints with filament 
	data

Usage:
Run simply the program on a computer with a suitable installation of 
python and associated dependencies.
The program generates information of location of base and tips of filamentous 
as Array files (.p) in the PLOT folder. 

plot_filamentous_flocks.py
Python script for ploting length, number and bacterial density of 
filamentous flock data

Usage:
Run simply the program on a computer with a suitable installation of 
python and associated dependencies.