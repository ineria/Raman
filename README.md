# Raman
Data analysis for Raman spectra


Raman_compounds_calibration_v2 extracts two peak values from Raman spectra and creates a ratio fraction I1/(I1+I2). 
  Intended to be used for getting a calibration equation for a two component mixture with varying amounts of component 1
  Plots I1/(I1+I2) vs %of1 and performs linear regression
  Reads in multiple files
  
Raman_ratiomap_data_analysis_v2 used to extract two peak values from Raman spectra and creates a ratio fraction I1/(I1+I2).
  2D map data needed
  averages data along y axis and plots 1D of ratio fractions
  Reads in multiple files
  
Raman_lineplot_with_baseline-removal_pick-peaks reads in txt files of Raman spectra
  performs a baseline correction using asymmetric least squares
  also able to find peaks and write to list
  will perform 2D line plot of spectra
  Reads in multiple files

Raman_plot_single_with_baseline-removal_pick-peaks 
  used to tune parameters for asymmetric least square baseline correction in Raman_lineplot_with_baseline-removal_pick-peaks
  can be used to tune to find peaks and write to list
  READS IN SINGLE FILE ONLY
