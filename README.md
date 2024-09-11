# About

This project has two main purposes
* Analyze Laser Pulses exibiting Random Lasing and Sponaneous emission to extra their Order Parameters
* Use that statisical profile to Generate Synthetic Pulses

# Algorithm

  Order Paramets are extracted using the algorithm outlines in [This Paper](https://onlinelibrary.wiley.com/doi/full/10.1002/lpor.202200314)
  Synthetic data is then made using a veriaty of methods:
  * Distribution Analysis: The PDF on each wavelength is replicated and then used to generate new shots
    - This method does a couple other cool things like restricting varience between wavelength and expanding the PDF parameters over random lasing wavelengths
 * Frequency Mixing: This method simply mixes the amplitudes of differant frequencies between pulses
 * Average + Error: This Method just adds some white gaussian noise ontop of each pulse

THe quality of each algorithm has been evaluated using the MSE between real and synthetic Order parameter histograms.

The code for theese methods is primarily in the **ControlledRunning** folder
 # Exclutions
 
This repository does not include the following due to space restrictions

* The raw data
* Much of the output figures
  
