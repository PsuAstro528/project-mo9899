# Astro 528 Project: Create models for Absorption Lines using Solar Spectrum Data from NEID instrument
## Mary Ogborn and Karthik Yadavalli

### Project Overview:
## A (Short) Scientific Interlude:
The NEID instrument takes spectral observations of the Sun (our sun, Sol) in the wavelength band 4550-4610 Angstroms. All of the absorption lines in this wavelength band are observed and the specific flux over 
this wavelength band is observed by NEID. NEID is a disperser, so the measured NEID fluxes follow a 'blaze' function background. This background blaze needs to be removed from the observed
fluxes before any fitting to absorption lines is done. 

## The Project:
The goal of this project is to find and fit to many absorption lines in this wavelength band. It does so by fitting a high-order Gauss-Hermite polynomial to the shape for each absorption line in the flux observations. 
In our final implementation, we fit to 53 known absorption lines in parallel. To find the "best" fit, the code uses a value we call 'loss', equivalent to the sum of the difference between observed flux and fitted flux over a window 
of wavelengths that includes the wavelength at which the code tries to fit. The code tries a series of wavelength centered on the wavelength manually fed in by the coder and returns the wavelength and evaluated loss of the wavelength
that gives the smallest value of loss. The code also returns the coefficients of the G-H polynomial fit for each predefined absorption line.

There are thus two windows-the window over which loss is calculated for a given wavelength fit, and the window over which different wavlengths fits are tried. We simply tried different window sizes for both and went with the choice of 
window size that produces the "best" fit (here, "best" fit is determined by eye rather than by evaluated loss or some other fit evaluation function). As explained above, before all of this fitting to G-H polynomials and finding absorption 
lines, the background blaze function needs to be removed. Fitting to a blaze function, removing that blaze background, fitting G-H coefficients, and outputting losses and G-H coefficients is all done by the 'Script (Main).jl' script.
In the Main Outputs directory, we plot the original NEID data with a fitted Blaze background, the NEID data with blaze removed, and the absorption lines that have been fitted by the code.

### Installation Instructions
### Testing:
For testing, we use the fitting algorithm to fit to two artifically created lines:
* Perfect Line - A line that is defined by a Gauss-Hermite fit is fed into the algorithm (which tries to fit a Gauss-Hermite polynomial), so the fit should be "perfect". 
* Terrible Line - A horizontal line is fed into the algorithm, so the first should be "terrible".

For both, the serial implementation and the parallel implementation is run. The test should be that both serial and parallel should produce the same output. In addition, the evaluated loss
should be very small for the "perfect" fit. This is all done in the 'Script (Testing).jl' script. This script also creates a digital plot of the perfect fits and terrible fits, and saves those
to the local directory 'Testing_plots', showing how good the "perfect" fits and how bad the "terrible" fits are. The results of testing is done 


### Benchmarking:
For benchmarking, we benchmark the serial and parallel implementation of fitting to the 53 lines. With four workers, we see a speed-up of roughly 2-2.5x each benchmarking. This is done on 
'Script (Timing).jl', and the benchmarking is shown in terminal. 

### Overview of Package Structure:
The main "meat" of this code is in the 'Script (Main).jl' script. The other scripts do testing and benchmarking, as outlined above. The plots from the testing script are output in the 'Testing_plots' directory.
Important codes and functions are in src. The