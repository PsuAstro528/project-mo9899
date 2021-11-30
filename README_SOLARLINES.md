# Astro 528 Project: Create models for Absorption Lines using Solar Spectrum Data from NEID instrument
## Mary Ogborn and Karthik Yadavalli

### Project Overview:
## A (Short) Scientific Interlude:
The NEID instrument takes spectral observations of the Sun (our sun, Sol) in the wavelength band 4550-4610 Angstroms. All of the absorption lines in this wavelength band are observed and the specific flux over 
this wavelength band is observed by NEID. 
## The Project:
The goal of this project is to find and fit to many absorption lines in this wavelength band. It does so by fitting a high-order Gauss-Hermite polynomial to the shape for each absorption line in the flux observations. 
In our final implementation, we fit to 53 known absorption lines in parallel. To find the "best" fit, the code uses a value we call 'loss', equivalent to the sum of the difference between observed flux and fitted flux over a window 
of wavelengths that includes the wavelength at which the code tries to fit. The code tries a series of wavelength centered on the wavelength manually fed in by the coder and returns the wavelength and evaluated loss of the wavelength
that gives the smallest value of loss. 

There are thus

### Installation Instructions
### Testing:
### Benchmarking:
### Overview of Package Structure:
