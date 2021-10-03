# Astro 528 Project: Create models for Absorption Lines using Solar Spectrum Data from NEID instrument
## Mary Ogborn and Karthik Yadavalli

### Goals:
-Create models for known solar absorption lines utilizing Gauss-Hermite polynomials
-Evaluate common shifts in absorption lines and determine what type of broadening or shift is at work utilizing the coefficients found of the Gauss-Hermite Polynomials

### Structure of the Code:
-Read in a fits file of data obtained from [here](https://neid.ipac.caltech.edu/search_solar.php)-tested that there is data with matching length values
-Exclude certain pixels from the dataset
-Fit model to a blaze function since we are dealing with diffraction data
-Divide our observation by our blaze function to view the continuum and all the absorption features
-Pick a range of wavelengths to start on the code
-Tested fitting one absorption line using the Gauss-Hermite Polynomials to the fourth order (fit_line_v0), with the Gauss-Hermite polynomial function at the bottom of the code
-Fit many absorption lines (fit_lines_v0)
-Examine the standard deviation for each line model
-Calculated a loss function for each fit-i.e. how the model compares to the actual data for these polynomials
-Print out loss functions and coefficients of the Hermite Polynomials
-tested that the number of lines found is equal to the number of lines we looked for
-tested that we found the expected number of Gauss-Hermite coefficients

### Next Steps:
-Write a function that calculates standard deviation for each line model
-Use the Hermite coefficients to examine the shifts in solar data and attribute to physical processes on the sun


