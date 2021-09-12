# Project Proposal
## Mary Ogborn and Karthik Yadavalli

In this project, we will fit a model to each of the absorption lines found in the NEID Solar data. These lines shift and we can measure shifts and shape changes that are common to many lines. (i.e. Doppler effect, width change from Earth's motion, and stellar variability.)

Our inputs are the approximately ~40,000 fits files of L1 Solar Data from the NEID instrument, which can be found at this [site](https://neid.ipac.caltech.edu/search_solar.php). The given data was all taken during the instrument's high resolution mode. 

The output of this project are models that fit over the absorption lines and can thus reveal stellar variability due to common shifts we find in the created models, which can then be used to determine properties such as radial velocity. 

One way to test this would be to inject a simulated absorption line then see how the model fits certain parameters versus what we had originally specified as the parameters for that simulated line. 

We could potentially parallalize over multiple spectra, but after discussion with Dr. Ford, it may be better to parallalize over multiple lines within a spectrum.

In order to effectively parallalize, we will be writing the code in julia as it is suited for high-performance numerical computing. It will also be helpful to write in julia as many of the exercises presented in class will be written in julia. We will begin by looking into the available functions and items related to the project in JuliaAstro in order to utilize already written codes. 