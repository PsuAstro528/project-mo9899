# Astro 528 Project: Create models for Absorption Lines using Solar Spectrum Data from NEID instrument
## Mary Ogborn and Karthik Yadavalli

### Goals:
* Create models for known solar absorption lines utilizing Gauss-Hermite polynomials
* Evaluate common shifts in absorption lines and determine what type of broadening or shift is at work utilizing the coefficients found of the Gauss-Hermite Polynomials
	* The first order G-H coefficient finds the depth of the line, the second order coefficient finds the asymmetry of the line, the third order finds the broadening of the line, and so on. This code will find all such coefficients to order N.
	* (So far the code has found G-H polynomials, but has not interpreted them yet)

### Structure of the Code:
* Read in a fits file of data obtained from [here](https://neid.ipac.caltech.edu/search_solar.php) on March-5-2021
* Tested that data exists with matching length values
* Exclude certain pixels near the edge of the detector from the dataset, where the observations are distorted or unreliable
* Remove blaze function background from the observed data, since we are dealing with diffraction data which follows a blaze pattern
* Divide our observation by our blaze function to view the continuum and all the absorption features
* Pick a range of wavelengths to start on the code (this has been extended to the full range available)
* There are three scripts-main, testing, and timing; main fits to a set of 53 lines in parallel and outputs the Gauss-Hermite polynomials from the fit
* Testing creates fake lines and fits to them in both serial and parallel, showing that both produce the same output. 
* Timing re-reads in all the data, and fits to all 53 lines, timing how long such a fit takes. We regularly see a speed up of ~2-2.5x from serial to parallel



### Next Steps:
* Simply need to create a presentation version of this code for our presentation.
* Will need to properly document everything, perhaps in a pdf file, to fully highlight everything that the code does, and will need to package the code, as in lab 9.



[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-f059dc9a6f8d3a56e377f745f24479a46679e63a5d9fe6f495e02850cd0d8118.svg)](https://classroom.github.com/online_ide?assignment_repo_id=5589173&assignment_repo_type=AssignmentRepo)
# Astro 528 [Class Project](https://psuastro528.github.io/project/)

## Project Goals:  
- Put software development skills learned during class and lab exercises into practice on a real world problem
- Gain experience optimizing real scientific code
- Gain experience parallelizing a real scientific code 
- Investigate how performance scales with problem size and number of processors
- Share lessons learned with the class

Prior to the [peer code review](https://psuastro528.github.io/project/code_reviews/), update this readme, so as to make it easy for your peer code reviewer to be helpful.  What should they know before starting the review?  Where should they look first?  

Remember to commit often and push your repository to GitHub prior to each project deadline.

## Class Project Schedule
- Project proposal (due Sept 13)
- Serial version of code (due Oct 4)
- Peer code review (due Oct 11)
- Parallel version of code (multi-core) (due Nov 1)
- Second parallel version of code (distributed-memory/GPU/cloud) (due Nov 18)
- Completed code, documentation, tests, packaging (optional) & reflection (due Dec 2)
- Class presentations (Nov 29 - Dec 9, [detailed schedule](https://github.com/PsuAstro528/PresentationsSchedule2021) )

