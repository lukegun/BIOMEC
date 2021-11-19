[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4983)


# BIOMEC
Bayesian Inference and Optimisation for the Monash electrochemical Simulator (MECSim) is the application developed by the
monash electrochemistry group with the assistance of Assosiate Professor Jie Zhang, Emeratious Professor Alan Bond and technical assistance from Dr Gareth Kennedy And Dr Martin Robinson.

It is an automatic plaform for parameterisation that uses mathematical optimisation and Bayesian Inference to calculate parameters involved in the electrochemical simulation.
Built around [MECSim](http://www.garethkennedy.net/MECSim.html) and first applied in the paper titled [A Comparison of Bayesian Inference Strategies for Parameterisation of Large Amplitude AC Voltammetry Derived from Total Current and Fourier Transformed Versions](https://chemistry-europe.onlinelibrary.wiley.com/doi/abs/10.1002/celc.202100391). BIOMEC allows for automated parameterisation of DC and FTAC voltammetry, allowing highly dimensional fits of the posterior distribution.

For an in depth tutorial on installation, application and analysis of BIOMEC and its outputs watch the four part series on [Youtube](https://www.youtube.com/watch?v=LjVesAtftog&list=PLqz7aW7nxQkpNyWkI8JXK4NhhEOlgs2h-).

BIOMEC uses [PINTS](https://github.com/pints-team/pints) for univariant Bayesian inference.

For information of current uses see the original [BIOMEC paper](https://chemistry-europe.onlinelibrary.wiley.com/doi/abs/10.1002/celc.202100391), the [original Bayesian inference paper](https://doi.org/10.1002/celc.201700678) for AC voltammetry or our most recent [featured article](https://doi.org/10.1039/D0CC07549C).

## Installing BIOMEC image
The code is run in a singularity container which works for Ubuntu/UNIX and MAC (untested) OS systems.
Singularity will need to be installed to use the image. Where the guide is seen in the following [website](https://sylabs.io/guides/3.6/user-guide/quick_start.html) or downloaded from connected singularity hub.

Once singularity has been installed, download the BIOMEC file and run the code to create the BIOMEC container (which should be around 580MB). 

```
$ sudo singularity build BIOMEC.simg Singularity.def
```
Once the image is built the imput file (input.txt) can be passed to the image by using the following command.
```
$ ./BIOMEC.simg input.txt
```
or 
```
$ singularity run BIOMEC.simg input.txt
```

This will generate and ouput file with plots and results once completed.

## Generating input files
inputwritter.py can guide users unfamilaur with generating input files to create an input file for the BIOMEC container, this program is contained in the BIOMEC_inputwritter.
Simply run the file using the following command and follow the prompts and an input file will be generated.
```
$ python3 inputwritter.py
```
The output of this file will then be of the form <input.txt> though other names will work.
It is important that a copy of the MECSim Master.inp file is present in the folder you run inputwritter.py as the MECSim input file is required for BIOMEC to run.

Once comfortable with writting the input file it is recommended to use any text editor. 


## Currently Supported Optimizable Parameters
These are the currently supported parameters that can be treated as varibles in BIOMEC for CMA-ES and Bayesian Inference calculations.

| Parameter  | Code # |
| ------------- | ------------- |
| Uncompensated Resistance  | 11  |
| Kinematic Viscosity  | 12  |
| Experimental noise (TCDS only) | 13 |
| Concentration | 21 |
| Diffussion coefficent | 22 |
| Forward Reaction rate | 31 |
| Backward Reaction rate | 32 |
| Formal Potential | 33 |
| Electron Transfer Rate | 34 |
| α or λ | 35 |
| Equilibrium Constant magnitude (func) | 41 |
| Equilibrium Constant θ (func) | 42 |
| Capacitance constants C<sub>0</sub>-C<sub>4</sub> | 51-55 |
| Capacitance Scalar Multiple | 56 |

For parameters that occur on repeatable lines in the MECSim input file, change the number after the parameter code from 1 to what repeat line the parameter of interest is on.

## Running BIOMEC
For an in depth tutorial on installation, application and analysis of BIOMEC and its outputs watch the four part series on [Youtube](https://www.youtube.com/watch?v=LjVesAtftog&list=PLqz7aW7nxQkpNyWkI8JXK4NhhEOlgs2h-).

PDF tutorial may be written up later.


## Supporting Code
 - BIOMEC_inputwritter: Basic terminal/ .exe code for guiding uses in writting the input files for BIOMEC
 - MCMC PLOTTER: code to plot the mcmc output chains from the Iter_log.txt to images 

## Known Issues
 - Custom waveforms have not been tested and Estart and End cannot equal zero.
 - Number of data poins in experimental data must be a multiple of two.

## Citing
Please, cite the original [BIOMEC paper](https://chemistry-europe.onlinelibrary.wiley.com/doi/abs/10.1002/celc.202100391) if you have used this package in a publication.

## License

BIOMEC analysis/python code is open source under the GPL-3.0 License, with MECSim developed by Gareth Kennedy and contaned in the mecsim.cpython-37m-x86_64-linux-gnu.so shared object is under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

## Get in touch
For Questions/Bugs Email me at luke.gundry1@monash.edu.
