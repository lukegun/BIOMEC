[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4983)


# BIOMEC
Bayesian Inference and Optimisation for the Monash electrochemical Simulator (MECSim) is the application developed by the
monash electrochemistry group (SOME LINK) of Assosiate Professor Jie Zhang and Emeratious Professor Alan Bond.

It is an automatic plaform for parameterisation that uses mathmatical optimisation and Bayesian Inference to calculate parameters involved in the electrochemical simulation.
Built around [MECSim](http://www.garethkennedy.net/MECSim.html) and first applied in the PAPER. BIOMEC allows for automated parameterisation of DC and FTAC voltammetery, allowing highly dimensional fits of the 

## Installing BIOMEC image
The code is run in a singularity container which works for Unbuntu/UNIX and maybe MAC OS systems and has been tested in version 3.1.
Singularity will need to be installed to use the image seen where the guide is seen in the following [website](https://sylabs.io/guides/3.6/user-guide/quick_start.html) or downloaded from conectd singularity hub.

Once singularity has been installed download the downloaded BIOMEC file and run the code to create the BIOMEC image

```
$ sudo singularity build BIOMEC.simg Singularity.def
```
Once the image is built the imput file (input.txt) can be passed to the image by using the following command
```
$ ./BIOMEC.simg input.txt
```
this will generate and ouput file with plots and results once completed


## Running BIOMEC
PDF tutorial or youube videos to come.


## Supporting Code
BIOMEC_inputwritter: Basic terminal/ .exe code for guiding uses in writting the input files for BIOMEC
MCMC PLOTTER: code to plot the mcmc output chains from the Iter_log.txt to images 

## Citing
Cite the following paper if you have used this package in a publication PAPER

## Known Issues
Custom waveforms have not been tested and Estart and End cannot equal zero.
Experimental data must be a multiple of two.

## License

BIOMEC analysis/python code is open source under GPL-3.0 License, with MECSim developed by Gareth kennedy and contaned in the mecsim.cpython-37m-x86_64-linux-gnu.so is closed source.

## Get in touch
For Questions/Bugs Email me at luke.gundry1@monash.edu.