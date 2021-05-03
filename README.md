[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4983)


# BIOMEC
Bayesian Inference and Optimisation for the Monash electrochemical Simulator (MECSim) is the application developed by the
monash electrochemistry group with the assistance of Assosiate Professor Jie Zhang, Emeratious Professor Alan Bond and technical assistance from Gareth Kennedy And Martin Robinson.

It is an automatic plaform for parameterisation that uses mathmatical optimisation and Bayesian Inference to calculate parameters involved in the electrochemical simulation.
Built around [MECSim](http://www.garethkennedy.net/MECSim.html) and first applied in the PAPER. BIOMEC allows for automated parameterisation of DC and FTAC voltammetery, allowing highly dimensional fits of the posteriour distrabution.

BIOMEC uses [PINTS](https://github.com/pints-team/pints) for univariant Bayesian inference.

For information of current uses see [original Bayesian inference paper](https://doi.org/10.1002/celc.201700678) or our most recent [featured article](https://doi.org/10.1039/D0CC07549C).

## Installing BIOMEC image
The code is run in a singularity container which works for Unbuntu/UNIX and MAC (untested) OS systems.
Singularity will need to be installed to use the image. Where the guide is seen in the following [website](https://sylabs.io/guides/3.6/user-guide/quick_start.html) or downloaded from connected singularity hub.

Once singularity has been installed, download the BIOMEC file and run the code to create the BIOMEC container (which should be around 580MB). 

```
$ sudo singularity build BIOMEC.simg Singularity.def
```
Once the image is built the imput file (input.txt) can be passed to the image by using the following command.
```
$ ./BIOMEC.simg input.txt
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

## Running BIOMEC
PDF tutorial or youtube videos to come.


## Supporting Code
 - BIOMEC_inputwritter: Basic terminal/ .exe code for guiding uses in writting the input files for BIOMEC
 - MCMC PLOTTER: code to plot the mcmc output chains from the Iter_log.txt to images 

## Known Issues
 - Custom waveforms have not been tested and Estart and End cannot equal zero.
 - Number of data poins in experimental data must be a multiple of two.

## Citing
Please, cite [Soon to be published paper]() if you have used this package in a publication, which will be updated once the paper is published.

## License

BIOMEC analysis/python code is open source under the GPL-3.0 License, with MECSim developed by Gareth Kennedy and contaned in the mecsim.cpython-37m-x86_64-linux-gnu.so shared object is under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

## Get in touch
For Questions/Bugs Email me at luke.gundry1@monash.edu.
