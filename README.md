# Resistance to CO<sub>2</sub> diffusion

The code calculates the resistance to CO<sub>2</sub> diffusion of a leaf based on 
the method of [Berghuijs et al. 2015](10.1016/j.plantsci.2015.06.022)
The code was used in the publication by [Retta et al. 2024](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.20136)

## app.py
A tkinter based gui for running the application. This app can be run to 
interact with the program using gui. The app allows selection of species, 
treatments, analysis type and the results are output in the consol and excel
files are save in the project folder.

## Main.py
Main program to calculate the resistances to CO2 transport by leaf 
microstructure and carry sensitivity analysis

## Assumed_values.py
The values of all assumed values are listed in this file

## Cell_component.py
A python class to calculate the resistance of a cell component such as 
cell wall, plasma membrance, cytosol.

## Cell.py
A python class that model cell as having surface and total resistance and 
calculate total resistance of components in a cell type, exposed length of 
mesophyll and chloroplast per leaf surface.

## Data.py
A python class that returns the parameters of the model classified according to
anatomical data, biochemical capacity and so on for each replicate plant.

## Electron.py
A python class to calculate the rate of electron transport

## Gas_exchange_measurement.py
A python class to extract gas exchange data for a selected species and growth and 
measurement conditions
The code is explaind in detail on [Github](https://github.com/Moges-Retta/FvCB_Parameters)

## Make_mesophyll.py
Make mesophyll tissue from inputs of anatomical and transport properties
The class allows construction of cell component from the inputs data or from
an object given as an input and the type of tissue the object belongs

## Names.py
An enum for names used in the code

## Photosynthesis.py
A python class to calculate the rate of photosynthesis for a selected species
given the measurement conditions

## Physical_constants.py
Physical constants used in the model

## Plant.py
A python class to calculate plant level traits

