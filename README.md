# Resistance to CO<sub>2</sub> diffusion
## Overview
The Mesophyll Conductance project focuses on studying the rate at which carbon dioxide (CO₂) moves from the intercellular air spaces in the leaf into the mesophyll cells. This process plays a critical role in photosynthesis, as it directly impacts CO₂ availability for the enzyme RuBisCO in the chloroplasts. The research aims to better understand the factors affecting mesophyll conductance in plants, which could have implications for improving crop yields and addressing food security challenges.

This repository contains the code, models, and data for the project, which can be used for analyzing and simulating mesophyll conductance in different plant species or under various environmental conditions.

## Features
Data Analysis: Scripts for processing experimental data related to mesophyll conductance measurements.
Mathematical Modeling: Tools and models to simulate mesophyll conductance in plants under varying environmental factors.
Visualizations: Functions to generate graphs and charts for interpreting the results of mesophyll conductance studies.
Parameter Estimation: Code to estimate key physiological parameters affecting mesophyll conductance.
Comparative Analysis: Utilities to compare mesophyll conductance across different plant species or experimental treatments.

The code calculates the resistance to CO<sub>2</sub> diffusion of a leaf based on 
the method of [Berghuijs et al. 2015](10.1016/j.plantsci.2015.06.022)
The code was used in the publication by [Retta et al. 2024](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.20136)

## Tech Stack
The project utilizes the following technologies and libraries:

Python: Main programming language for data analysis and modeling.
NumPy: For numerical computations.
Pandas: For data manipulation and analysis.
Matplotlib/Seaborn: For visualizing results.
SciPy: For scientific computing and modeling.

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

## Contact
If you have any questions or need further information, feel free to reach out:

Author: Moges Retta
Email: [moges.retta@gmail.com]
GitHub: Moges-Retta
