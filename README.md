# Emergence of Urban Growth Patterns from Human Movements

## Empirical City Data

We use three empirical datasets collected from cities in USA, GB and Berlin region. The USA and GB datasets are originally released by US Census Bureau (https://www.census.gov/) and Statistical Office of the European Union (https://ec.europa.eu/info/departments/eurostat-european-statistics_en), respectively. We use the preprocessed datasets extracted from USA 2000 Census and GB 1991 Cenus, which are released in previous research [1]. The Berlin data is extracted from the settlement area distributions in Berlin region in 1910, 1920 and 1945, which are originally acquired from Wikipedia (http://en.wikipedia.org/wiki/Boroughs_and_localities_of_Berlin) and available in the Berlin folder.
[1] Rozenfeld H D, Rybski D, Andrade J S, et al. Laws of population growth[J]. Proceedings of the National Academy of Sciences, 2008, 105(48): 18702-18707.

## Collective Mobility Model
The main.cpp is the source code of a c++ implementation. "main" function specifies the model paramenters and initiates the simulation. "run" and "update" functions maintain the information throughout the simulatoin. Class "t_user" is the implementation of simulated agents, while classes "SortedArray" and "AliasTable" implement the improved sampling algorithms described in SI.

## Enviornment
The code is run and tested on Visual Studio 2012 platform. No external dependency or installation is needed. 

## How to run

1) Please open a new project in Visual Studio 2012 platform, and import main.cpp as the source file.
2) Specificy the model parameters of N, alpha, rho0 in main function.
3) Run the code with "debug" or "release" function in Visual Studio platform. No input data is needed. 
4) The output will be stored in an automaticly created folder in current working path. "p[#epochs].dat" is the simulated urban population distribution.

The typical running time is 12 hours on a workstation with 40 cores of 3.6GHz Intel i7 processor for 20,000 epochs simulation.
