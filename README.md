# Emergence of Urban Growth Patterns from Human Movements

## Empirical City Data

We use three empirical datasets collected from cities in USA, GB and Berlin region. The USA and GB datasets are originally released by US Census Bureau (https://www.census.gov/) and Statistical Office of the European Union (https://ec.europa.eu/info/departments/eurostat-european-statistics_en), respectively. We use the preprocessed datasets extracted from USA 2000 Census and GB 1991 Cenus, which are available available in (https://hmakse.ccny.cuny.edu/software-and-data/). The Berlin data is extracted from the settlement area distributions in Berlin region in 1910, 1920 and 1945, which are originally acquired from Wikipedia (http://en.wikipedia.org/wiki/Boroughs_and_localities_of_Berlin) and available in the Berlin folder.

## Collective Mobility Model
The main.cpp is a c++ implementation source code. "main" function specifies the model paramenters and initiates the simulation. "run" and "update" functions maintain the information throughout the simulatoin. "t_user" is the implementation of each simulated agents. Classes "SortedArray" and "AliasTable" implement the improved sampling algorithms described in SI of the manuscript.

## Enviornment
The code is run and tested on Visual Studio 2012 platform. No external dependency or installation is needed. 

## How to run
The code can run with "debug" and "release" features in Visual Studio platform. No input data is needed. The typical running time is 12 hours on a workstation with 40 cores of 3.6GHz Intel i7 processor for 20,000 epochs simulation.
