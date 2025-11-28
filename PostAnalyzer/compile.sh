#!/bin/bash

# compile code (produces two executables)
g++ -o ttbarMakeHist `root-config --cflags --libs` -lMathMore -std=c++11 ttbarMakeHist.cxx
g++ -o ttbarMakePlots `root-config --cflags --libs` -std=c++11 ttbarMakePlots.cxx

# create needed directories if do not exist yet
mkdir -p data mc hist plots
