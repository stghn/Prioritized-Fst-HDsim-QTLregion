#!/bin/bash

# Make alias command work in bash script or bashrc file
shopt -s expand_aliases

QMSim2 -o <<AA > HD_2000QTL_h2040_param_5rep.log
HD_2000QTL_h2040_param_5rep
AA


