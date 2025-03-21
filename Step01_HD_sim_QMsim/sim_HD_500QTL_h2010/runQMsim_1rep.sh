#!/bin/bash

# Make alias command work in bash script or bashrc file
shopt -s expand_aliases

QMSim2 -o <<AA > HD_500QTL_h2010_param_1rep.log
HD_500QTL_h2010_param_1rep
AA



