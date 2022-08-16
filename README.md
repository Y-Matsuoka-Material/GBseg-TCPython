# GBseg-TCPython
A demo code for calculating grain boundary segregation with Hillert's parallel tangent law using TC-Python.

## Installation

This code itself does not require installation. However, it depends on [Thermo-Calc](https://thermocalc.com/) and [TC-Python](https://thermocalc.com/products/software-development-kits/tc-python/), so please install TC-Python first.

## Run

Running calculate.py will calculate and display the grain boundary segregation of Cantor alloy.

By default, the grain interior is set as FCC single-phase, but you can run the calculation considering other equilibrium phases by reversing the comment-outs in lines 13 and 14.

## Modification

Grain boundary segregation for other alloys and conditions can also be calculated by modifying the code.

You can change target alloy, phases in grain interior, database and temperature range of the calculation by editing lines between 7 and 23.

## Reference

If you use this code in your work, please cite the reference bellow.

[TO BE FILLED] 
