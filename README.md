# xppToolbox

A small collection of tools used for

1. converting XPP ode files to other formats
2. plotting XPP bifurcation diagrams

Download this folder and add it to your Matlab path to easily use these tools in your projects.

## 1. converting XPP ode files

A parser for XPP files, and code generation tools to create:

* Matlab m-file or mex-file versions of the ODE system
* OpenCL version that runs with clODE (link to repo)
* XML file compatible with [VFGEN](http://www.warrenweckesser.net/vfgen/)
* VCML file compatible with the ODE mode of [Virtual Cell](http://vcell.org/)

Also included are minor updates to tools for writing data back to XPP files, based on XPP-Matlab by Rob Clewley

## 2. plotting XPP bifurcation diagrams and nullclines

Functions plotxppaut1/2 read data output from XPPAUT's auto/file/write-points option and plot them in Matlab. The line properties can be set independently for each type of solution plotted using a simple structure, which can be generated using plotxppautset1/2. 

plotnullclines reads data saved from XPP's nullcline menu and plots them in Matlab. The x and y nullcline appearance can be controlled using Matlab's line properties.

Inspired by PlotXppaut by Mohammad S. Imtiaz
