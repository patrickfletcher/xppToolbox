# xppToolbox

A small collection of tools used for

1. converting XPP ode files to other formats
2. plotting XPP bifurcation diagrams

Download this folder and add it to your Matlab path to easily use these tools in your projects.

## 1. Converting XPP ode files

A parser for XPP files, and code generation tools. `parseODEfile` reads an ODE file and collects all the information contained into a structure. The code converter tools make use of this to convert the ODE file definition to the following formats:

* Matlab m-file or mex-file versions of the ODE system. Use `ode2m` or `ode2mex`
* OpenCL version that runs with clODE (link to repo pending). Use `ode2cl`
* XML file compatible with [VFGEN](http://www.warrenweckesser.net/vfgen/). Use `ode2vfgen`
* VCML file compatible with the ODE mode of [Virtual Cell](http://vcell.org/). Use `ode2vcml`

These tools can be used alone, or via the `xppConverter` convenience routine.

A novel syntax can be used in the ODE file to specify a range of values for a parameter or initial condition. Simply append a set of square brackets with the low and high value of the range separated by a comma: `par p1=2.0[1.0,3.0]`. The XPPaut parser ignores this, but note that it seems to append the digits following preceding the comma into the parameter value, which is undesireable. This syntax should be used only with the Matlab parser.

An Euler method solver (`ode_euler`) is also included that allows solving ODEs with weiner variables in Matlab. See `demo_ode2m` for an example.

Also included is `package4XPP`, which easy packaging of parameter and initial condition vectors for writing back to ODE files using `ChangeXPPodeFile` and `ChangeXPPsetFile` from [XPP-Matlab](http://www2.gsu.edu/~matrhc/XPP-Matlab.html) by Rob Clewley

## 2. Plotting XPP bifurcation diagrams and nullclines

Functions `plotxppaut1` and `plotxppaut2` read and plot one- or two-parameter bifurcation diagram data saved using XPPAUT's auto/file/write-points. The line properties can be set independently for each type of solution plotted using a simple structure of parameters, which can be generated using `plotxppautset1` or `plotxppautset2`.

`plotnullclines` reads data saved from XPP's nullcline menu and plots them in Matlab. The x and y nullcline appearance can be controlled using Matlab's line properties.

Inspired by PlotXppaut by Mohammad S. Imtiaz, which can be found on the [XPPaut website](http://www.math.pitt.edu/~bard/xpp/xpp.html)