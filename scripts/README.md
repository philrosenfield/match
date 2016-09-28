# Scripts for MATCH

(These functions are totally seperate from the top level directory)

A bit of what you can do:

## asts.py
Using the match fake file, you may be interested to know the 90% completeness
fraction.

- Load the ast file:
`ast = asts.AST(filename)`
- Compute the completeness
`ast.completeness(interpolate=True, combined_filters=True)`
- Compute the completeness fraction:
`ast.get_completeness_fraction(0.9)`

- Make a completeness plot with the 90% completeness marked:
`ast.completeness_plot(comp_fracs=[0.9])`

- Make the typical magin-magout plots:
`ast.magdiff_plot()`

You could also run it from the command line:
`python -m match.scripts.asts.py -h`


## calcsfh_parallel.py
depreciated

## diagnostics.py
depreciated

## fileio.py
utilities for opening, reading, writing files.

## graphics.py
utilities for visualizing match outputs.

## likelihood.py
stellar_prob is a work in progress, lots of dumb warnings are thrown due to order of operations.

## match_param.py
Make `calcsfh` pararameter files and plot the CMD match will be considering as data.
Broken: Click on a cmd to set mag limits, or use completeness fraction, also can make basic exclude gates.

It's set to read fits tables or match photometery. If it's a fits table, the
filters passed must exist exactly in the fits file.

See example notebooks for some tutorials.

## ssp.py
SSP: Class for visualizations and statistics on `calcsfh` in `ssp` mode. Development in progress.

## utils.py
Common utilities that are not file I/O or graphics.

## cmd.py
CMD: class to read `.cmd` files, an automatic product from `calcsfh`

## sfh.py
SFH: class to read `calcsfh` terminal output, `zcombine` and `zcmerge` output files

## vary_matchparam.py
Write bash or slurm scripts to run a parameter sweep with template calcsfh param, photometry, and fake file. Examples coming... 

## bash
Here are some basic bash scripts to run MATCH.calcsfh

### write_hybridMC_script.sh
Used to run many hybridMC calls on all sfh files in a directory

### theworks.sh
Run all MATCH programs after a calcsfh (with mcdata flag) run
