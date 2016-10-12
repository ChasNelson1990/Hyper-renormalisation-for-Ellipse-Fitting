# Hyper-renormalisation-for-Ellipse-Fitting
This repository contains a MATLAB implementation of *Hyper-renormalization: Non-minimization approach for geometric estimation* by Kanatani, Al-Sharadqah, Chernov and Sugaya (2014). Available here: http://dx.doi.org/10.2197/ipsjtcva.6.143.

- Implementation by **Carl J. Nelson**

***Our thanks go to Professor Kenichi Kanatani (Okayama University, Japan) and Dr. Yasuyuki Sugaya (Toyohashi University of Technology, Japan) for their assistance in understanding and implementing this method.***

## Overview
This repository provides a MATLAB implementation for fitting ellipses to a set of data points in 2D.

This project contains one MATLAB function:
1. Main Hyper-Renormalisation code, see Kanatani et al for details (hyperRenormalisation.m)

For general usage:
- Open MATLAB, right click project directory and add to path.
- From the command prompt, call `data=hyperRenormalisation(X);` to run on data points `X` (and n-by-2 array) with default options.
- Alternatively, define the parameters, as discussed in the paper and the code preamble.

## License
This implemented code is licensed under the GNU General Public License Version 3.
- For alternative licenses, please contact *carl.nelson@durham.ac.uk*.

The original concepts and ideas remain the property of the authors: Kanatani, Al-Sharadqah, Chernov and Sugaya.
