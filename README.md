# Regress-Bivariate

This MATLAB code produces a weighted least squares fit of a straight line to a set of points with error in both coordinates.

### Features

- It can handle bivariate regression where the errors in both coordinates are correlated 
- It is capable of performing force-fit regression
- For details on the algorithm used in the code, see *Unified Equations for the slope, intercept, and standard errors of the best straight line* ([York et al., 2004](https://aapt.scitation.org/doi/abs/10.1119/1.1632486))

### Reference to Cite

Thirumalai, K., A. Singh, and R. Ramesh (2011), *A MATLAB code to perform weighted linear regression with (correlated or uncorrelated) errors in bivariate data*, Journal of the Geological Society of India, 77(4), 377–380, doi: [10.1007/s12594-011-0044-1](https://link.springer.com/article/10.1007/s12594-011-0044-1) 
[[PDF](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.464.8784&rep=rep1&type=pdf)].

## Input

The code accepts an Excel spreadsheet or MATLAB array with the following format as input:
```
Column 1 | Column 2 | Column 3 | Column 4 | Column 5    |
X data   | sigX     | Y data   | sigY     | r^2(sigmas) | 
```
where, `sigX` and `sigY` are the errors on `X` and `Y` respectively

**Notes:**
- Column 5 is optional 
- Force-fit regression works with `sigX` and/or `sigY` set to `10e-9`.

## Output

- Graph of Y vs. X data (with errorbars)
- Comparison of Weighted Linear Regression vs. Simple Linear Regression
- Errors on slope and intercept of the line of best-fit

## Test Code

In it’s current version the code will run *as is* by using the data in `Pearson.xls` as input.

## Versions

- Regress-Bivariate Version 1.0.0 | May 30, 2018 | (Initial Deployment)

## Author

Written by [Kaustubh Thirumalai](mailto:kaustubh@arizona.edu) at the University of Texas Institute for Geophysics in July, 2009; currently at the Department of Geosciences, University of Arizona.

## Acknowledgments

- Arvind Singh (Physical Research Laboratory)
- R. Ramesh (Physical Research Laboratory)
