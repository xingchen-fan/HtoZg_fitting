# HtoZg_fitting

# Environment
The code is tested unfer `CMSSW_12_6_0` with `ROOT 6.24` and `Python 3.9.14`. 

Make sure `RooGaussStepBernstein` is properly installed! In my case, it is installed in the `combine`, specificly in `CMSSW_PATH/src/HiggsAnalysis/CombinedLimit/src/HZGRooPdfs.cxx`.

# User Guide
The functions are defined in `Utilities/bkg_functions_class.py`. A typical definition of a model is like this: 
```
model = #PDF#Class(x, gauss_mu, cat, sigma_init, step_init, ...#PDF# specific arguments...)
```
where `x` is the fitting axis, `gauss_mu` is the mean of the convolution Gaussian, `cat` is the category name, `sigma_init` is the initial value of Gaussian sigma, and `step_init` is the step position initial value. Some default values are assigned already, but feel free to play around with them.

Then define a test statistics (either `RooChi2Var` or `RooNLLVar`) using a `RooDataHist` and `model.pdf`. A typical minimization of certain test statistics is like this:
```
Minimizer_#STAT#(STAT, printLevel, eps, offset, strategy)
```

## F Test
In the `F_test.py`, please change the sample directory. F test of multiple categories can happen in one run, and the print out will become messy. So please modify the file so that you do one category at a time. 

One single F test is defined as
```
singleBernFTest(x, gauss_mu, histogram, cat, method, e_type, eps, offset, strategy)
```
where `x` is the fitting axis, `gauss_mu` is the mean of the convolution Gaussian, `histogram` is the `RooDataHist` to be fit, `cat` is the category name, `e_type` is the error type, `eps` is the EDM, `offset` controls whether using the offset, and `strategy` should be 0.

When calling the F test, please specify the fitting method (Chi2 or NLL).
```
bash$ python3 F_test.py Chi2
```
Output has both the test statistics and the P-values.
