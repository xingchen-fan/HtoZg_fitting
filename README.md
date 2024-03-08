# HtoZg_fitting

# Environment
The code is tested unfer `CMSSW_12_6_0` with `ROOT 6.24` and `Python 3.9.14`. 

Make sure `RooGaussStepBernstein` is properly installed! In my case, it is installed in the `combine`, specifiacally in `CMSSW_PATH/src/HiggsAnalysis/CombinedLimit/src/HZGRooPdfs.cxx`.

# User Guide
The minimization is highly integrated due to some problem with pyroot. It is defined in `Utilities/bkg_functions_fit.py`. The typical minimization of certain background pdf is like this:
```
#PDF#Minimization(x = fit axis, gauss_mu = convolution gaussian mu, histogram, ...PDF specific parameters...,
                  printLevel = fitter print level, eps = EDM, offSet = offset of test statistics, strategy = please use 0)
```
Some default values are assigned already, but feel free to play around with them.

## F Test
In the `F_test.py`, please change the sample directory. F test of multiple categories can happen in one run, and the print out will become messy. So please modify the file so that you do one category at a time. When calling the F test, please specify the fitting method (Chi2 or NLL).
```
#bash: python3 F_test.py Chi2
```
