# HtoZg_fitting

# Environment
The code is tested under `CMSSW_12_6_0` with `ROOT 6.24`, `Python 3.9.14` and `combine v9.1.0`. 
Should work without `CMSSW` or `combine` but plain `pyroot`, but require an extra step.

Make sure `RooGaussStepBernstein` is properly installed!

### Option 1: 
  It is installed in the `combine v9.1.0`, specificly in `CMSSW_PATH/src/HiggsAnalysis/CombinedLimit/src/HZGRooPdfs.cxx`.

### Option 2: 
  No dependence on `combine` or `CMSSW`. Generate a c++ dictionary `.so` file, `Utilities/HZGRooPdfs_cxx.so` using the following line ([reference](https://root.cern/manual/io_custom_classes/#generating-dictionaries))
  ```
  root[] .L HZGRooPdfs.cxx+
  ```
Please run it with your version of `ROOT` to overwirte the `.so` file. The `.so` and `.d` files in this repo are generated under `ROOT 6.30`, thus not working under other ROOT versions.

# User Guide

## Sample Reader
Right now, we have two types of files: `.dat` for Cornell MC samples and `.root` for Peking data samples. The readers are defined in `Utilities/sample_reader.py`. Please change the `dir` when using them.

Cornell `.dat` files have specific columns for relevant variables, make sure assign them to the correct `RooRealVar`. Similarly, Peking `.root` files have specific branch names, the asignment is hard-coded in `sample_reader.py`.

## Background Functions
The functions are defined in `Utilities/bkg_functions_class.py`. A typical definition of a model is like this: 
```
bkg_model = #PDF#Class(x, gauss_mu, cat, sigma_init, step_init, ...#PDF# specific arguments...)
```
where `x` is the fitting axis, `gauss_mu` is the mean of the convolution Gaussian, `cat` is the category name, `sigma_init` is the initial value of Gaussian sigma, and `step_init` is the step position initial value. Some default values are assigned already, but feel free to play around with them.

## Signal Functions
The functions are defined in `Utilities/sig_functions_class.py`. A typical definition of a model is like this:
```
sig_model = #PDF#Class(x, MH, cat, ...#PDF# specific arguments...)
```
where `x` is the fitting axis, `MH` is the higgs mass variable, and `cat` is the category name.

## Minimizer
Define a test statistics (either `RooChi2Var` or `RooNLLVar`) using a `RooDataHist` and `model.pdf`. A typical minimization of certain test statistics is like this:
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

## Spurious Signal Test
In the `Spurious_signal_test.py`, please change the sample directory. Chi2 fit with SumW2 and Poisson options are used in this test. By default, the signal model is Double-sided Crystal Ball (DSCB), and the signal yield with the errors are reported. Both the signal fit and S+B fit plots are saved. 

