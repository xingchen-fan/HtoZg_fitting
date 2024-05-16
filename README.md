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

Note: On the LXPLUS machine, the ROOT doesn't have the FFT package properly installed and when running the fit, you may see some errors accordingly. I suggest you use the ROOT in any CMSSW environment. 

# User Guide

## Sample Reader
Right now, we have two types of files: `.dat` for Cornell MC samples and `.root` for Peking data samples. The readers are defined in `Utilities/sample_reader.py`. Please change the `dir` when using them.

Cornell `.dat` files have specific columns for relevant variables, make sure assign them to the correct `RooRealVar`. Similarly, Peking `.root` files have specific branch names, the asignment is hard-coded in `sample_reader.py`.

The current Cornell samples are available at the [link](https://cernbox.cern.ch/s/EIUbbYia6JCpfC6). `sample_reader.readDat` has two directory inputs, first one for the MC and second for the data. If both MC and data samples are stored in the same directory, you still need to provide two same directories.

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
singleFTestSidebandNLL(x, pdf_list, histogram, cat, eps, offset, strategy, range, className)
```
where `x` is the fitting axis, `pdf_list` is the list of the functions that you wnat to do F test over, `histogram` is the `RooDataHist` to be fit, `cat` is the category name, `eps` is the EDM, `offset` controls whether using the offset, and `strategy` should be 0, `range` specicies the range to plot, and `className` is the string added to the title of the multi-function plot output file.

Note: Sideband fit of FFT evaluated function is not possible, thus the full range fit is conducted regardless of `range` 

When calling the F test, NLL fit is used by default
```
bash$ python3 F_test.py CAT XLOW FUNC
```
where `CAT` is the category name (u1, u2, u3 or u4 for now), `XLOW` is the lower boundary of mllg fitting range and `FUNC` is the family of the functions (bern, pow, exp or lau for now).

Output has both the test statistics and the P-values, as well as a multi-function plot.

## Spurious Signal Test
In the `Spurious_signal_test.py`, please change the sample directory. Chi2 fit with SumW2 and Poisson options are used in this test. By default, the signal model is Double-sided Crystal Ball (DSCB), and the signal yield with the errors are reported. Both the signal fit and S+B fit plots are saved.

## Bias Test
Two approaches are available for the envelope bias test: using the `combine`, or using `Bias_test/Bias_test.py`. 
### Using `Bias_test/Bias_test.py`
```
bash$ python3 Bias_test.py N_toys
```
where you decided how many toy samples generated for each member of the profile.

The steps can be summarized as the following: 1. Define varables. 2. Read samples. 3. Define signal model and fix it. 4. Define the seed bkg profile functions. 5. Define the fit bkg profile functions. 6. Fit the seed to the sample histogram. 7. Gnerate toy samples. 7. Discerete profiling fit. 8. NLL scan using the envelope.

After reading the samples, select which category to do the test in the script. The seed functions are used to generate toy samples, and they are fit to the histogram of the category accoding to the selection. 

The fit functions are the RooAbsPdf objects that are used to fit the toy samples. A single discrete profiling fit using the fit functions with the fixed signal model will report the best-fit signal yield. Then, an NLL scan over the signal strength will give an NLL curve for the envelope. The error should be obtained from the curve.

Note: This script is still experimental. 

### Using `combine`
Please run the script `Make_combine_workspaces/make_workspace_bias_test.py`.

In the script, firstly define the variables, read the samples and select the category you want to test. Similarly, define the signal and background functions. After fitting the signal function and background functions, output the two workspaces (one signal and one background). Write the `combine` datacard.

Run the `combine` command to generate toys for a certain member of the envelope:
```
bash$ combine DATACARD_BIAS_CAT.txt -M GenerateOnly --setParameters pdfindex_CAT=MEMBER_FUNC_ID --toysFrequentist -t 1000 --expectSignal 0 --saveToys -m 125 --freezeParameters pdfindex_CAT -n .CAT.MEMBER_FUNC
```
Next, run the `combine` command to get the pull distribution using the toys generated previously:
```
bash$ combine DATACARD_BIAS_CAT.txt -M FitDiagnostics -m 125 --toysFile higgsCombine.CAT.MEMBER_FUNC.GenerateOnly.mH125.123456.root  -t 1000 --rMin -4 --rMax 4 --cminDefaultMinimizerStrategy=0 -n .CAT.MEMBER_FUNC
```

## Reference Log
I upload some results when I run the scripts so that you can comapre with.

