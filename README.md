# HtoZg_fitting

# Environment
The code is tested under `CMSSW_14_0_4/` and `CMSSW_14_1_0_pre4`, and `combine v10.1.0`. 
Should work without `CMSSW` or `combine` but plain `pyroot`, but require an extra step.

Certain functions are not installed in the default `ROOT`. They are:

* `RooGaussStepBernstein`: Defined in `Utilities/HZGRooPdfs.cxx`
* `RooGaussStepBernsteinRange`: Only used when the range of plotting is different from evaluation. Defined in `Utilities/RooGaussStepBernstein.cxx`.
* `AsymGenGaussian`: Defined in `Utilities/AsymGenGaussian.cxx`.
* `EXModGaus`: Defined in `Utilities/EXModGaus.cxx`.
* `ModGaus`: Defined in `Utilities/ModGaus.cxx`.

Two options to install them.

### Option 1: 
  Install them in the `combine v10.1.0`, specificly in their own `cxx` file at `CMSSW_PATH/src/HiggsAnalysis/CombinedLimit/src/`.

### Option 2: 
  Thier `cxx` and header files locate at `Utilities/`. No dependence on `combine` or `CMSSW`. Generate a c++ dictionary `.so` file, `Utilities/HZGRooPdfs_cxx.so` using the following line ([reference](https://root.cern/manual/io_custom_classes/#generating-dictionaries))
  ```
  root[] .L HZGRooPdfs.cxx+
  ```
Please run it with your version of `ROOT` to overwirte the `.so` file. The `.so` and `.d` files in this repo are generated under `ROOT 6.30`, thus not working under other ROOT versions.

Depdending on the option, you may or may not need to include the library at the beginning of each python script. Please comment them out if you are using option 1.
```
ROOT.gInterpreter.AddIncludePath('../Utilities/HZGRooPdfs.h')
ROOT.gSystem.Load('../Utilities/HZGRooPdfs_cxx.so')
```
> [!NOTE]
> On the LXPLUS machine, the default ROOT doesn't have the FFT package properly installed and when running the fit, you may see some errors accordingly. I suggest you use the ROOT in any CMSSW environment. 

# User Guide

## Sample Reader
Right now, we have two types of files: `.dat` for Cornell MC samples and `.root` for Peking data samples. The readers are defined in `Utilities/sample_reader.py`. Please change the `dir` when using them.

Cornell `.dat` files have specific columns for relevant variables, make sure assign them to the correct `RooRealVar`. Similarly, Peking `.root` files have specific branch names, the asignment is hard-coded in `sample_reader.py`.

The current Cornell samples are available at the [link](https://cernbox.cern.ch/s/EIUbbYia6JCpfC6). `sample_reader.readDat()` requires two directory inputs, first one for the MC and second for the data. If both MC and data samples are stored in the same directory, you still need to provide two same directories.

Update 3/27/25: Run2+3 data samples are read from Rui's dir by `readRuiROOTggFdata` and `readRuiROOTVBFdata`, and signal MC samples are read by `readRuiROOTggFSignalggF`, `readRuiROOTggFSignalVBF`, `readRuiROOTVBFSignalggF` and `readRuiROOTVBFSignalVBF`. They all share the same pattern of usage:
```
readRuiROOTggFdata(x, direct='', year='', bdt1=0, bdt2=0, bdt3=0)
```
Samples are read by the year and each reader object has nCategory of `RooDataHist` for data and nCategory * 2 of `RooDataHist` for signal.

## Background Functions
The functions are defined in `Utilities/bkg_functions_class.py`. A typical definition of a model is like this: 
```
bkg_model = #PDF#Class(x, gauss_mu, cat, sigma_init, step_init, ...#PDF# specific arguments...)
```
where `x` is the fitting axis, `gauss_mu` is the mean of the convolution Gaussian, `cat` is the category name, `sigma_init` is the initial value of Gaussian sigma, and `step_init` is the step position initial value. Some default values are assigned already, but feel free to play around with them.

Note that every bkg function has a sideband version blinding (120, 130) as a member, `self.SBpdf`.

To further automate the process, you can also choose to use the `profileClass` where multiple background functions are defined as members and their essential initial values are assigned through a configuration file.
```
profile = profileClass(x, mu_gauss, cat='', config='')
```

To add your own background model, follow the steps:

1. Add the `.cxx` and `.h` files at `Utilities/`
  
2. Do `.L YourFunc.cxx+` in ROOT, and include the `.so` library in the python script that you want to use the function.

3. Add the class in `Utilities/bkg_functions_class.py` just like the other functions.

4. Add it as a attribute of `profileClass` in `Utilities/profile_class.py`

5. If you want to set the initial values from the config file, you need to edit the config file so that it includes the keys of the parameters.

> [!NOTE]
> This method can **NOT** add your function into `combine`, thus not useable by any `combine` test.
> To do so, replace the first two steps with the [Option 1](https://github.com/xingchen-fan/HtoZg_fitting/edit/main/README.md#option-1).
   
## Signal Functions
The functions are defined in `Utilities/sig_functions_class.py`. A typical definition of a model is like this:
```
sig_model = #PDF#Class(x, MH, cat, ...#PDF# specific arguments...)
```
where `x` is the fitting axis, `MH` is the higgs mass variable (can be a `RooFormulaVar`), and `cat` is the category name.

So far, only Double Sided Crystal Ball (DSCB) function is used. Two versions of DSCB are defined, `DSCB_sys_Class` and `DSCB_Class`.

* `DSCB_sys_Class`:
  ```
  DSCB_sys_Class(x, mu_form, sigmaL_form, sigmaR_form, dMH, sigmaL, sigmaR, cat = "", nL_init = 3, nL_bond = 100, nR_init = 4, nR_bond = 100, alphaL_init = 0.5, alphaR_init = 0.5, di_sigma = False)
  ```
  To include Michael's systematics implementation, `mu_form`, `sigmaL_form` and `sigmaR_form` take any forms of `RooFormulaVar` or a simple `RooRealVar`. Different category/production mode/lepton share the same `MH` variable provided externally in `mu_form`.
  To be able to set `dMH` and `sigmaL` (`sigmaR`) without systematics to certain values and set as constant, there are separate arguements for them.
  And there is a feature to set the parameter values to the post-fit values, `assignVal()`, using a configuration file.
  `assignValModular` set the variables to values provided from the config file. The fitting of the signal model for difference systematic alternatives happens in `Signal_model_preparation/signal_fit_sys.py`. Post-fit parameter values stored in `Config/DSCB_config_sys.json`.

* `DSCB_Class`:
  ```
  DSCB_Class(x, MH, cat = "", sigmaL_init = 1.2, sigmaR_init = 1.2, nL_init = 4, nL_bond = 100, nR_init = 4, nR_bond = 100, alphaL_init = 0.5, alphaR_init = 0.5, di_sigma = False)
  ```
  Different category/production mode/lepton share the same `MH` variable provided by `MH`.
  Each year(era), category, lepton flavor and production mode can have its own signal model. The parameter values stored in the config file are passed to the model by `assignVal`:
  ```
  assignVal(self, config='', year="", cat="", lep="", prod="", sys="nominal", debug = False)
  ```
> [!NOTE]
> For our final analysis, we will only have separate per lepton flavor `lep` and per category `cat` DSCB. Leaving the other arguments blank is expected if the fianl model is a double (el + mu) DSCB per category.

Double DSCB fitting to the nominal MC happens in `Signal_model_preparation/signal_fit_doubleDSCB.py`. Post-fit parameter values are parsered to `Config/config_DSCB_double_flat.json`

To study the effect of having individual signal model per lepton flavor, per year (era), per production mode (ggF, VBF) and per category, class `combineSignal` is defined based on `DSCB_Class` So far, we only include two production modes, ggF and VBF. In total, we define 32 individual DSCB models (Run 2). To use the combined signal model,
```
sig_model = combineSignal(x, MH, cat='', config='')
```
To fit the 32 signal models, go to the `Signal_model_preparation/Multimodel_scripts/` folder and `./run.sh` where the fitting results are save to a log in `/logs` and then parsed to a config file, `Config/config_DSCB_flat.json`.

## Minimizer
Define a test statistics (either `RooChi2Var` or `RooNLLVar`) using a `RooDataHist` and `model.pdf`. A typical minimization of certain test statistics is like this:
```
Minimizer_#STAT#(STAT, printLevel, eps, offset, strategy)
```

## Configuration File
Configurations `.json` files of background functions and signal functions can be found at `Config/`. They contain the post-fit parameter values of the models. The struction is explained as the following.
### Background config file
Each category has its own multi-layer dictionary, for example, for ggf1:
```
"ggf1":{
  "Range":[97, 162],
  "Bins": 260,
  "Chi2":[],
  "SST":[],
  "FT":[],
  "CMSBias":[],
  "bern2":{ ... },
  "bern3":{ ... },
  ...
}
```
The value of `Range` is the mllg fitting range of this category. `Bins` is the number of bins used in the binned-fit. Lists after the key `Chi2`, `SST`, `FT` and `CMSBias` contain the names of the background functions that pass the corresponding test. And the contents of the function keys, such as `bern2`, `bern3`, are the post-fit values of the models. They are used as initial values for the next fit.
### Signal config file
We split the signal model into individual ones depending on the lepton flavor and BDT category (or production mode or year(era)). A typical block in the DSCB config file looks like this:
```
{
   "Htozg_el_cat_ggf1_nominal": {
    "disigma": true,
    "nexp": 0.8886373972111211,
    "sigmaL": 1.65766119136833,
    "sigmaR": 0.823308580027148,
    "nL": 1.993337870785072,
    "nR": 18.37252806042449,
    "alphaL": 1.5050571118870453,
    "alphaR": 1.021482204450084,
    "dMH": 0.0707401603039021
  },
  ...
}

```
`disigma` controls whether there are two sigmas in the DSCB. If false, only `sigmaL` will be used and `sigmaR` is ignored.

## Data Sideband Fit Test
Laurent series is a family that we need to take extra care of, as we are not sure about what powers to use for specific category. So please run `lau2_power_finder.py`, `lau3_power_finder.py` and `lau4_power_finder.py` before running the full test.
```
bash$ ./lau2_power_finder.py -c ggf1 -lo -10 -hi -5 -NLL 1 -d 0 -f 0
```
This `/lau2_power_finder.py` command looks for the best combination of two powers between -10 and -5 for category ggf1. By default, `-NLL 1` indicates using NLL fit, `-d 0` indicates not using di-Gaussian in the convolution, `-f 0` indicates not fixing sigma of the Gaussian in the convolution.
```
bash$ ./lau3_power_finder.py -c ggf1 -lo -9 -hi -6
```
This `/lau3_power_finder.py` command looks for the best combination of three powers, given -9 and -6 are the best two powers, for category ggf1.
```
bash$ ./lau4_power_finder.py -c ggf1 -lo -9 -mi -8 -hi -6
```
This `/lau4_power_finder.py` command looks for the best combination of four powers, given -9, -8 and -6 are the best three powers, for category ggf1.
Be aware that this process takes a LONG time. Values provided in the config file are found through this process, unless the data sample is changed, there is no need to run these again.

To run the sideband NLL fit for category ggf1 with the config file `chi2_config_xgboost_nodrop.json`, run
```
bash$ ./Chi2_test.py -c ggf1 -con ../Config/chi2_config_xgboost_nodrop.json ggf1.log 2>&1
```
After the fit, the post-fit valures are parsered to the config file through `write_config_file`. Fitting results are save into `ggf1.log`. 

Alternative to `write_config_file`, post-fit values can be written into the config file using
```
bash$ ./config_parser.py -c ggf1 -con ../Config/chi2_config_xgboost_nodrop.json -log ggf1.log
```
Add the models that pass the test to the `"Chi2"` value of ggf1 in the config file.
> [!TIP]
> When running `Chi2_test.py` for multiple categories simultaneously, `write_config_file` will cause conflict in the config file, please use the `config_parser.py` method instead!

## Spurious Signal Test
The path of the high stats MC histograms is hard-coded in `Spurious_signal_test_hist.py`. Contact me if you want to do this test and the up-to-date path will be provided. To run the test in ggf1, for example,
```
bash$ ./Spurious_signal_test_hist.py -c ggf1 -conB ../Config/config_xgboost_final_final.json -conS ../Config/config_DSCB_double_flat.json
```
The signal model used is the combination of the two individual DSCB models. Background models tested are selected based on the `"Chi2"` value in the config file.

## F Test
One single F test is defined as
```
singleFTestSidebandNLL(x, pdf_list, histogram, cat, eps, offset, strategy, range, className, sideBand = True, FT = True)
```
where `x` is the fitting axis, `pdf_list` is the list of the functions that you wnat to do F test over, `histogram` is the `RooDataHist` to be fit, `cat` is the category name, `eps` is the EDM, `offset` controls whether using the offset, and `strategy` should be 0, `range` specicies the range to plot, and `className` is the string added to the title of the multi-function plot output file. `sideBand` and `FT` are always set to true.

When calling the F test, NLL fit is used by default. A command to run the F test in ggf1 is
```
bash$ ./F_test.py -c ggf1 -con ../Config/config_xgboost_final_final.json
```
Background models tested are selected based on the `"SST"` value in the config file.Output has both the test statistics and the P-values, as well as a multi-function plot.

## Make Workspaces 
To operate the next test, we need to make toys with `combine`, thus need datacards and worksapces. To generate a background workspace for ggf1, for example,
```
bash$ ./make_bkg_workspace.py -c ggf1 -a 0 -con ../Config/config_xgboost_final_final.json -conS ../Config/config_DSCB_double_flat.json -t FT
```
where `-a` means whether you want to generate and store an Asimov histogram in the workspace. This helps us to cross check with `combine` expected significance calculation.
`-t` indicates for what specific background models you want to create the workspace for. For the bias test, we want to include all the functions passing the F test, thus `-t FT`.
But depending on the bias test (signal injection test) result, new workspaces may be necessary with `-t CMSBias`.

To generate a signal workspace, run
```
bash$ ./make_signal_workspace.py -c ggf1 -conB ../Config/config_xgboost_final_final.json -conS ../Config/config_DSCB_double_flat.json
```
To generate a signal workspace including all the systematics (preferably don't create bkg workspace with this) for 4 ggF and 4 VBF categories, run
```
bash$ ./make_final_combine_workspace.py -s ../Config/DSCB_config_sys.json -b ../Config/config_xgboost_final_final.json -d DATACARD.txt
```
## CMS Bias Test and Signal Injection Test
Two approaches are available for the envelope bias test: using the `combine`, or using `Bias_test/condor/Bias_test_combineToy_condor.py`.
In either approach, we use `combine` to generate toys. So both signal and background worksapces, as well as a datacard, are expected:
```
bash$ combine DATACARD_BIAS_CAT.txt -M GenerateOnly --setParameters pdfindex_CAT=MEMBER_FUNC_ID --toysFrequentist -t 1000 --expectSignal X --saveToys -m 125 --freezeParameters pdfindex_CAT -n .Xsig.CAT.MEMBER_FUNC
```
1000 toys for category `CAT` and bkg function `MEMBER_FUNC` as the true shape with `X` signal injected are generated and saved. Run this command in `Make_combine_workspaces/` to save the toys in that folder.

### Using `Bias_test/condor/Bias_test_combineToy_condor.py`
For example, in ggF1, we want to run tests for 0, 1 and 5 folds of signals for bern4 and pow1 bkg functions:
```
bash$ ./submitter.sh "ggf1" "(bern4 pow1)" "(0 1 5)"
```
where we use LXPLUS condor system to run the jobs. Each job includes only 5 toys.

A folder `Bias_test/condor/MEMBER_FUNC_CAT_Xsig/` (i.e. `Bias_test/condor/pow1_ggf1_0sig/`) will be created to store the corresponding 1000 toy outputs.
And `Bias_test/bias_heatmap.py` and `Bias_test/read_condor_output.c` can get the heat map and pull distributions out of the output files.

To run a single toy(s), you can run 
```
bash$ ./Bias_test_combineToy_condor.py -n 1 -i 1 -f pow1 -c ggf1 -conB ../Config/config_xgboost_final_final.json -conS ../Config/config_DSCB_double_flat.json -s 0
```
This run the bias test on the 6th toy of ggF1 pow1 toys with 0 signal injected. Asumming each condor job have 5 toys, `-n 1` means we will run 1 toy and the toy is the first after the first iteration (5 jobs) `-i 1`.

The steps in the bias test can be summarized as the following: 
  1. Define varables.
  2. Read samples.
  3. Define signal model and fix it.
  4. Define the seed bkg profile functions.
  5. Define the fit bkg profile functions.
  6. Discerete profiling fit.
  7. NLL scan using the envelope.

The fit functions are the RooAbsPdf objects that are used to fit the toy samples. A single discrete profiling fit using the fit functions with the fixed signal model will report the best-fit signal yield. Then, an NLL scan over the signal strength will give an NLL curve for the envelope. The error should be obtained from the curve.

### Using `combine`
Run the `combine` command to get the pull distribution using the toys generated previously:
```
bash$ combine DATACARD_BIAS_CAT.txt -M FitDiagnostics -m 125 --toysFile higgsCombine.Xsig.CAT.MEMBER_FUNC.GenerateOnly.mH125.123456.root -t 1000 --rMin -4 --rMax 4 --cminDefaultMinimizerStrategy=0 -n .Xsig.CAT.MEMBER_FUNC
```
 
> [!TIP]
> An output file `fitDiagnostics.Xsig.CAT.MEMBER_FUNC.root` will be craeted.
> Post-fit signal strengths of 1000 toys are stored in the branch `r` of the tree `tree_fit_sb`. Only plot the ones with `fit_status > 0`. ([Reference](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part3/nonstandard/#roomultipdf-conventional-bias-studies))

## Plotting Features
In the `Multi_plot/` folder, there are three important plotting scripts:
* `multi_plot.py`: Plot selected bkg models with data sideband, controlled by the arguments.
* `money_plot.py`: Plot a combined S+B and data of all categories over a given range (toy only before unblinding).
* `plot_signal_errorbar.py`: Plot observed signal strength of all categories with errorbars (toy only before unblinding).
> [!NOTE]
> When `matplotlib` is called, please use a CMSSW environment where the library is correctly installed. You may need to switch to a different environment from where you are running the rest of the scripts.
## Reference Log
I upload some results when I run the scripts so that you can comapre with.

