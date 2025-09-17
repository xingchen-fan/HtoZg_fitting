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
> This method will not add your function into `combine`, thus not useable by any `combine` test.
   
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
  `assignValModular` set the variables to values provided from the config file.

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

To study the effect of having individual signal model per lepton flavor, per year (era), per production mode (ggF, VBF) and per category, class `combineSignal` is defined based on `DSCB_Class` So far, we only include two production modes, ggF and VBF. In total, we define 32 individual DSCB models (Run 2). To use the combined signal model,
```
sig_model = combineSignal(x, MH, cat='', config='')
```
To fit the signal MC samples, go to the `Signal_model_preparation/` folder and `./run.sh` where the fitting results are save to a log in `/logs` and then parsed to a config file. Right now we only use DSCB model.

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
  "Range":97,
  "Chi2":[],
  "SST":[],
  "FT":[],
  "CMSBias":[],
  "bern2":{ ... },
  "bern3":{ ... },
  ...
}
```
The value of `Range` is the starting mllg value of the fitting range of this category. Lists after the key `Chi2`, `SST`, `FT` and `CMSBias` contain the names of the background functions that pass the corresponding test. And the contents of the function keys, such as `bern2`, `bern3`, are the post-fit values of the models. They are used as initial values for the next fit.
### Signal config file
We split the signal model into individual ones depending on the year(era), lepton flavor and production mode. A typical block in the DSCB config file looks like this:
```
{
    "ggf1": {
        "2016": {
            "el": {
                "disigma ggf": 1,
                "MH ggf": 124.99,
                "nexp ggf": 0.31,
                "sigmaL ggf": 0.84726,
                "sigmaR ggf": 0.39915,
                "nL ggf": 5.327,
                "nR ggf": 6.3908,
                "alphaL ggf": 0.55437,
                "alphaR ggf": 0.53537,
              "mu": {
              "mu":{
                "disigma vbf": 1,
                "MH vbf": 124.99,
                "nexp vbf": 0.05,
                "sigmaL vbf": 0.8124,
                "sigmaR vbf": 0.43476,
                "nL vbf": 5.5343,
                "nR vbf": 6.6618,
                "alphaL vbf": 0.5494,
                "alphaR vbf": 0.55964
            },
...
        }
    }
}

```
Notise that each key has the production mode name in the end, `ggf` or `vbf`. Whether two sigmas are used in the model is controled by `disigma` and `nexp` shows the expected signal events in this category, year, lepton channel and production mode. The rest of parameter values is obvious as they are the post-fit values.
## Chi^2 Test (Data Sideband Fit Test)
Laurent series is a family that we need to take extra care of, as we are not sure about what powers to use for specific category. So please run `lau2_power_finder.py`, `lau3_power_finder.py` and `lau4_power_finder.py` before running the full test.
```
bash$ ./lau2_power_finder.py -c ggf1 -lo -10 -hi -5
```
This `/lau2_power_finder.py` command looks for the best combination of two powers between -10 and -5 for category ggf1.
```
bash$ ./lau3_power_finder.py -c ggf1 -lo -9 -hi -6
```
This `/lau3_power_finder.py` command looks for the best combination of three powers, given -9 and -6 are the best two powers, for category ggf1.
```
bash$ ./lau4_power_finder.py -c ggf1 -lo -9 -mi -8 -hi -6
```
This `/lau4_power_finder.py` command looks for the best combination of four powers, given -9, -8 and -6 are the best three powers, for category ggf1.
Be aware that this process takes a LONG time. Values provided in the config file are found through this process, unless the data sample is changed, there is no need to run these again.

To run the sideband Chi^2 fit for category ggf1 with the config file `chi2_config_xgboost_nodrop.json`, run
```
bash$ ./Chi2_test.py -c ggf1 -con ../Config/chi2_config_xgboost_nodrop.json ggf1.log 2>&1
```
After the fit, the post-fit valures are parsered to the config file through `write_config_file`. Fitting results are save into `ggf1.log`. 

Alternative to `write_config_file`, post-fit values can be written into the config file using
```
bash$ ./config_parser.py -c ggf1 -con ../Config/chi2_config_xgboost_nodrop.json -log ggf1.log
```
Add the models that pass the test to the `"Chi2"` value of ggf1 in the config file.

## Spurious Signal Test
To make it easier for people to run this test, the simulations with extended DY are saved as histograms at `Data/`. To run the test in ggf1, for example,
```
bash$ ./Spurious_signal_test_hist.py -c ggf1 -conB ../Config/chi2_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
```
The signal model used is the combination of the individual models described previously by defining a `combineSignal` object. Background models tested are selected based on the `"Chi2"` value in the config file.

## F Test
One single F test is defined as
```
singleFTestSidebandNLL(x, pdf_list, histogram, cat, eps, offset, strategy, range, className, sideBand = True, FT = True)
```
where `x` is the fitting axis, `pdf_list` is the list of the functions that you wnat to do F test over, `histogram` is the `RooDataHist` to be fit, `cat` is the category name, `eps` is the EDM, `offset` controls whether using the offset, and `strategy` should be 0, `range` specicies the range to plot, and `className` is the string added to the title of the multi-function plot output file. `sideBand` and `FT` are always set to true.

When calling the F test, NLL fit is used by default. A command to run the F test in ggf1 is
```
bash$ ./F_test.py -c ggf1 -con ../Config/chi2_config_xgboost_nodrop.json
```
Background models tested are selected based on the `"SST"` value in the config file.Output has both the test statistics and the P-values, as well as a multi-function plot.
## Make Workspaces 
To operate the next test, we need to make toys with `combine`, thus need datacards and worksapces. To generate a workspace for ggf1, for example,
```
bash$ ./make_workspace_bias_test.py -c ggf1 -a 0 -con ../Config/chi2_config_xgboost_nodrop.json
```
where `-a` means whether you want to generate and store an Asimov histogram in the workspace. This helps us to cross check with `combine` expected significance calculation.

To generate the signal workspace only, run
```
bash$ ./make_signal_workspace.py -c ggf1 -conB ../Config/chi2_config_xgboost_nodrop.json -conS ../Config/config_DSCB.json
```

## CMS Bias Test
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

