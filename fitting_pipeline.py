#!/usr/bin/env python3
import os
import sys
import ROOT
import argparse
import json
import subprocess
import time
import logging

#---------------------------How to use----------------------------#
#Run this from the main HtoZg_fitting directory. 
#Performs the whole fitting process end-to-end (ideally. . .)
#CURRENTLY INCOMPLETE



dataSignalBlinding = True
mcSignalBlinding = True
asimovData = False
asimovMC = False

genNewSignalFcts = True
tsStart = time.perf_counter()

#inputs for asimov vs real for both data and MC?

#----------------------Generate signal models-----------------------#
#Split by lepton flavor, and by BDT category.
#This script can be run as a subprocess inside other function, or directly from command line.
settingsFile = '/afs/cern.ch/work/j/jgrassi/HtoZg_fitting/Config/pipeline_settings.json'
jfile = open(settingsFile)
configs = json.load(jfile)
sigSettings = configs["signal_settings"]
logFile = sigSettings["log_file"]

print(logFile)
logging.basicConfig(filename=logFile, level=logging.INFO)

timeStamp = str(time.time())
logging.info("Beginning signal model generation")

fitArgs = ['./Signal_model_preparation/signal_fit_jdg.py','-j', settingsFile]
with open(logFile, 'a') as logfile:
    subprocess.run(fitArgs, stdout=logfile, stderr=logfile)
logging.info("...........................")

tsGenSig = time.perf_counter()
logging.info(f"Generated signal models in {tsGenSig-tsStart:0.3f} seconds")
logging.info("...........................")


"""
#----------------------Generate background models--------------------#
#Split by BDT category, using convolution model, many models generated.
testLaurent = False
if testLaurent == True:
    logging.info("Beginning to test powers for Laurent background models")
    lauArgs = ['./Chi2_test/lau2_power_finder.py', '-sr', '110', '-tr', '65', '-hi', '-2', '-lo', '-6', '-NLL', '0', '-ts', '1', '-mt', '2', '-js', settingsFile] 
    with open(logFile, 'a') as logfile:
        subprocess.run(fitArgs, stdout=logfile, stderr=logfile)

    tsTestLau = time.perf_counter()
    logging.info("...........................")
    logging.info(f"Finished testing powers for Laurent background models in {tsTestLau-tsGenSig:0.3f} seconds")
    logging.info("...........................")



tsGenBkg = time.perf_counter()
logging.info(f"Generated bkg models in {tsGenBkg-tsGenSig:0.3f} seconds")

#Goodness of fit test (is this always using the signal blinded data?)
tsGoodFit = time.perf_counter()
logging.info(f"Ran goodness of fit test in {tsGoodFit-tsGenBkg:0.3f} seconds")

#F-test (always on signal blinded data?)
tsFTest = time.perf_counter()
logging.info(f"Ran F-test in {tsFTest-tsGoodFit:0.3f} seconds")

#bias test (Toy data sample and MC) on spurious signal and envelope bias. 
tsBias = time.perf_counter()
logging.info(f"Ran bias tests in {tsBias-tsFTest:0.3f} seconds")

#envelope fit (on either asimov or true sample)
tsEnvelope = time.perf_counter()
logging.info(f"Generated signal models in {tsEnvelope-tsBias:0.3f} seconds")
"""


