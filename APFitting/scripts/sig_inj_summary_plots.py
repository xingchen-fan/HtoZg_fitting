import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def nice_hist(column,hist_name,xlab="",nbins=50,mu_equals_one=-1):
  fig = plt.figure(figsize=(10,8))
  average_column = np.mean(column)
  std_dev = np.std(column)
  plt.hist(column,bins=nbins,range=[int(min(column)), int(max(column))], histtype=u'step')
  plt.axvline(average_column,color='b')
  plt.axvline(average_column + std_dev, color='b',linestyle='--')
  plt.axvline(average_column - std_dev, color='b',linestyle='--')

  if(mu_equals_one > 0):
    plt.axvline(mu_equals_one,color='r')

  print(hist_name + ": " + str(average_column) + "+-" + str(std_dev))
  plt.xlabel(xlab)
  plt.ylabel("Toys/" + str( (max(column)-min(column))/nbins ) )
  fig.savefig("./output/siginjection_summary/" + hist_name)

  return average_column,std_dev

def line_plot(x, y, y_unc, hist_name, x_label='',y_label=''):
  fig = plt.figure(figsize=(10,8))

  plt.plot(x,y,linestyle='-',color='r')
  plt.plot(x,y,'o',color='r')

  plt.plot(x,y+y_unc,linestyle='--',color='r')
  plt.plot(x,y-y_unc,linestyle='--',color='r')
  plt.plot(x,x,linestyle='-',color='b')

  plt.xlabel(x_label, fontsize=24)
  plt.ylabel(y_label, fontsize=24)
  fig.savefig("./output/siginjection_summary/" + hist_name)

def line_bkg_plot(x, y, y_unc, x_inj,x_inj_unc, hist_name, x_label='',y_label=''):
  fig = plt.figure(figsize=(10,8))

  plt.plot(x,y,linestyle='-',color='r')
  plt.plot(x,y,'o',color='r')

  plt.plot(x,y+y_unc,linestyle='--',color='r')
  plt.plot(x,y-y_unc,linestyle='--',color='r')
  plt.axhline(x_inj,x[0],x[-1],linestyle='-',color='b')
  plt.axhline(x_inj+x_inj_unc,x[0],x[-1],linestyle='-',color='b')
  plt.axhline(x_inj-x_inj_unc,x[0],x[-1],linestyle='-',color='b')

  plt.xlabel(x_label, fontsize=24)
  plt.ylabel(y_label, fontsize=24)
  fig.savefig("./output/siginjection_summary/" + hist_name)


if __name__=="__main__":
              #[ggF,   VBF, ttH_lep, WH_3l, ZH_MET, ttH_had]
  sig_yields = {'ggF':167.6, 'VBF':37.0, 'ttH_lep':0.40, 'WH_3l':0.85, 'ZH_MET':0.51, 'ttH_had':0.67}

  #if not os.listdir("./output/siginjection_summary"):
  #  os.mkdir("./output/siginjection_summary")
  #else:
  #  print("Directory \"./output/siginjection_summary\" exists.")
 
  #NsigMC,Nsiginj_mean,Nsiginj_std_dev,Nbkginj_mean,Nbkginj_std_dev,Nsigfit_mean,Nsigfit_std_dev,
  #Nsig_err_mean,Nsig_err_std_dev,Nbkgfit_mean,Nbkgfit_std_dev,Nbkgfit_blinded_mean,Nbkgfit_blinded_std_dev,
  #Nbkgfit_err_mean,Nbkgfit_err_std_dev,Nbonly_fit_mean,Nbonly_fit_std_dev,NBonly_err_mean,
  #NBonly_err_std_dev,Nll_splusb_mean,Nll_splusb_std_dev,Nll_bkgonly_mean,Nll_bkgonly_std_dev
  
  categories=['VBF', 'ttH_lep', 'WH_3l', 'ZH_MET', 'ttH_had']#'ggF', 

  for cat in categories:
    #With only sig region generated
    #summary_name = "./siginjection_csv/csv_siginjection_summary_" + cat + ".csv"
    #output_name = 'sig_injection_fit_yields_bsig_' + cat + '.pdf'

    #With full range generated
    summary_name = './siginjection_csv/csv_siginjection_ffrgen_summary_' + cat + '.csv'
    output_name = 'sig_injection_fit_yields_ffrgen_' + cat + '.pdf'

    df = pd.read_csv(summary_name) 
    line_plot(df['NsigMC'], df['Nsigfit_mean'], df['Nsigfit_std_dev'], output_name, '$N_{S,MC}$', 'Average $N_{S,fit}$')

    print(df['Nsigfit_mean'])
    print(df['Nsig_err_mean'])

    #line_plot(df['NsigMC']/sig_yields[cat], df['Nsigfit_mean']/df['NsigMC'], df['Nsigfit_std_dev']/df['NsigMC'], 'sig_injection_relative_fit_yields_' + cat + '.pdf', '$N_{S,MC}$', 'Average $N_{S,fit}$/$N_{S,MC}$')
    #line_plot(df['NsigMC']/sig_yields[cat], (df['Nsigfit_mean']-df['NsigMC'])/df['NsigMC'], df['Nsigfit_std_dev']/df['NsigMC'], 'sig_injection_relative_diff_fit_yields_' + cat + '.pdf', '$N_{S,MC}$', '(Average $N_{S,fit}$ - $N_{S,MC}$)/$N_{S,MC}$')

    #line_bkg_plot(df['NsigMC']/sig_yields[cat], df['Nbkgfit_blinded_mean'], df['Nbkgfit_blinded_std_dev'], df['Nbkginj_mean'], df['Nbkginj_std_dev'], 'sig_injection_bkg_fit_yields_' + cat + '.pdf','$N_{bkg}$')
     
  

 # make_plots("csv_ttH_had_agg_sig_strength_1_extended.csv","test_ttH_had_agg_1",'ttH_had',1)




  #make_plots("csv_ttH_had_agg_sig_strength_5_extended.csv","ttH_had_agg_5") 
  #dataframe_test = pd.read_csv("./siginjection_csv/csv_ttH_had_agg_sig_strength_1_extended.csv")
  #print(dataframe_test.head())
  #variables = list(dataframe_test.columns.values)
  #print(variables)
  #nice_hist(dataframe_test['Nsigfit'],"test_hisogram",variables[0])

  #for var in variables:
  #  if(var=="Nsigfit" or var=="Nsiginj"):
  #    nice_hist(dataframe_test[var],"test_hisogram_" + var,var,50,sig_yields['ttH_had'])
  #  else:
  #    nice_hist(dataframe_test[var],"test_hisogram_" + var,var)

  #nice_hist()
  
