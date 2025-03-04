import ROOT

# Cornell sample in dat format

class readDat: #Cornell MC sampleS
    def __init__(self, list, dir="", data_dir = "/afs/cern.ch/user/f/fanx/rui_newdat/"):
        x = list[0]
        y = list[1]
        w = list[2]
        bdt = list[3]
        year = list[4]
        lep = list[5]
        ph_eta = list[6]
        nlep = list[7]
        njet = list[8]

        # Read dat samples
        #bkg_run2 = ROOT.RooDataSet.read("SMZg_deathvalley_v3_untagged.dat,DY_deathvalley_v3_untagged.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "", dir)
        #sig_run2 = ROOT.RooDataSet.read("FullSig_deathvalley_v3_untagged.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "", dir)
        #tot_run2 = ROOT.RooDataSet.read("SMZg_deathvalley_v3_untagged.dat,DY_deathvalley_v3_untagged.dat,FullSig_deathvalley_v3_untagged.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "", dir)
        bkg_run2 = ROOT.RooDataSet.read("SMZg_ggF1_newbdt_trig.dat,SMZg_ggF2_newbdt_trig.dat,SMZg_ggF3_newbdt_trig.dat,SMZg_ggF4_newbdt_trig.dat,DY_ggF1_newbdt_trig.dat,DY_ggF2_newbdt_trig.dat,DY_ggF3_newbdt_trig.dat,DY_ggF4_newbdt_trig.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet),"", dir)
        sig_run2 = ROOT.RooDataSet.read("Signal_ggF1_newbdt_trig.dat,Signal_ggF2_newbdt_trig.dat,Signal_ggF3_newbdt_trig.dat,Signal_ggF4_newbdt_trig.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "", dir)
        tot_run2 = ROOT.RooDataSet.read("SMZg_ggF1_newbdt_trig.dat,SMZg_ggF2_newbdt_trig.dat,SMZg_ggF3_newbdt_trig.dat,SMZg_ggF4_newbdt_trig.dat,DY_ggF1_newbdt_trig.dat,DY_ggF2_newbdt_trig.dat,DY_ggF3_newbdt_trig.dat,DY_ggF4_newbdt_trig.dat,Signal_ggF1_newbdt_trig.dat,Signal_ggF2_newbdt_trig.dat,Signal_ggF3_newbdt_trig.dat,Signal_ggF4_newbdt_trig.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "", dir)
        bkg_SM = ROOT.RooDataSet.read("SMZg_ggF1_newbdt_trig_FullSample.dat,SMZg_ggF2_newbdt_trig_FullSample.dat,SMZg_ggF3_newbdt_trig_FullSample.dat,SMZg_ggF4_newbdt_trig_FullSample.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet),"", dir)


        data = ROOT.RooDataSet.read("data_ggF1_newbdt_trig.dat,data_ggF2_newbdt_trig.dat,data_ggF3_newbdt_trig.dat,data_ggF4_newbdt_trig.dat", ROOT.RooArgList(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "", data_dir)


        bdt1 = -0.5 #-0.36
        bdt2 = 0.02 #-0.06
        bdt3 = 0.34 #0.08
        bdt4 = 0.66 #0.15

        # Make RooDataSet
        # u1_bkg_run2_el = ROOT.RooDataSet("u1_bkg_run2_el", "u1_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        # u2_bkg_run2_el = ROOT.RooDataSet("u2_bkg_run2_el", "u2_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        # u3_bkg_run2_el = ROOT.RooDataSet("u3_bkg_run2_el", "u3_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        # u4_bkg_run2_el = ROOT.RooDataSet("u4_bkg_run2_el", "u4_bkg_run2_el", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        # u1_sig_run2_el = ROOT.RooDataSet("u1_sig_run2_el", "u1_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        # u2_sig_run2_el = ROOT.RooDataSet("u2_sig_run2_el", "u2_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        # u3_sig_run2_el = ROOT.RooDataSet("u3_sig_run2_el", "u3_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        # u4_sig_run2_el = ROOT.RooDataSet("u4_sig_run2_el", "u4_sig_run2_el", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        # u1_tot_run2_el = ROOT.RooDataSet("u1_tot_run2_el", "u1_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        # u2_tot_run2_el = ROOT.RooDataSet("u2_tot_run2_el", "u2_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        # u3_tot_run2_el = ROOT.RooDataSet("u3_tot_run2_el", "u3_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        # u4_tot_run2_el = ROOT.RooDataSet("u4_tot_run2_el", "u4_tot_run2_el", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep < 1 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        # u1_bkg_run2_mu = ROOT.RooDataSet("u1_bkg_run2_mu", "u1_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        # u2_bkg_run2_mu = ROOT.RooDataSet("u2_bkg_run2_mu", "u2_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        # u3_bkg_run2_mu = ROOT.RooDataSet("u3_bkg_run2_mu", "u3_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        # u4_bkg_run2_mu = ROOT.RooDataSet("u4_bkg_run2_mu", "u4_bkg_run2_mu", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        # u1_sig_run2_mu = ROOT.RooDataSet("u1_sig_run2_mu", "u1_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        # u2_sig_run2_mu = ROOT.RooDataSet("u2_sig_run2_mu", "u2_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        # u3_sig_run2_mu = ROOT.RooDataSet("u3_sig_run2_mu", "u3_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        # u4_sig_run2_mu = ROOT.RooDataSet("u4_sig_run2_mu", "u4_sig_run2_mu", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        # u1_tot_run2_mu = ROOT.RooDataSet("u1_tot_run2_mu", "u1_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        # u2_tot_run2_mu = ROOT.RooDataSet("u2_tot_run2_mu", "u2_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        # u3_tot_run2_mu = ROOT.RooDataSet("u3_tot_run2_mu", "u3_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        # u4_tot_run2_mu = ROOT.RooDataSet("u4_tot_run2_mu", "u4_tot_run2_mu", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "lep > 0 && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        u1_bkg_run2 = ROOT.RooDataSet("u1_bkg_run2", "u1_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        u2_bkg_run2 = ROOT.RooDataSet("u2_bkg_run2", "u2_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        u3_bkg_run2 = ROOT.RooDataSet("u3_bkg_run2", "u3_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        u4_bkg_run2 = ROOT.RooDataSet("u4_bkg_run2", "u4_bkg_run2", bkg_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        u1_sig_run2 = ROOT.RooDataSet("u1_sig_run2", "u1_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt >" + str(bdt4), "w")
        u2_sig_run2 = ROOT.RooDataSet("u2_sig_run2", "u2_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        u3_sig_run2 = ROOT.RooDataSet("u3_sig_run2", "u3_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        u4_sig_run2 = ROOT.RooDataSet("u4_sig_run2", "u4_sig_run2", sig_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")

        u1_tot_run2 = ROOT.RooDataSet("u1_tot_run2", "u1_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        u2_tot_run2 = ROOT.RooDataSet("u2_tot_run2", "u2_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4) , "w")
        u3_tot_run2 = ROOT.RooDataSet("u3_tot_run2", "u3_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3) , "w")
        u4_tot_run2 = ROOT.RooDataSet("u4_tot_run2", "u4_tot_run2", tot_run2, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2) , "w")

        u1_bkg_SM = ROOT.RooDataSet("u1_bkg_SM", "u1_bkg_SM", bkg_SM, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt4), "w")
        u2_bkg_SM = ROOT.RooDataSet("u2_bkg_SM", "u2_bkg_SM", bkg_SM, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4), "w")
        u3_bkg_SM = ROOT.RooDataSet("u3_bkg_SM", "u3_bkg_SM", bkg_SM, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3), "w")
        u4_bkg_SM = ROOT.RooDataSet("u4_bkg_SM", "u4_bkg_SM", bkg_SM, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2), "w")


        u1_data = ROOT.RooDataSet("u1_data", "u1_data", data, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt4))
        u2_data = ROOT.RooDataSet("u2_data", "u2_data", data, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt3) +" && bdt < " + str(bdt4))
        u3_data = ROOT.RooDataSet("u3_data", "u3_data", data, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt2) +" && bdt < " + str(bdt3))
        u4_data = ROOT.RooDataSet("u4_data", "u4_data", data, ROOT.RooArgSet(x, y, bdt, w, year, lep, ph_eta, nlep, njet), "CMS_hzg_mass > 100. && nlep <= 2 && njet < 2 && bdt > " + str(bdt1) +" && bdt < " + str(bdt2))

        self.data_hist_u1_SM =  ROOT.RooDataHist("data_hist_u1_SM", "data_hist_u1_SM", x, u1_bkg_SM)
        self.data_hist_u2_SM =  ROOT.RooDataHist("data_hist_u2_SM", "data_hist_u2_SM", x, u2_bkg_SM)
        self.data_hist_u3_SM =  ROOT.RooDataHist("data_hist_u3_SM", "data_hist_u3_SM", x, u3_bkg_SM)
        self.data_hist_u4_SM =  ROOT.RooDataHist("data_hist_u4_SM", "data_hist_u4_SM", x, u4_bkg_SM)

        self.data_hist_untagged1 = ROOT.RooDataHist("data_hist_untagged1", "data_hist_untagged1", x, u1_tot_run2)
        self.data_hist_untagged2 = ROOT.RooDataHist("data_hist_untagged2", "data_hist_untagged2", x, u2_tot_run2)
        self.data_hist_untagged3 = ROOT.RooDataHist("data_hist_untagged3", "data_hist_untagged3", x, u3_tot_run2)
        self.data_hist_untagged4 = ROOT.RooDataHist("data_hist_untagged4", "data_hist_untagged4", x, u4_tot_run2)

        self.data_hist_untagged1_bkg = ROOT.RooDataHist("data_hist_untagged1_bkg", "data_hist_untagged1_bkg", x, u1_bkg_run2)
        self.data_hist_untagged2_bkg = ROOT.RooDataHist("data_hist_untagged2_bkg", "data_hist_untagged2_bkg", x, u2_bkg_run2)
        self.data_hist_untagged3_bkg = ROOT.RooDataHist("data_hist_untagged3_bkg", "data_hist_untagged3_bkg", x, u3_bkg_run2)
        self.data_hist_untagged4_bkg = ROOT.RooDataHist("data_hist_untagged4_bkg", "data_hist_untagged4_bkg", x, u4_bkg_run2)

        self.data_hist_untagged1_sig = ROOT.RooDataHist("data_hist_untagged1_sig", "data_hist_untagged1_sig", x, u1_sig_run2)
        self.data_hist_untagged2_sig = ROOT.RooDataHist("data_hist_untagged2_sig", "data_hist_untagged2_sig", x, u2_sig_run2)
        self.data_hist_untagged3_sig = ROOT.RooDataHist("data_hist_untagged3_sig", "data_hist_untagged3_sig", x, u3_sig_run2)
        self.data_hist_untagged4_sig = ROOT.RooDataHist("data_hist_untagged4_sig", "data_hist_untagged4_sig", x, u4_sig_run2)
        
        self.data_u1 = ROOT.RooDataHist("data_u1", "data_u1", x, u1_data)
        self.data_u2 = ROOT.RooDataHist("data_u2", "data_u2", x, u2_data)
        self.data_u3 = ROOT.RooDataHist("data_u3", "data_u3", x, u3_data)
        self.data_u4 = ROOT.RooDataHist("data_u4", "data_u4", x, u4_data)

    def numCheck(self):
        print("# bin1 bkg = ", self.data_hist_untagged1_bkg.sumEntries())
        print("# bin2 bkg = ", self.data_hist_untagged2_bkg.sumEntries())
        print("# bin3 bkg = ", self.data_hist_untagged3_bkg.sumEntries())
        print("# bin4 bkg = ", self.data_hist_untagged4_bkg.sumEntries())

    def dataNumCheck(self):
        print("# bin1 data = ", self.data_u1.sumEntries())
        print("# bin2 data = ", self.data_u2.sumEntries())
        print("# bin3 data = ", self.data_u3.sumEntries())
        print("# bin4 data = ", self.data_u4.sumEntries())


class readRoot: #Beijing data sample
    def __init__(self, x, dir=""):
        file = ROOT.TFile(dir)
        tree = file.Get("test")
        bdt = 0.
        mllg = 0.

        nentries = tree.GetEntries()

        bin1 = ROOT.RooDataSet("bin1", "bin1", ROOT.RooArgSet(x))
        bin2 = ROOT.RooDataSet("bin2", "bin2", ROOT.RooArgSet(x))
        bin3 = ROOT.RooDataSet("bin3", "bin3", ROOT.RooArgSet(x))
        bin4 = ROOT.RooDataSet("bin4", "bin4", ROOT.RooArgSet(x))

        for i in range(nentries):
            tree.GetEntry(i)
            bdt = tree.bdt_score_t
            mllg = tree.H_mass
            if (i%1000 == 0): print(bdt)
            x.setVal(mllg)
            if (mllg > 105. and mllg < 170. and bdt > 0. and bdt <= 0.29):
                bin1.add(ROOT.RooArgSet(x))
            if (mllg > 105. and mllg < 170. and bdt > 0.29 and bdt <= 0.57):
                bin2.add(ROOT.RooArgSet(x))
            if (mllg > 105. and mllg < 170. and bdt > 0.57 and bdt <= 0.73):
                bin3.add(ROOT.RooArgSet(x))
            if (mllg > 105. and mllg < 170. and bdt > 0.73 and bdt <= 1.):
                bin4.add(ROOT.RooArgSet(x))

        self.data_hist_bin1 = ROOT.RooDataHist("data_hist_bin1", "data_hist_bin1", x, bin1)
        self.data_hist_bin2 = ROOT.RooDataHist("data_hist_bin2", "data_hist_bin2", x, bin2)
        self.data_hist_bin3 = ROOT.RooDataHist("data_hist_bin3", "data_hist_bin3", x, bin3)
        self.data_hist_bin4 = ROOT.RooDataHist("data_hist_bin4", "data_hist_bin4", x, bin4)

    def numCheck(self):
        print("# bin1 bkg = ", self.data_hist_bin1.sumEntries())
        print("# bin2 bkg = ", self.data_hist_bin2.sumEntries())
        print("# bin3 bkg = ", self.data_hist_bin3.sumEntries())
        print("# bin4 bkg = ", self.data_hist_bin4.sumEntries())
        

class readWsp: #Zebing core func ggf samples
    def __init__(self, x, dir='', dataname=''):
        file = ROOT.TFile(dir)
        wsp = file.Get("CMS_hzg_workspace")
        data = wsp.data(dataname)
        moddata = ROOT.RooDataSet("data_zebing_"+dataname, "data_zebing"+dataname, data, ROOT.RooArgSet(x), "CMS_hzg_mass > "+ str(x.getMin())+ " && CMS_hzg_mass < " + str(x.getMax()))
        self.hist = ROOT.RooDataHist(dataname, dataname, x, moddata)
    
    def numCheck(self):
        print ('n events = ', self.hist.sumEntries())

class readWspMC: #Zebing core func signal samples
    def __init__(self, x , dir='', dataname=''):
        file = ROOT.TFile(dir)
        wsp = file.Get("CMS_hzg_workspace")
        data = wsp.data(dataname)
        #moddata = ROOT.RooDataSet("mc_zebing_"+dataname, "mc_zebing"+dataname, data, ROOT.RooArgSet(x, w), "", "CMS_hzg_weight")
        self.hist = ROOT.RooDataHist("mc_hist_zebing"+dataname, "mc_hist_zebing"+dataname, x, data)
    def numCheck(self):
        print ('n events = ', self.hist.sumEntries())

