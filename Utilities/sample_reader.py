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


def readRuiROOTdata(x, direct='', bdt_low=0, bdt_high=0, cat=''):
    chain = ROOT.TChain('TreeB')
    chain.Add(direct + 'data_2016APV_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2016_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2017_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2018_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2022EE_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2022_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2023_pinnacles_fullrange.root')
    chain.Add(direct + 'data_2023BPix_pinnacles_fullrange.root')
    hist_TH1 = ROOT.TH1F(cat+'_th1f', cat+'_th1f', 340, 95, 180)
    for entry in chain:
        if entry.bdt_score_test > bdt_low and entry.bdt_score_test < bdt_high:
            if cat == 'ggf1' or cat == 'ggf2' or cat == 'ggf3' or cat == 'ggf4' or cat == 'ggF1' or cat == 'ggF2' or cat == 'ggF3' or cat == 'ggF4':
                if entry.met < 90:
                    hist_TH1.Fill(entry.llphoton_refit_m)
            else:
                hist_TH1.Fill(entry.llphoton_refit_m)
    hist_data = ROOT.RooDataHist('hist_data', 'hist_data', x, hist_TH1)
    return hist_data
class readRuiROOTggFdata:
    def __init__(self, x, direct='', bdt1=0, bdt2=0, bdt3=0):
        chain = ROOT.TChain('outtree')
        # Old files are called data_YEAR_pinnacles_ggf_fixed.root
        # Tree is called TreeB
        # BDT score branch called bdt_score_test
        """
        chain.Add(direct + 'data_2016APV_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2016_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2017_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2018_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2022EE_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2022_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2023_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'data_2023BPix_pinnacles_ggf_fixed.root')
        """
        chain.Add(direct + 'data_2016APV_output.root')
        chain.Add(direct + 'data_2016_output.root')
        chain.Add(direct + 'data_2017_output.root')
        chain.Add(direct + 'data_2018_output.root')
        chain.Add(direct + 'data_2022EE_output.root')
        chain.Add(direct + 'data_2022_output.root')
        chain.Add(direct + 'data_2023_output.root')
        chain.Add(direct + 'data_2023BPix_output.root')

        hist1 = ROOT.TH1F('ggf1_th1f', 'ggf1_th1f', 340, 95, 180)
        hist2 = ROOT.TH1F('ggf2_th1f', 'ggf2_th1f', 340, 95, 180)
        hist3 = ROOT.TH1F('ggf3_th1f', 'ggf3_th1f', 340, 95, 180)
        hist4 = ROOT.TH1F('ggf4_th1f', 'ggf4_th1f', 340, 95, 180)
        hist1EL_TH1 = ROOT.TH1F('ggf1el_th1f', 'ggf1el_th1f', 340, 95, 180)
        hist2EL_TH1 = ROOT.TH1F('ggf2el_th1f', 'ggf2el_th1f', 340, 95, 180)
        hist3EL_TH1 = ROOT.TH1F('ggf3el_th1f', 'ggf3el_th1f', 340, 95, 180)
        hist4EL_TH1 = ROOT.TH1F('ggf4el_th1f', 'ggf4el_th1f', 340, 95, 180)
        hist1MU_TH1 = ROOT.TH1F('ggf1mu_th1f', 'ggf1mu_th1f', 340, 95, 180)
        hist2MU_TH1 = ROOT.TH1F('ggf2mu_th1f', 'ggf2mu_th1f', 340, 95, 180)
        hist3MU_TH1 = ROOT.TH1F('ggf3mu_th1f', 'ggf3mu_th1f', 340, 95, 180)
        hist4MU_TH1 = ROOT.TH1F('ggf4mu_th1f', 'ggf4mu_th1f', 340, 95, 180)
        for entry in chain:
            if entry.ll_lepid == 11 and entry.njet < 2:
                if entry.BDT_score > bdt1 and entry.met < 90:
                    hist1EL_TH1.Fill(entry.llphoton_refit_m)
                    hist1.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.met < 90:
                    hist2EL_TH1.Fill(entry.llphoton_refit_m)
                    hist2.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.met < 90:
                    hist3EL_TH1.Fill(entry.llphoton_refit_m)
                    hist3.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.met < 90:
                    hist4EL_TH1.Fill(entry.llphoton_refit_m)
                    hist4.Fill(entry.llphoton_refit_m)
            elif entry.ll_lepid == 13 and entry.njet < 2:
                if entry.BDT_score > bdt1 and entry.met < 90:
                    hist1MU_TH1.Fill(entry.llphoton_refit_m)
                    hist1.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.met < 90:
                    hist2MU_TH1.Fill(entry.llphoton_refit_m)
                    hist2.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.met < 90:
                    hist3MU_TH1.Fill(entry.llphoton_refit_m)
                    hist3.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.met < 90:
                    hist4MU_TH1.Fill(entry.llphoton_refit_m)
                    hist4.Fill(entry.llphoton_refit_m)

        self.ggf1_el = ROOT.RooDataHist('hist_ggf1_data_el', 'hist_ggf1_data_el', x, hist1EL_TH1)
        self.ggf2_el = ROOT.RooDataHist('hist_ggf2_data_el', 'hist_ggf2_data_el', x, hist2EL_TH1)
        self.ggf3_el = ROOT.RooDataHist('hist_ggf3_data_el', 'hist_ggf3_data_el', x, hist3EL_TH1)
        self.ggf4_el = ROOT.RooDataHist('hist_ggf4_data_el', 'hist_ggf4_data_el', x, hist4EL_TH1)
        
        self.ggf1_mu = ROOT.RooDataHist('hist_ggf1_data_mu', 'hist_ggf1_data_mu', x, hist1MU_TH1)
        self.ggf2_mu = ROOT.RooDataHist('hist_ggf2_data_mu', 'hist_ggf2_data_mu', x, hist2MU_TH1)
        self.ggf3_mu = ROOT.RooDataHist('hist_ggf3_data_mu', 'hist_ggf3_data_mu', x, hist3MU_TH1)
        self.ggf4_mu = ROOT.RooDataHist('hist_ggf4_data_mu', 'hist_ggf4_data_mu', x, hist4MU_TH1)

        self.ggf1 = ROOT.RooDataHist('hist_ggf1_data', 'hist_ggf1_data', x, hist1)
        self.ggf2 = ROOT.RooDataHist('hist_ggf2_data', 'hist_ggf2_data', x, hist2)
        self.ggf3 = ROOT.RooDataHist('hist_ggf3_data', 'hist_ggf3_data', x, hist3)
        self.ggf4 = ROOT.RooDataHist('hist_ggf4_data', 'hist_ggf4_data', x, hist4)

class readRuiROOTVBFdata:
    def __init__(self, x, direct='', bdt1=0, bdt2=0, bdt3=0):
        # Old files, BDT branch is called BDT_score_2j
        chain = ROOT.TChain('outtree')
        chain.Add(direct + 'data_2016APV_output.root')
        chain.Add(direct + 'data_2016_output.root')
        chain.Add(direct + 'data_2017_output.root')
        chain.Add(direct + 'data_2018_output.root')
        chain.Add(direct + 'data_2022EE_output.root')
        chain.Add(direct + 'data_2022_output.root')
        chain.Add(direct + 'data_2023_output.root')
        chain.Add(direct + 'data_2023BPix_output.root')
        hist1_TH1 = ROOT.TH1F('vbf1_th1f', 'vbf1_th1f', 360, 90, 180)
        hist2_TH1 = ROOT.TH1F('vbf2_th1f', 'vbf2_th1f', 360, 90, 180)
        hist3_TH1 = ROOT.TH1F('vbf3_th1f', 'vbf3_th1f', 360, 90, 180)
        hist4_TH1 = ROOT.TH1F('vbf4_th1f', 'vbf4_th1f', 360, 90, 180)
        hist12_TH1 = ROOT.TH1F('vbf12_th1f', 'vbf12_th1f', 360, 90, 180)
        
        hist1_TH1_el = ROOT.TH1F('vbf1_th1f_el', 'vbf1_th1f_el', 360, 90, 180)
        hist2_TH1_el = ROOT.TH1F('vbf2_th1f_el', 'vbf2_th1f_el', 360, 90, 180)
        hist3_TH1_el = ROOT.TH1F('vbf3_th1f_el', 'vbf3_th1f_el', 360, 90, 180)
        hist4_TH1_el = ROOT.TH1F('vbf4_th1f_el', 'vbf4_th1f_el', 360, 90, 180)
        hist12_TH1_el = ROOT.TH1F('vbf12_th1f_el', 'vbf12_th1f_el', 360, 90, 180)
        hist1_TH1_mu = ROOT.TH1F('vbf1_th1f_mu', 'vbf1_th1f_mu', 360, 90, 180)
        hist2_TH1_mu = ROOT.TH1F('vbf2_th1f_mu', 'vbf2_th1f_mu', 360, 90, 180)
        hist3_TH1_mu = ROOT.TH1F('vbf3_th1f_mu', 'vbf3_th1f_mu', 360, 90, 180)
        hist4_TH1_mu = ROOT.TH1F('vbf4_th1f_mu', 'vbf4_th1f_mu', 360, 90, 180)
        hist12_TH1_mu = ROOT.TH1F('vbf12_th1f_mu', 'vbf12_th1f_mu', 360, 90, 180)
        
        for entry in chain:
            if entry.ll_lepid == 11 and entry.njet > 1 and entry.nbdfm == 0:
                if entry.BDT_score > bdt1:
                    hist1_TH1.Fill(entry.llphoton_refit_m)
                    hist1_TH1_el.Fill(entry.llphoton_refit_m)
                    hist12_TH1.Fill(entry.llphoton_refit_m)
                    hist12_TH1_el.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 :
                    hist2_TH1.Fill(entry.llphoton_refit_m)
                    hist2_TH1_el.Fill(entry.llphoton_refit_m)
                    hist12_TH1.Fill(entry.llphoton_refit_m)
                    hist12_TH1_el.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 :
                    hist3_TH1.Fill(entry.llphoton_refit_m)
                    hist3_TH1_el.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > -1 and entry.BDT_score < bdt3 :
                    hist4_TH1.Fill(entry.llphoton_refit_m)
                    hist4_TH1_el.Fill(entry.llphoton_refit_m)
            if entry.ll_lepid == 13 and entry.njet > 1 and entry.nbdfm == 0:
                if entry.BDT_score > bdt1:
                    hist1_TH1.Fill(entry.llphoton_refit_m)
                    hist1_TH1_mu.Fill(entry.llphoton_refit_m)
                    hist12_TH1.Fill(entry.llphoton_refit_m)
                    hist12_TH1_mu.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 :
                    hist2_TH1.Fill(entry.llphoton_refit_m)
                    hist2_TH1_mu.Fill(entry.llphoton_refit_m)
                    hist12_TH1.Fill(entry.llphoton_refit_m)
                    hist12_TH1_mu.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 :
                    hist3_TH1.Fill(entry.llphoton_refit_m)
                    hist3_TH1_mu.Fill(entry.llphoton_refit_m)
                elif entry.BDT_score > -1 and entry.BDT_score < bdt3 :
                    hist4_TH1.Fill(entry.llphoton_refit_m)
                    hist4_TH1_mu.Fill(entry.llphoton_refit_m)

        self.vbf1 = ROOT.RooDataHist('hist_vbf1_data', 'hist_vbf1_data', x, hist1_TH1)
        self.vbf2 = ROOT.RooDataHist('hist_vbf2_data', 'hist_vbf2_data', x, hist2_TH1)
        self.vbf3 = ROOT.RooDataHist('hist_vbf3_data', 'hist_vbf3_data', x, hist3_TH1)
        self.vbf4 = ROOT.RooDataHist('hist_vbf4_data', 'hist_vbf4_data', x, hist4_TH1)
        self.vbf12 = ROOT.RooDataHist('hist_vbf12_data', 'hist_vbf12_data', x, hist12_TH1)
        
        self.vbf1_el = ROOT.RooDataHist('hist_vbf1_data_el', 'hist_vbf1_data_el', x, hist1_TH1_el)
        self.vbf2_el = ROOT.RooDataHist('hist_vbf2_data_el', 'hist_vbf2_data_el', x, hist2_TH1_el)
        self.vbf3_el = ROOT.RooDataHist('hist_vbf3_data_el', 'hist_vbf3_data_el', x, hist3_TH1_el)
        self.vbf4_el = ROOT.RooDataHist('hist_vbf4_data_el', 'hist_vbf4_data_el', x, hist4_TH1_el)
        self.vbf12_el = ROOT.RooDataHist('hist_vbf12_data_el', 'hist_vbf12_data_el', x, hist12_TH1_el)
        
        self.vbf1_mu = ROOT.RooDataHist('hist_vbf1_data_mu', 'hist_vbf1_data_mu', x, hist1_TH1_mu)
        self.vbf2_mu = ROOT.RooDataHist('hist_vbf2_data_mu', 'hist_vbf2_data_mu', x, hist2_TH1_mu)
        self.vbf3_mu = ROOT.RooDataHist('hist_vbf3_data_mu', 'hist_vbf3_data_mu', x, hist3_TH1_mu)
        self.vbf4_mu = ROOT.RooDataHist('hist_vbf4_data_mu', 'hist_vbf4_data_mu', x, hist4_TH1_mu)
        self.vbf12_mu = ROOT.RooDataHist('hist_vbf12_data_mu', 'hist_vbf12_data_mu', x, hist12_TH1_mu)
        
class readRuiROOTggFSignal:
    def __init__(self, x, direct='', bdt1=0, bdt2=0, bdt3=0):
        year = ['2016', '2016APV', '2017', '2018', '2022', '2022EE', '2023', '2023BPix']
        prod = ['GGF', 'VBF', 'ZH', 'WH', 'ttH']
        chain = ROOT.TChain('outtree')
        for y in year:
            for p in prod:
                chain.Add(direct +p+'_'+y+'_output.root')
        hist1_TH1_el = ROOT.TH1F('ggfSig1_th1f_el', 'ggfSig1_th1f_el', 340, 95, 180)
        hist2_TH1_el = ROOT.TH1F('ggfSig2_th1f_el', 'ggfSig2_th1f_el', 340, 95, 180)
        hist3_TH1_el = ROOT.TH1F('ggfSig3_th1f_el', 'ggfSig3_th1f_el', 340, 95, 180)
        hist4_TH1_el = ROOT.TH1F('ggfSig4_th1f_el', 'ggfSig4_th1f_el', 340, 95, 180)
        hist1_TH1_mu = ROOT.TH1F('ggfSig1_th1f_mu', 'ggfSig1_th1f_mu', 340, 95, 180)
        hist2_TH1_mu = ROOT.TH1F('ggfSig2_th1f_mu', 'ggfSig2_th1f_mu', 340, 95, 180)
        hist3_TH1_mu = ROOT.TH1F('ggfSig3_th1f_mu', 'ggfSig3_th1f_mu', 340, 95, 180)
        hist4_TH1_mu = ROOT.TH1F('ggfSig4_th1f_mu', 'ggfSig4_th1f_mu', 340, 95, 180)
        hist1_TH1 = ROOT.TH1F('ggfSig1_th1f', 'ggfSig1_th1f', 340, 95, 180)
        hist2_TH1 = ROOT.TH1F('ggfSig2_th1f', 'ggfSig2_th1f', 340, 95, 180)
        hist3_TH1 = ROOT.TH1F('ggfSig3_th1f', 'ggfSig3_th1f', 340, 95, 180)
        hist4_TH1 = ROOT.TH1F('ggfSig4_th1f', 'ggfSig4_th1f', 340, 95, 180)
        nsig1_el = 0
        nsig1_mu = 0
        nsig2_el = 0
        nsig2_mu = 0
        nsig3_el = 0
        nsig3_mu = 0
        nsig4_el = 0
        nsig4_mu = 0
        for entry in chain:
            if entry.BDT_score > bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist1_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist1_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig1_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist1_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist1_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig1_mu+= entry.weight_corr
            elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist2_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist2_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig2_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist2_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist2_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig2_mu+= entry.weight_corr
            elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist3_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist3_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig3_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist3_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist3_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig3_mu += entry.weight_corr
            elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist4_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist4_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig4_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist4_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist4_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig4_mu += entry.weight_corr
        self.ggf1El = ROOT.RooDataHist('hist_ggf1_sig_el', 'hist_ggf1_sig_el', x, hist1_TH1_el)
        self.ggf2El = ROOT.RooDataHist('hist_ggf2_sig_el', 'hist_ggf2_sig_el', x, hist2_TH1_el)
        self.ggf3El = ROOT.RooDataHist('hist_ggf3_sig_el', 'hist_ggf3_sig_el', x, hist3_TH1_el)
        self.ggf4El = ROOT.RooDataHist('hist_ggf4_sig_el', 'hist_ggf4_sig_el', x, hist4_TH1_el)
        self.ggf1Mu = ROOT.RooDataHist('hist_ggf1_sig_mu', 'hist_ggf1_sig_mu', x, hist1_TH1_mu)
        self.ggf2Mu = ROOT.RooDataHist('hist_ggf2_sig_mu', 'hist_ggf2_sig_mu', x, hist2_TH1_mu)
        self.ggf3Mu = ROOT.RooDataHist('hist_ggf3_sig_mu', 'hist_ggf3_sig_mu', x, hist3_TH1_mu)
        self.ggf4Mu = ROOT.RooDataHist('hist_ggf4_sig_mu', 'hist_ggf4_sig_mu', x, hist4_TH1_mu)
        self.ggf1 = ROOT.RooDataHist('hist_ggf1_sig', 'hist_ggf1_sig', x, hist1_TH1)
        self.ggf2 = ROOT.RooDataHist('hist_ggf2_sig', 'hist_ggf2_sig', x, hist2_TH1)
        self.ggf3 = ROOT.RooDataHist('hist_ggf3_sig', 'hist_ggf3_sig', x, hist3_TH1)
        self.ggf4 = ROOT.RooDataHist('hist_ggf4_sig', 'hist_ggf4_sig', x, hist4_TH1)
        self.nggf1_el = nsig1_el
        self.nggf1_mu = nsig1_mu
        self.nggf2_el = nsig2_el
        self.nggf2_mu = nsig2_mu
        self.nggf3_el = nsig3_el
        self.nggf3_mu = nsig3_mu
        self.nggf4_el = nsig4_el
        self.nggf4_mu = nsig4_mu
        
class readRuiROOTggFSignalggF:
    def __init__(self, x, direct='', year='', bdt1=0, bdt2=0, bdt3=0):
        chain = ROOT.TChain('outtree')
        chain.Add(direct + 'GGF_'+year+'_output.root')
        #chain.Add(direct + 'VBF_'+year+'_pinnacles_ggf_fixed.root')
        hist1_TH1_el = ROOT.TH1F('ggfSig1_th1f_ggfel', 'ggfSig1_th1f_ggfel', 340, 95, 180)
        hist2_TH1_el = ROOT.TH1F('ggfSig2_th1f_ggfel', 'ggfSig2_th1f_ggfel', 340, 95, 180)
        hist3_TH1_el = ROOT.TH1F('ggfSig3_th1f_ggfel', 'ggfSig3_th1f_ggfel', 340, 95, 180)
        hist4_TH1_el = ROOT.TH1F('ggfSig4_th1f_ggfel', 'ggfSig4_th1f_ggfel', 340, 95, 180)
        hist1_TH1_mu = ROOT.TH1F('ggfSig1_th1f_ggfmu', 'ggfSig1_th1f_ggfmu', 340, 95, 180)
        hist2_TH1_mu = ROOT.TH1F('ggfSig2_th1f_ggfmu', 'ggfSig2_th1f_ggfmu', 340, 95, 180)
        hist3_TH1_mu = ROOT.TH1F('ggfSig3_th1f_ggfmu', 'ggfSig3_th1f_ggfmu', 340, 95, 180)
        hist4_TH1_mu = ROOT.TH1F('ggfSig4_th1f_ggfmu', 'ggfSig4_th1f_ggfmu', 340, 95, 180)
        for entry in chain:
            if entry.BDT_score > bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist1_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist1_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist2_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist2_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist3_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist3_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist4_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist4_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
        self.ggf1El = ROOT.RooDataHist('hist_ggf1_sig_ggf_el', 'hist_ggf1_sig_ggf_el', x, hist1_TH1_el)
        self.ggf2El = ROOT.RooDataHist('hist_ggf2_sig_ggf_el', 'hist_ggf2_sig_ggf_el', x, hist2_TH1_el)
        self.ggf3El = ROOT.RooDataHist('hist_ggf3_sig_ggf_el', 'hist_ggf3_sig_ggf_el', x, hist3_TH1_el)
        self.ggf4El = ROOT.RooDataHist('hist_ggf4_sig_ggf_el', 'hist_ggf4_sig_ggf_el', x, hist4_TH1_el)
        self.ggf1Mu = ROOT.RooDataHist('hist_ggf1_sig_ggf_mu', 'hist_ggf1_sig_ggf_mu', x, hist1_TH1_mu)
        self.ggf2Mu = ROOT.RooDataHist('hist_ggf2_sig_ggf_mu', 'hist_ggf2_sig_ggf_mu', x, hist2_TH1_mu)
        self.ggf3Mu = ROOT.RooDataHist('hist_ggf3_sig_ggf_mu', 'hist_ggf3_sig_ggf_mu', x, hist3_TH1_mu)
        self.ggf4Mu = ROOT.RooDataHist('hist_ggf4_sig_ggf_mu', 'hist_ggf4_sig_ggf_mu', x, hist4_TH1_mu)

class readRuiROOTggFSignalVBF:
    def __init__(self, x, direct='', year='', bdt1=0, bdt2=0, bdt3=0):
        chain = ROOT.TChain('output')
        #chain.Add(direct + 'GGF_'+year+'_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'VBF_'+year+'_output.root')
        hist1_TH1_el = ROOT.TH1F('ggfSig1_th1f_vbf_el', 'ggfSig1_th1f_vbf_el', 340, 95, 180)
        hist2_TH1_el = ROOT.TH1F('ggfSig2_th1f_vbf_el', 'ggfSig2_th1f_vbf_el', 340, 95, 180)
        hist3_TH1_el = ROOT.TH1F('ggfSig3_th1f_vbf_el', 'ggfSig3_th1f_vbf_el', 340, 95, 180)
        hist4_TH1_el = ROOT.TH1F('ggfSig4_th1f_vbf_el', 'ggfSig4_th1f_vbf_el', 340, 95, 180)
        hist1_TH1_mu = ROOT.TH1F('ggfSig1_th1f_vbf_mu', 'ggfSig1_th1f_vbf_mu', 340, 95, 180)
        hist2_TH1_mu = ROOT.TH1F('ggfSig2_th1f_vbf_mu', 'ggfSig2_th1f_vbf_mu', 340, 95, 180)
        hist3_TH1_mu = ROOT.TH1F('ggfSig3_th1f_vbf_mu', 'ggfSig3_th1f_vbf_mu', 340, 95, 180)
        hist4_TH1_mu = ROOT.TH1F('ggfSig4_th1f_vbf_mu', 'ggfSig4_th1f_vbf_mu', 340, 95, 180)
        for entry in chain:
            if entry.BDT_score > bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist1_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist1_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist2_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist2_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist3_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist3_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.met < 90 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist4_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist4_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
        self.ggf1El = ROOT.RooDataHist('hist_ggf1_sig_vbf_el', 'hist_ggf1_sig_vbf_el', x, hist1_TH1_el)
        self.ggf2El = ROOT.RooDataHist('hist_ggf2_sig_vbf_el', 'hist_ggf2_sig_vbf_el', x, hist2_TH1_el)
        self.ggf3El = ROOT.RooDataHist('hist_ggf3_sig_vbf_el', 'hist_ggf3_sig_vbf_el', x, hist3_TH1_el)
        self.ggf4El = ROOT.RooDataHist('hist_ggf4_sig_vbf_el', 'hist_ggf4_sig_vbf_el', x, hist4_TH1_el)
        self.ggf1Mu = ROOT.RooDataHist('hist_ggf1_sig_vbf_mu', 'hist_ggf1_sig_vbf_mu', x, hist1_TH1_mu)
        self.ggf2Mu = ROOT.RooDataHist('hist_ggf2_sig_vbf_mu', 'hist_ggf2_sig_vbf_mu', x, hist2_TH1_mu)
        self.ggf3Mu = ROOT.RooDataHist('hist_ggf3_sig_vbf_mu', 'hist_ggf3_sig_vbf_mu', x, hist3_TH1_mu)
        self.ggf4Mu = ROOT.RooDataHist('hist_ggf4_sig_vbf_mu', 'hist_ggf4_sig_vbf_mu', x, hist4_TH1_mu)

class readRuiROOTVBFSignal:
    def __init__(self, x, direct='', bdt1=0, bdt2=0, bdt3=0):
        year = ['2016', '2016APV', '2017', '2018', '2022', '2022EE', '2023', '2023BPix']
        prod = [ 'ZH', 'WH', 'ttH']
        chain = ROOT.TChain('outtree')
        
        for y in year:
            for p in prod:
                chain.Add(direct +p+'_'+y+'_output.root')
        
        chain.Add(direct + 'sig.root')
        hist1_TH1_el = ROOT.TH1F('vbfSig1_th1f_el', 'vbfSig1_th1f_el', 340, 95, 180)
        hist2_TH1_el = ROOT.TH1F('vbfSig2_th1f_el', 'vbfSig2_th1f_el', 340, 95, 180)
        hist3_TH1_el = ROOT.TH1F('vbfSig3_th1f_el', 'vbfSig3_th1f_el', 340, 95, 180)
        hist4_TH1_el = ROOT.TH1F('vbfSig4_th1f_el', 'vbfSig4_th1f_el', 340, 95, 180)
        hist1_TH1_mu = ROOT.TH1F('vbfSig1_th1f_mu', 'vbfSig1_th1f_mu', 340, 95, 180)
        hist2_TH1_mu = ROOT.TH1F('vbfSig2_th1f_mu', 'vbfSig2_th1f_mu', 340, 95, 180)
        hist3_TH1_mu = ROOT.TH1F('vbfSig3_th1f_mu', 'vbfSig3_th1f_mu', 340, 95, 180)
        hist4_TH1_mu = ROOT.TH1F('vbfSig4_th1f_mu', 'vbfSig4_th1f_mu', 340, 95, 180)
        hist1_TH1 = ROOT.TH1F('vbfSig1_th1f', 'vbfSig1_th1f', 340, 95, 180)
        hist2_TH1 = ROOT.TH1F('vbfSig2_th1f', 'vbfSig2_th1f', 340, 95, 180)
        hist3_TH1 = ROOT.TH1F('vbfSig3_th1f', 'vbfSig3_th1f', 340, 95, 180)
        hist4_TH1 = ROOT.TH1F('vbfSig4_th1f', 'vbfSig4_th1f', 340, 95, 180)
        nsig1_el = 0
        nsig1_mu = 0
        nsig2_el = 0
        nsig2_mu = 0
        nsig3_el = 0
        nsig3_mu = 0
        nsig4_el = 0
        nsig4_mu = 0
        for entry in chain:
            if entry.BDT_score > bdt1 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist1_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist1_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig1_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist1_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist1_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig1_mu+= entry.weight_corr
            elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist2_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist2_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig2_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist2_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist2_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig2_mu+= entry.weight_corr
            elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist3_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist3_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig3_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist3_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist3_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig3_mu += entry.weight_corr
            elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11:
                    hist4_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist4_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig4_el += entry.weight_corr
                elif entry.ll_lepid == 13:
                    hist4_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    hist4_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
                    nsig4_mu += entry.weight_corr
        self.vbf1El = ROOT.RooDataHist('hist_vbf1_sig_el', 'hist_vbf1_sig_el', x, hist1_TH1_el)
        self.vbf2El = ROOT.RooDataHist('hist_vbf2_sig_el', 'hist_vbf2_sig_el', x, hist2_TH1_el)
        self.vbf3El = ROOT.RooDataHist('hist_vbf3_sig_el', 'hist_vbf3_sig_el', x, hist3_TH1_el)
        self.vbf4El = ROOT.RooDataHist('hist_vbf4_sig_el', 'hist_vbf4_sig_el', x, hist4_TH1_el)
        self.vbf1Mu = ROOT.RooDataHist('hist_vbf1_sig_mu', 'hist_vbf1_sig_mu', x, hist1_TH1_mu)
        self.vbf2Mu = ROOT.RooDataHist('hist_vbf2_sig_mu', 'hist_vbf2_sig_mu', x, hist2_TH1_mu)
        self.vbf3Mu = ROOT.RooDataHist('hist_vbf3_sig_mu', 'hist_vbf3_sig_mu', x, hist3_TH1_mu)
        self.vbf4Mu = ROOT.RooDataHist('hist_vbf4_sig_mu', 'hist_vbf4_sig_mu', x, hist4_TH1_mu)
        self.vbf1 = ROOT.RooDataHist('hist_vbf1_sig', 'hist_vbf1_sig', x, hist1_TH1)
        self.vbf2 = ROOT.RooDataHist('hist_vbf2_sig', 'hist_vbf2_sig', x, hist2_TH1)
        self.vbf3 = ROOT.RooDataHist('hist_vbf3_sig', 'hist_vbf3_sig', x, hist3_TH1)
        self.vbf4 = ROOT.RooDataHist('hist_vbf4_sig', 'hist_vbf4_sig', x, hist4_TH1)
        self.nvbf1_el = nsig1_el
        self.nvbf1_mu = nsig1_mu
        self.nvbf2_el = nsig2_el
        self.nvbf2_mu = nsig2_mu
        self.nvbf3_el = nsig3_el
        self.nvbf3_mu = nsig3_mu
        self.nvbf4_el = nsig4_el
        self.nvbf4_mu = nsig4_mu
                    
class readRuiROOTVBFSignalggF:
    def __init__(self, x, direct='', year='', bdt1=0, bdt2=0, bdt3=0):
        chain = ROOT.TChain('outtree')
        chain.Add(direct + 'GGF_'+year+'_output.root')
        #chain.Add(direct + 'VBF_'+year+'_pinnacles_ggf_fixed.root')
        hist1_TH1_el = ROOT.TH1F('vbfSig1_th1f_ggf_el', 'vbfSig1_th1f_ggf_el', 340, 95, 180)
        hist2_TH1_el = ROOT.TH1F('vbfSig2_th1f_ggf_el', 'vbfSig2_th1f_ggf_el', 340, 95, 180)
        hist3_TH1_el = ROOT.TH1F('vbfSig3_th1f_ggf_el', 'vbfSig3_th1f_ggf_el', 340, 95, 180)
        hist4_TH1_el = ROOT.TH1F('vbfSig4_th1f_ggf_el', 'vbfSig4_th1f_ggf_el', 340, 95, 180)
        hist1_TH1_mu = ROOT.TH1F('vbfSig1_th1f_ggf_mu', 'vbfSig1_th1f_ggf_mu', 340, 95, 180)
        hist2_TH1_mu = ROOT.TH1F('vbfSig2_th1f_ggf_mu', 'vbfSig2_th1f_ggf_mu', 340, 95, 180)
        hist3_TH1_mu = ROOT.TH1F('vbfSig3_th1f_ggf_mu', 'vbfSig3_th1f_ggf_mu', 340, 95, 180)
        hist4_TH1_mu = ROOT.TH1F('vbfSig4_th1f_ggf_mu', 'vbfSig4_th1f_ggf_mu', 340, 95, 180)
        for entry in chain:
            if entry.BDT_score > bdt1 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist1_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist1_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist2_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist2_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist3_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist3_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist4_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist4_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
        self.vbf1El = ROOT.RooDataHist('hist_vbf1_sig_ggf_el', 'hist_vbf1_sig_ggf_el', x, hist1_TH1_el)
        self.vbf2El = ROOT.RooDataHist('hist_vbf2_sig_ggf_el', 'hist_vbf2_sig_ggf_el', x, hist2_TH1_el)
        self.vbf3El = ROOT.RooDataHist('hist_vbf3_sig_ggf_el', 'hist_vbf3_sig_ggf_el', x, hist3_TH1_el)
        self.vbf4El = ROOT.RooDataHist('hist_vbf4_sig_ggf_el', 'hist_vbf4_sig_ggf_el', x, hist4_TH1_el)
        self.vbf1Mu = ROOT.RooDataHist('hist_vbf1_sig_ggf_mu', 'hist_vbf1_sig_ggf_mu', x, hist1_TH1_mu)
        self.vbf2Mu = ROOT.RooDataHist('hist_vbf2_sig_ggf_mu', 'hist_vbf2_sig_ggf_mu', x, hist2_TH1_mu)
        self.vbf3Mu = ROOT.RooDataHist('hist_vbf3_sig_ggf_mu', 'hist_vbf3_sig_ggf_mu', x, hist3_TH1_mu)
        self.vbf4Mu = ROOT.RooDataHist('hist_vbf4_sig_ggf_mu', 'hist_vbf4_sig_ggf_mu', x, hist4_TH1_mu)

class readRuiROOTVBFSignalVBF:
    def __init__(self, x, direct='', year='', bdt1=0, bdt2=0, bdt3=0):
        chain = ROOT.TChain('outtree')
        #chain.Add(direct + 'GGF_'+year+'_pinnacles_ggf_fixed.root')
        chain.Add(direct + 'VBF_'+year+'_output.root')
        chain.Add(direct + 'ZH_'+year+'_output.root')
        chain.Add(direct + 'WH_'+year+'_output.root')
        chain.Add(direct + 'ttH_'+year+'_output.root')
        hist1_TH1_el = ROOT.TH1F('vbfSig1_th1f_vbf_el', 'vbfSig1_th1f_vbf_el', 340, 95, 180)
        hist2_TH1_el = ROOT.TH1F('vbfSig2_th1f_vbf_el', 'vbfSig2_th1f_vbf_el', 340, 95, 180)
        hist3_TH1_el = ROOT.TH1F('vbfSig3_th1f_vbf_el', 'vbfSig3_th1f_vbf_el', 340, 95, 180)
        hist4_TH1_el = ROOT.TH1F('vbfSig4_th1f_vbf_el', 'vbfSig4_th1f_vbf_el', 340, 95, 180)
        hist1_TH1_mu = ROOT.TH1F('vbfSig1_th1f_vbf_mu', 'vbfSig1_th1f_vbf_mu', 340, 95, 180)
        hist2_TH1_mu = ROOT.TH1F('vbfSig2_th1f_vbf_mu', 'vbfSig2_th1f_vbf_mu', 340, 95, 180)
        hist3_TH1_mu = ROOT.TH1F('vbfSig3_th1f_vbf_mu', 'vbfSig3_th1f_vbf_mu', 340, 95, 180)
        hist4_TH1_mu = ROOT.TH1F('vbfSig4_th1f_vbf_mu', 'vbfSig4_th1f_vbf_mu', 340, 95, 180)
        for entry in chain:
            if entry.BDT_score > bdt1 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist1_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist1_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt2 and entry.BDT_score < bdt1 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist2_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist2_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > bdt3 and entry.BDT_score < bdt2 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist3_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist3_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif entry.BDT_score > -1 and entry.BDT_score < bdt3 and entry.weight_corr < 0.5:
                if entry.ll_lepid == 11: hist4_TH1_el.Fill(entry.llphoton_refit_m, entry.weight_corr)
                elif entry.ll_lepid == 13: hist4_TH1_mu.Fill(entry.llphoton_refit_m, entry.weight_corr)
        self.vbf1El = ROOT.RooDataHist('hist_vbf1_sig_vbf_el', 'hist_vbf1_sig_vbf_el', x, hist1_TH1_el)
        self.vbf2El = ROOT.RooDataHist('hist_vbf2_sig_vbf_el', 'hist_vbf2_sig_vbf_el', x, hist2_TH1_el)
        self.vbf3El = ROOT.RooDataHist('hist_vbf3_sig_vbf_el', 'hist_vbf3_sig_vbf_el', x, hist3_TH1_el)
        self.vbf4El = ROOT.RooDataHist('hist_vbf4_sig_vbf_el', 'hist_vbf4_sig_vbf_el', x, hist4_TH1_el)
        self.vbf1Mu = ROOT.RooDataHist('hist_vbf1_sig_vbf_mu', 'hist_vbf1_sig_vbf_mu', x, hist1_TH1_mu)
        self.vbf2Mu = ROOT.RooDataHist('hist_vbf2_sig_vbf_mu', 'hist_vbf2_sig_vbf_mu', x, hist2_TH1_mu)
        self.vbf3Mu = ROOT.RooDataHist('hist_vbf3_sig_vbf_mu', 'hist_vbf3_sig_vbf_mu', x, hist3_TH1_mu)
        self.vbf4Mu = ROOT.RooDataHist('hist_vbf4_sig_vbf_mu', 'hist_vbf4_sig_vbf_mu', x, hist4_TH1_mu)

class readRuiROOTSignal:
    def __init__(self, x, direct='', sampParam = ['ggf1','el','2016APV','ggf'], bdtBins=[0,0,0], rangeInfo=[340,95,180]):
        chain = ROOT.TChain('outtree')

        cat = sampParam[0]
        flav = sampParam[1]
        year = sampParam[2]
        prod = sampParam[3]

        bins = rangeInfo[0]
        lower = rangeInfo[1]
        upper = rangeInfo[2]

        metCut = 90
        if flav == 'el': lep = 11
        elif flav == 'mu': lep = 13

        if prod == 'ggf':
            chain = ROOT.TChain('outtree')
            chain.Add(direct + 'GGF_'+year+'_output.root')
        elif prod == 'vbf':
            chain = ROOT.TChain('outtree')
            chain.Add(direct + 'VBF_'+year+'_output.root')
            chain.Add(direct + 'ZH_'+year+'_output.root')
            chain.Add(direct + 'WH_'+year+'_output.root')
            chain.Add(direct + 'ttH_'+year+'_output.root')

        label = cat + '_th1f_' + prod + '_' + flav + '_' + year
        histlabel = 'hist' + label
        hist_TH1 = ROOT.TH1F(label, label, bins, lower, upper)
        for entry in chain:
            if entry.ll_lepid != lep:
                continue
            if cat == 'vbf1' and entry.BDT_score > bdtBins[0] and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'vbf2' and entry.BDT_score > bdtBins[1] and entry.BDT_score < bdtBins[0] and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'vbf3' and entry.BDT_score > bdtBins[2] and entry.BDT_score < bdtBins[1] and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'vbf4' and entry.BDT_score > -1 and entry.BDT_score < bdtBins[2] and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'ggf1' and entry.BDT_score > bdtBins[0] and entry.met < metCut and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'ggf2' and entry.BDT_score > bdtBins[1] and entry.BDT_score < bdtBins[0] and entry.met < metCut and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'ggf3' and entry.BDT_score > bdtBins[2] and entry.BDT_score < bdtBins[1] and entry.met < metCut and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)
            elif cat == 'ggf4' and entry.BDT_score > -1 and entry.BDT_score < bdtBins[2] and entry.met < metCut and entry.weight_corr < 0.5:
                hist_TH1.Fill(entry.llphoton_refit_m, entry.weight_corr)

        self.catflav = ROOT.RooDataHist(histlabel, histlabel, x, hist_TH1)



class readPico: #draw_pico datacard format
    def __init__(self, x, directory="", cat="", proc="", syst=""):
        """Initializes the following RooDataHist attributes for use with 
        HtoZg_fitting utilities: 
          data_hist_<CATNAME> 
        where <CAT> is ggf1, ggf2, ggf3, ggf4, vbf1, vbf2, vbf3, vh3l, vhmet,
        tthhad, or tthlep

        Args:
          x: three body invariant mass RooAbsReal
          directory: name of rawdata root file output by draw_pico
        
        """

        #TODO CAT_NAMES should be moved to some global location
        self.default_var = x
        root_file = ROOT.TFile(directory)
        weight = ROOT.RooRealVar("weight", "", -50.0, 50.0)
        is_data = (proc == "data_obs")
        dataname = "{}_cat_{}".format(proc, cat)
        if not is_data: dataname = f"mcdata_{dataname}_{syst}"
        print ("Read ",dataname)
        lowx = x.getMin()
        highx = x.getMax()
        mllg_cut = "mllg_cat_{0}>{1}&&mllg_cat_{0}<{2}".format(
                        cat, lowx, highx)
        ws_name = f"WS_{proc}_cat_{cat}"
        mllg_vars = getattr(root_file, ws_name).var(f"mllg_cat_{cat}")
        datasets_raw = getattr(root_file, ws_name).data(dataname).reduce(mllg_cut)
        setattr(self, dataname, self.convert_datahist_variable(dataname, datasets_raw, mllg_vars, x))
        '''
        CAT_NAMES = ["ggf1", "ggf2", "ggf3", "ggf4", "vbf1", "vbf2", "vbf3",
                     "vbf4", "vh3l", "vhmet", "tthhad", "tthlep"]
        PROCS = ["data_obs", "Htozg_el", "Htozg_mu", "Htomm"]
        SYST_NAMES = ["nominal", "CMS_scale_eUp", "CMS_scale_eDown", 
                      "CMS_res_eUp", "CMS_res_eDown", "CMS_scale_gUp", 
                      "CMS_scale_gDown", "CMS_res_gUp", "CMS_res_gDown", 
                      "CMS_scale_mUp", "CMS_scale_mDown"]

        self.default_var = x
        mllg_vars = []
        datasets_raw = []
        root_file = ROOT.TFile(directory)
        weight = ROOT.RooRealVar("weight", "", -50.0, 50.0)
        for proc in PROCS:
            for cat_name in CAT_NAMES:
                for syst in SYST_NAMES:
                    is_data = (proc == "data_obs")
                    dataname = "{}_cat_{}".format(proc, cat_name)
                    if is_data:
                        if syst != "nominal":
                            continue
                    else:
                        dataname = f"mcdata_{dataname}_{syst}"
                        print ("Read ", f"mcdata_{dataname}_{syst}")
                    #get RooDataSet from pico ROOT file and apply cut
                    #hack to apply range cut
                    lowx = x.getMin()
                    highx = x.getMax()
                    mllg_cut = "mllg_cat_{0}>{1}&&mllg_cat_{0}<{2}".format(
                        cat_name, lowx, highx)
                    ws_name = f"WS_{proc}_cat_{cat_name}"
                    mllg_vars.append(getattr(root_file, ws_name)
                        .var(f"mllg_cat_{cat_name}"))
                    datasets_raw.append(getattr(root_file, ws_name)
                        .data(dataname).reduce(mllg_cut))

                    #create RooDataHist and save as an attribute of this class
                    setattr(self, dataname, self.convert_datahist_variable(
                        dataname, datasets_raw[-1], mllg_vars[-1], x))
        '''
        
        
    def convert_datahist_variable(self, name, dataset, x_original, x):
        """Converts a RooDataSet into a RooDataHist with variable x

        Args:
          name: name of final histogram
          dataset: RooDataSet with arbitrary variable
          x_original: variable for original RooDataSet
          x: variable to use for final RooDataHist

        Returns: 
          RooDataHist representation of dataset with x as variable
        """
        lowx = x.getMin()
        highx = x.getMax()
        nbins = x.getBinning().numBins()
        dummy_hist = ROOT.TH1D("","",nbins,lowx,highx)
        dataset.fillHistogram(dummy_hist, ROOT.RooArgList(x_original))
        return ROOT.RooDataHist(name, name, x, dummy_hist)

    def numCheck(self):
        print("# bin1 bkg = ", self.hist_data_obs_cat_ggf1.sumEntries())
        print("# bin2 bkg = ", self.hist_data_obs_cat_ggf2.sumEntries())
        print("# bin3 bkg = ", self.hist_data_obs_cat_ggf3.sumEntries())
        print("# bin4 bkg = ", self.hist_data_obs_cat_ggf4.sumEntries())
