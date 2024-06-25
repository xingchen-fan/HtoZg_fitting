#include "RooRealVar.h"
#include "RooDataSet.h"
#include <cmath>
#include <sstream>
using namespace RooFit;

// Toy1 script, simple two samples: fitTo and chi2FitTo, two options, please choose as intended. Please remember chi2FitTo doesn't have asymptotic error option!
// ROOT VERSION MUST BE HIGHER THAN 6.28 !!!

template <typename T>
std::string to_string_with_precision(T a_value, int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

void weightsample_test_toy1(){
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    TRandom *generator = new TRandom();
    RooRealVar x("x", "x", -1., 1.);

    RooRealVar b("b", "b", -0.5, -20, 20);
    RooRealVar a("a", "a", 0.5, -5, 20.);
    RooRealVar c("c", "c", 0., -5, 5);
    RooRealVar d("d", "d", 0., -10, 10);
    RooRealVar e("e", "e", 0.5, -10., 10);
    RooRealVar w1("w1", "w1", 1);
    RooRealVar w2("w2", "w2", 2);
//Simple two samples------------------------
    RooGenericPdf line1("line1", "line1", "@1*@0 + 0.5", RooArgList(x, a));
    RooGenericPdf line2("line2", "line2", "@1*@0 + 0.5", RooArgList(x, b));

    RooGenericPdf fitline("fitline", "fitline", "@1*@0 + 0.5", RooArgList(x, c));
    
    TH1D *h_sumw2 = new TH1D("hsumw2", "hsumw2", 60, -3., 3.);
    TH1D *h_offsumw2 = new TH1D("hoffsumw2", "hoffsumw2", 60, -3., 3.);
    TH1D *h_asym = new TH1D("hasym", "hasym", 60, -3., 3.);
    TH1D *h_dummy = new TH1D("hdummy", "hdummy", 60, -3., 3.);
    
    TH1D *h_dummy_cen = new TH1D("hdummy_cen", "hdummy_cen", 60, -0.03, 0.03);
    TH1D *h_weight_cen = new TH1D("hweight_cen", "hweight_cen", 60, -0.03, 0.03);
    
    double err_sumw2 = 0.;
    double err_offsumw2 = 0.;
    double err_asym = 0.;
    double err_dummy = 0.;

    double sumw2 = 0.;
    double offsumw2 = 0.;
    double asym = 0.;
    double dummy = 0.;


    int N = 10000;
    for (int i(0); i < N; i++){
      int N1 = (int)generator->Gaus(10000, 100);
      int N2 = (int)generator->Gaus(5000, 71);
      int N3 = (int)generator->Gaus(20000, 141);
//Simple 2 samples--------------
        auto data1 = line1.generate(x, NumEvents(N1));
        data1->addColumn(w1);
        RooDataSet wdata1("wdata1", "wdata1", data1,  RooArgSet(x, w1), "", "w1");
        auto data2 = line2.generate(x, NumEvents(N2));
        data2->addColumn(w2);
        RooDataSet wdata2("wdata2", "wdata2", data2, RooArgSet(x, w2),"", "w2");
        wdata1.append(wdata2);

        x.setBins(2);
        auto hist = RooDataHist("hist", "hist", x, wdata1);
        auto dummy = fitline.generateBinned(x, NumEvents(N3));
        if (i==1) std::cout << "Sum = " << hist.sumEntries() << std::endl << std::endl;

        // auto res = fitline.fitTo(hist, Save(true), Strategy(0),Offset(false), SumW2Error(true), PrintLevel(-1));
        auto res = fitline.chi2FitTo(hist, Save(true), Strategy(0),Offset(false),  PrintLevel(-1), DataError(RooAbsData::SumW2));
        if (i == 1) res->Print("V");
        h_weight_cen->Fill(c.getVal());
        h_sumw2->Fill(c.getVal()/c.getError());
        err_sumw2 += c.getError();
        d.setVal(0);
        c.setVal(0);
        d.setError(0);
        // res = fitline.fitTo(hist, Save(true), Strategy(0),Offset(false), SumW2Error(false), PrintLevel(-1));
        res = fitline.chi2FitTo(hist, Save(true), Strategy(0),Offset(false), PrintLevel(-1), DataError(RooAbsData::Poisson));
        if (i == 1) res->Print("V");
        h_offsumw2->Fill(c.getVal()/c.getError());
        err_offsumw2 += c.getError();
        c.setVal(0);
        d.setVal(0);
        d.setError(0);
        res = fitline.fitTo(hist, Save(true), Strategy(0),Offset(false), AsymptoticError(true), PrintLevel(-1));
        if (i == 1) res->Print("V");
        h_asym->Fill(c.getVal()/c.getError());
        err_asym += c.getError();
        c.setVal(0);
        d.setVal(0);
        d.setError(0);
        // res = fitline.fitTo(*dummy, Save(true), Strategy(0),Offset(false), PrintLevel(-1));
        res = fitline.chi2FitTo(*dummy, Save(true), Strategy(0),Offset(false), PrintLevel(-1));
        if (i == 1) res->Print("V");
        h_dummy_cen->Fill(c.getVal());
        h_dummy->Fill(c.getVal()/c.getError());
        err_dummy += c.getError();
        c.setVal(0);
        d.setVal(0);
        d.setError(0);

        if (i == 3){
            TCanvas *c2 = new TCanvas("c2", "c2", 400, 400);
            c2->cd();
            RooPlot *xframe = x.frame(Title("Weight Test"));
            hist.plotOn(xframe);
            xframe->Draw();
            c2->Draw();
        }
    }

    err_sumw2 /= N;
    err_offsumw2 /= N;
    err_asym /= N;
    err_dummy /= N;

    std::cout << "SumW2 err = " << err_sumw2 << ", no SumW2 err = " << err_offsumw2 << ", asymptotic err = " << err_asym << ", unweighted err = " << err_dummy << std::endl;
    // RooPlot *xframe = x.frame(Title("Weight Test"));
    // hist.plotOn(xframe);
    // fitline.plotOn(xframe,LineColor(kRed));
    h_sumw2->SetLineColor(1);
    h_sumw2->SetLineWidth(2);
    h_offsumw2->SetLineColor(2);
    h_offsumw2->SetLineWidth(2);
    // h_asym->SetLineColor(3);
    // h_asym->SetLineWidth(2);
    h_dummy->SetLineColor(4);
    h_dummy->SetLineWidth(2);
    auto legend = new TLegend(0.1,0.6,0.4,0.9);
    legend->AddEntry(h_sumw2,("SumW2 RMS = " + to_string_with_precision(h_sumw2->GetRMS())).c_str(),"l");
    // legend->AddEntry(h_asym,("Asymptotic RMS = " + to_string_with_precision(h_asym->GetRMS())).c_str(),"l");
    legend->AddEntry(h_offsumw2,("#splitline{Poisson RMS = " + to_string_with_precision(h_offsumw2->GetRMS()) + "}{}").c_str(),"l");
    legend->AddEntry(h_dummy,("Unweighted RMS = " + to_string_with_precision(h_dummy->GetRMS())).c_str(),"l");
    legend->SetTextSize(0.02);

    auto legend1 = new TLegend(0.1,0.8,0.4,0.9);
    legend1->AddEntry(h_dummy_cen, ("Unweighted RMS = " + to_string_with_precision(h_dummy_cen->GetRMS(), 4)).c_str(), "l");
    legend1->AddEntry(h_weight_cen, ("Weighted RMS = " + to_string_with_precision(h_weight_cen->GetRMS(), 4)).c_str(), "l");
    h_dummy_cen->SetLineColor(4);
    h_dummy_cen->SetLineWidth(2);
    h_weight_cen->SetLineColor(2);
    h_weight_cen->SetLineWidth(2);


    // h_dummy->Fit("gaus");
    // h_asym->Fit("gaus");
    TCanvas *c0 = new TCanvas("c0", "c0", 800, 800);
    gStyle->SetOptStat(0);
    c0->cd();
    // h_asym->Draw();
    // h_asym->SetTitle("Pull");
    // h_asym->GetXaxis()->SetTitle("Pull (a/#sigma_{a})");
    // h_asym->GetXaxis()->SetTitleOffset(1.);
    h_sumw2->Draw();
    h_dummy->Draw("same");
    h_dummy->SetTitle("Pull");
    h_dummy->GetXaxis()->SetTitle("Pull (d/#sigma_{d})");
    h_dummy->GetXaxis()->SetTitleOffset(1.);
    
    h_offsumw2->Draw("same");
    legend->Draw("same");
    c0->Draw();

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();
    h_dummy_cen->SetTitle("Central Value");
    h_dummy_cen->Draw();
    h_weight_cen->Draw("same");
    legend1->Draw("same");
    c1->Draw();

}
