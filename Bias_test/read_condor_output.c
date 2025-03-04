#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

void read_condor_output(TString func, TString cat){
  auto h = new TH1D("h", "h", 80, -4, 4);
  int ncovered = 0;
  for (int i(0); i < 200; i++){
    TString number = to_string(i);
    ifstream file("condor/"+func+"_" + cat + "/output"+number+".txt");
    string line;
    if (file.is_open()) {
      while (getline(file, line)) {
	string pull = line.substr(0,4);
	string cover = line.substr(0,7);
	if (cover == "covered"){
	  ncovered += stoi(line.substr(11,1));
	}
	if (pull == "pull"){
	  int brac = line.find("[");
	  int comma = line.find(",");
	  string a_pull = line.substr(brac+1, comma-brac-1);
	  //cout << a_pull <<endl;
	  h->Fill(stod(a_pull));
	  bool theEnd = false;
	  int left = 0;
	  while(!theEnd){
	    left = comma;
	    comma = line.find(",", left+1);
	    a_pull = line.substr(left+2, comma-left-2);
	    //cout << a_pull <<endl;
	    h->Fill(stod(a_pull));
	    if (line.find(",", comma+1) == string::npos) theEnd = true;
	  }
	  a_pull = line.substr(comma+2, line.find("]") - comma - 2);
	  //cout << a_pull <<endl;
	  h->Fill(stod(a_pull));
	}
      }
      file.close();
    }
  }
  auto latex = new TLatex();
  latex->SetTextSize(0.04);
  auto c = new TCanvas("c", "c", 500, 500);
  h->Fit("gaus");
  h->GetXaxis()->SetTitle("Pull");
  h->SetTitle(cat + " " + func +" Pull");
  gStyle->SetOptFit(1);
  latex->DrawLatexNDC(0.15, 0.8, ("Coverage = " + to_string_with_precision(ncovered/h->GetEntries(), 2)).c_str());
  c->SaveAs("plots/"+func+"_" + cat+ "_pull.pdf");
}
