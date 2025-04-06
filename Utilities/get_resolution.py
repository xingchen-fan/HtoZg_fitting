
import ROOT
from argparse import ArgumentParser

CATEGORIES = ["ggf1", "ggf2", "ggf3", "ggf4", "vbf1", "vbf2", "vbf3",
              "vbf4", "vh3l", "vhmet", "tthhad", "tthlep"]
PROCS = ["Htozg_el", "Htozg_mu"]

def parse_args():
    """Parses arguments

    Returns:
      namespace with parsed arguments
    """
    parser = ArgumentParser(description = "Get resolutions")
    parser.add_argument('-i', '--input_file', help = 'rawdata root filename')
    args = parser.parse_args()
    return args

def get_mean_and_resolution(hist: ROOT.TH1D, name: str):
    """Finds the mean and the 68.27% interval around the mean from a histogram,
    and prints the result

    Args:
      hist: data to perform calculation on
      name: name of category to print out
    """
    norm = hist.Integral(0,hist.GetNbinsX()+1)
    percentile = 0.0
    mean = None
    res_lo = None
    res_hi = None
    for ibin in range(hist.GetNbinsX()+2):
        #print(f"bin: {ibin}, percentile: {percentile}")
        bin_yield = hist.GetBinContent(ibin)
        next_percentile = percentile + bin_yield/norm
        if percentile < 0.15865 and next_percentile > 0.15865:
            res_lo = hist.GetXaxis().GetBinCenter(ibin)
            #print("res_lo")
        if percentile < 0.5 and next_percentile > 0.5:
            mean = hist.GetXaxis().GetBinCenter(ibin)
            #print("mean")
        if percentile < 0.84135 and next_percentile > 0.84135:
            res_hi = hist.GetXaxis().GetBinCenter(ibin)
            #print("res_hi")
        percentile = next_percentile
    res = (res_hi-res_lo)/2.0 #/2?
    print(f"Category {name}: {mean}, {res}")

def get_resolutions(input_filename: str):
    """Prints resolutions for signal distributions from a given root file

    Args:
      input_filename: name of rawdata root file
    """
    input_file = ROOT.TFile(input_filename)
    for category in CATEGORIES:
        for proc in PROCS:
            var_name = f"mllg_cat_{category}"
            ws_name = f"WS_{proc}_cat_{category}"
            dataset_name = f"mcdata_{proc}_cat_{category}_nominal"
            dataset = getattr(input_file, ws_name).data(dataset_name)
            var = getattr(input_file, ws_name).var(var_name)
            #0.01 GeV precision
            hist = ROOT.TH1D(f"hist_{proc}_{category}","",1000,120.0,130.0)
            dataset.fillHistogram(hist, ROOT.RooArgList(var))
            #f"{var_name}>120.0&&{var_name}<130.0")
            get_mean_and_resolution(hist, ws_name)

def main():
    """Prints resolutions for signal distributions from a given root file
    """
    args = parse_args()
    get_resolutions(args.input_file)

if __name__=="__main__":
    main()
