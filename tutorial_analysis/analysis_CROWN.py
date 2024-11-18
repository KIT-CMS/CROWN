import ROOT
 
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
 
# Create dataframe from Ntuple files
path="../build/bin/out_mm.root"
df = ROOT.RDataFrame("ntuple", path)
 
# Make histogram of dimuon mass spectrum. Note how we can set titles and axis labels in one go.
h = df.Histo1D(("Muon_InvMass", "Dimuon mass;m_{#mu#mu} (GeV);N_{Events}", 30000, 0.25, 300), "Muon_InvMass")
  
# Produce plot
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetTextFont(42)
c = ROOT.TCanvas("c", "", 800, 700)
c.SetLogx(); c.SetLogy()
 
h.SetTitle("")
h.GetXaxis().SetTitleSize(0.04)
h.GetYaxis().SetTitleSize(0.04)
h.Draw()
 
label = ROOT.TLatex(); label.SetNDC(True)
label.DrawLatex(0.175, 0.740, "#eta")
label.DrawLatex(0.205, 0.775, "#rho,#omega")
label.DrawLatex(0.270, 0.740, "#phi")
label.DrawLatex(0.400, 0.800, "J/#psi")
label.DrawLatex(0.415, 0.670, "#psi'")
label.DrawLatex(0.485, 0.700, "Y(1,2,3S)")
label.DrawLatex(0.755, 0.680, "Z")
label.SetTextSize(0.040); label.DrawLatex(0.100, 0.920, "#bf{CMS Open Data}")
label.SetTextSize(0.030); label.DrawLatex(0.630, 0.920, "#sqrt{s} = 8 TeV, L_{int} = 11.6 fb^{-1}")
 
c.SaveAs("dimuon_spectrum_CROWN.pdf")
