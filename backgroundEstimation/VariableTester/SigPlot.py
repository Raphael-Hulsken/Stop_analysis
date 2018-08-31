import ROOT
###############################
# The goal is to read signifance
# for stop & neutralino masses
# in a txt file and produce 
# a plot
###############################



###############################
#  TO-DO
#  - could become  parameters
#  - binning is hard-coded
###############################

ofilename = "significance.txt"
rootfilename = "sig.root"

plot = ROOT.TH2D("sig","",200,200,1200,200,0,1000)

# read file and fill the plot
for line in open(ofilename):
  values = line.split()
  print values
  plot.SetBinContent(plot.GetXaxis().FindBin(float(values[0])), plot.GetYaxis().FindBin(float(values[1])), float(values[2]))


c = ROOT.TCanvas()
#c.SetName("sig")
#c.cd()
plot.Draw("COLZ")

rfile = ROOT.TFile(rootfilename,"RECREATE")
c.Write()
plot.Write()
rfile.Close()

