import ROOT
###############################
# The goal is compute the ratio
# of 2 significance/limits plots
# for stop & neutralino masses
###############################



###############################
#  TO-DO
#  - could become  parameters
#  - binning is hard-coded
###############################

ifilename_num = "sig.root"
ifilename_denom = "Ref/sig.root" 
ofilename = "ratio.root"


###############################
# Numerator
###############################
ifile_num = ROOT.TFile(ifilename_num,"OPEN")
h2_num = ifile_num.Get("sig");

###############################
# Denominator
###############################
ifile_denom = ROOT.TFile(ifilename_denom,"OPEN")
h2_denom = ifile_denom.Get("sig");

h2_ratio = h2_num.Clone()
h2_ratio.Divide(h2_denom)

c = ROOT.TCanvas()
#c.SetName("sig")
#c.cd()
h2_ratio.Draw("COLZ")

rfile = ROOT.TFile(ofilename,"RECREATE")
c.Write()
h2_ratio.Write()
rfile.Close()

