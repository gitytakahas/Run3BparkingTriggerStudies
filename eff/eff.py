from ROOT import TFile, TH2F, TCanvas, gStyle, gROOT, kTRUE
import copy, math
gROOT.SetBatch(kTRUE)

def addOverflows(his):
    if his.GetSumw2N() == 0 : his.Sumw2(kTRUE)
    xbins = his.GetXaxis().GetNbins()
    ybins = his.GetYaxis().GetNbins()
    entries = his.GetEntries()
    for xbin in range(1,xbins+1):
        content = his.GetBinContent(xbin,ybins)
        outlier = his.GetBinContent(xbin,ybins+1)
        his.SetBinContent(xbin,ybins,content+outlier)
        his.SetBinContent(xbin,ybins+1,0.)
    content = his.GetBinContent(xbin,ybins)
    outlier = his.GetBinContent(xbins+1,ybins+1)
    his.SetBinContent(xbins,ybins,content+outlier)
    his.SetBinContent(xbins+1,ybins+1,0.)
    his.SetEntries(entries)
    return his

def createPdf(his,canvas,eff=True):
    canvas.cd()
    his.Draw("TEXT COLZ")
    his.SetStats(0)
    if eff: 
        gStyle.SetPaintTextFormat("4.2f");
        his.SetMinimum(1.e-6)
        his.SetMaximum(1.)
        his.SetMarkerSize(1.6)
    else: 
        gStyle.SetPaintTextFormat(".0f");
        his.SetMinimum(0.1)
        his.SetMarkerSize(1.)
    canvas.SetLogz()
    return canvas

# numer
print "numer"
files_available = 79.
files_processed = 77.
numer_file = TFile('root/numer.root')
numer_histo1 = numer_file.Get('numer_gen_lead')
numer_histo2 = numer_file.Get('numer_gen_sub')
numer_histo = numer_histo1.Clone("numer")
numer_histo.Add(numer_histo2)
numer_histo = addOverflows(numer_histo)
numer_histo.SetTitle("numer")
numer_histo.Scale(files_available/files_processed)

# denom
print "denom"
DAS_entries = 31585558.
denom_file = TFile('root/denom.root')
denom_histo = denom_file.Get('denom')
denom_histo = addOverflows(denom_histo)
denom_histo2 = denom_histo.Clone("denom_unweighted")
denom_histo.Scale(DAS_entries/denom_histo.GetEntries())

# eff
print "eff"
eff_histo = numer_histo.Clone("eff")
eff_histo.Divide(numer_histo,denom_histo)
eff_histo.SetTitle("Analysis efficiency")
#print "minimum  eff: ", eff_histo.GetMinimum(1.e-9)

# output
eff_file = TFile('root/eff.root','RECREATE')
numer_histo1.Write()
numer_histo2.Write()
numer_histo.Write()
denom_histo.Write()
eff_histo.Write()

numer_canvas = TCanvas()
numer_canvas = createPdf(numer_histo,numer_canvas,eff=False)
numer_canvas.SaveAs("plots/numer.pdf")

numer_canvas1 = TCanvas()
numer_canvas1 = createPdf(numer_histo1,numer_canvas1,eff=False)
numer_canvas1.SaveAs("plots/numer1.pdf")

numer_canvas2 = TCanvas()
numer_canvas2 = createPdf(numer_histo2,numer_canvas2,eff=False)
numer_canvas2.SaveAs("plots/numer2.pdf")

denom_canvas = TCanvas()
denom_canvas = createPdf(denom_histo,denom_canvas,eff=False)
denom_canvas.SaveAs("plots/denom.pdf")

denom_canvas2 = TCanvas()
denom_canvas2 = createPdf(denom_histo2,denom_canvas2,eff=False)
denom_canvas2.SaveAs("plots/denom_unweighted.pdf")

eff_canvas = TCanvas()
eff_canvas = createPdf(eff_histo,eff_canvas)
eff_canvas.SaveAs("plots/eff.pdf")

eff_file.Close()
