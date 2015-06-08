from ROOT import gStyle, TFile, TTree, TH1D, TCanvas
import lhcbStyle
import sys

class plot:
	name = ""
	lowerBoundary = 0
	upperBoundary = 0
	nBins = 1
	xAxisLabel = 0
	yAxisLabel = 0
	def binWidth(self):
		return str((self.upperBoundary-self.lowerBoundary)/float(self.nBins))

lhcbStyle.setLHCbStyle()
gStyle.SetOptStat(1)

filename = sys.argv[1]

f = TFile(filename, "READ")
t = f.Get("clusterAnalysis")

plots = []
cwm = plot()
cwm.name = "chargeWeightedMean"
cwm.lowerBoundary = 0
cwm.upperBoundary = 90
cwm.nBins = 50
cwm.xAxisLabel = "Charge weighted mean"
cwm.yAxisLabel = "Events / " + cwm.binWidth()
plots.append(cwm)

hwm = plot()
hwm.name = "hitWeightedMean"
hwm.lowerBoundary = 0
hwm.upperBoundary = 90
hwm.nBins = 50
hwm.xAxisLabel = "Hit weighted mean"
hwm.yAxisLabel = "Events / " + hwm.binWidth()
plots.append(hwm)

sc = plot()
sc.name = "sumCharge"
sc.lowerBoundary = 0
sc.upperBoundary = 100
sc.nBins = 50
sc.xAxisLabel = "Cluster charge"
sc.yAxisLabel = "Events / " + sc.binWidth()
plots.append(sc)

mc = plot()
mc.name = "maxCharge"
mc.lowerBoundary = 0
mc.upperBoundary = 40
mc.nBins = 40
mc.xAxisLabel = "Seed charge"
mc.yAxisLabel = "Events / " + mc.binWidth()
plots.append(mc)

cs = plot()
cs.name = "clusterSize"
cs.lowerBoundary = 0
cs.upperBoundary = 10
cs.nBins = 10
cs.xAxisLabel = "Cluster Size"
cs.yAxisLabel = "Events / " + cs.binWidth()
plots.append(cs)

can = TCanvas()
can.SaveAs("clusterPlots.pdf[")
for p in plots:
	hist = TH1D("hist", "", p.nBins, p.lowerBoundary, p.upperBoundary)
	t.Draw(p.name + ">>hist")
	hist.GetXaxis().SetTitle(p.xAxisLabel)
	hist.GetYaxis().SetTitle(p.yAxisLabel)
	hist.SetStats(True)
	hist.Draw()
	can.SaveAs("clusterPlots.pdf")
can.SaveAs("clusterPlots.pdf]")