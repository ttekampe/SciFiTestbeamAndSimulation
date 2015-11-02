import sys, argparse

class plot:
	name = ""
	lowerBoundary = 0
	upperBoundary = 0
	nBins = 1
	xAxisLabel = 0
	yAxisLabel = 0
	def binWidth(self):
		return str((self.upperBoundary-self.lowerBoundary)/float(self.nBins))


#parser = argparse.ArgumentParser(description='Plot cluster properties from data and simulation.')
#parser.add_argument('-d', '--data', type=str)
#parser.add_argument('-s', '--simulation')
#parser.add_argument("-p", "--position")

#cfg = parser.parse_args()


def reweightHist(hist, factor):
	for i in xrange(hist.GetNbinsX()):
		hist.SetBinContent(i, hist.GetBinContent(i) * factor)

import lhcbStyle
from ROOT import gStyle, TFile, TTree, TH1D, TCanvas, kRed, kOpenCircle, gPad
#ROOT.PyConfig.IgnoreCommandLineOptions = True
def doPlots(f, f_sim, name):
	plots = []
	cwm = plot()
	cwm.name = "chargeWeightedMean"
	cwm.lowerBoundary = 0
	cwm.upperBoundary = 90
	cwm.nBins = 180
	cwm.xAxisLabel = "Charge weighted mean"
	cwm.yAxisLabel = "Events / " + cwm.binWidth()
	#plots.append(cwm)

	hwm = plot()
	hwm.name = "hitWeightedMean"
	hwm.lowerBoundary = 0
	hwm.upperBoundary = 90
	hwm.nBins = 180
	hwm.xAxisLabel = "Hit weighted mean"
	hwm.yAxisLabel = "Events / " + hwm.binWidth()
	#plots.append(hwm)

	sc = plot()
	sc.name = "sumCharge"
	sc.lowerBoundary = 0
	sc.upperBoundary = 80
	sc.nBins = 80
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
	#can.SaveAs("clusterPlots_Pos" + cfg.position + ".pdf[")
	can.SaveAs(name + ".pdf[")
	for p in plots:
		hist = TH1D(p.name, p.name, p.nBins, p.lowerBoundary, p.upperBoundary)
		hist_sim = TH1D(p.name + "_sim", p.name + "_sim", p.nBins, p.lowerBoundary, p.upperBoundary)
		#if p.name=="sumCharge":
		#	reweightHist(hist_sim, 1.2)

		#distance_from_track = 0 +/- 400

		t.Draw(p.name + ">>" + p.name, "distance_from_track<300&&distance_from_track>-300", "")
		t_sim.Draw(p.name + ">>" + p.name + "_sim")


		hist.GetXaxis().SetTitle(p.xAxisLabel)
		hist.GetYaxis().SetTitle(p.yAxisLabel)
		hist.SetStats(True)

		hist_sim.GetXaxis().SetTitle(p.xAxisLabel)
		hist_sim.GetYaxis().SetTitle(p.yAxisLabel)
		hist_sim.SetStats(True)

		hist_sim.SetMarkerColor(kRed)
		hist_sim.SetLineColor(kRed)
		hist_sim.SetMarkerStyle(kOpenCircle)

		hist.SetName("Data")
		hist_sim.SetName("Simulation")

		hist.Draw("e")
		can.Update()
		statbox = hist.GetListOfFunctions().FindObject('stats')

		#"""Moves the statbox to a new location.  Uses the NDC coordinates where the canvas is 1 wide and tall, (0,0) is the bottom left."""

		if not p.name == "sumCharge":
			statbox.SetY1NDC(0.4)
			statbox.SetY2NDC(0.6)
		else:
			statbox.SetX1NDC(0.5)
			statbox.SetX2NDC(0.7)

		hist_sim.Draw("e")
		can.Update()
		statbox_sim = hist_sim.GetListOfFunctions().FindObject('stats')
		statbox_sim.SetTextColor(kRed)

		#if hist.GetMaximum() > hist_sim.GetMaximum():
		if hist.GetMaximum()/hist.GetEntries() > hist_sim.GetMaximum()/hist_sim.GetEntries():
			hist.DrawNormalized("e")
			gPad.Modified()
			hist_sim.DrawNormalized("esames")
		else:
			hist_sim.DrawNormalized("e")
			gPad.Modified()
			hist.DrawNormalized("esames")


		#hist.Draw("esames")
		#can.SaveAs("clusterPlots_Pos" + cfg.position + ".pdf")
		can.SaveAs(name + ".pdf")

	#hist = f.Get("clusterShape")
	#hist_sim = f_sim.Get("clusterShape")
	#
	#hist.GetXaxis().SetRangeUser(0, 10)
	#hist_sim.GetXaxis().SetRangeUser(0, 10)
	#
	#hist.GetYaxis().SetRangeUser(0, 50)
	#hist_sim.GetYaxis().SetRangeUser(0, 50)
	#
	#hist.GetXaxis().SetTitle("Channel relative to seed")
	#hist.GetYaxis().SetTitle("Adc response")
	#hist.SetName("Data")
	#hist.SetStats(True)
	#
	#hist_sim.GetXaxis().SetTitle("Channel relative to seed")
	#hist_sim.GetYaxis().SetTitle("Adc response")
	#hist_sim.SetName("Simulation")
	#hist_sim.SetStats(True)
	#
	#hist_sim.SetMarkerColor(kRed)
	#hist_sim.SetLineColor(kRed)
	#hist_sim.SetMarkerStyle(kOpenCircle)
	#
	#
	#
	#hist.Draw("e")
	#can.Update()
	#statbox = hist.GetListOfFunctions().FindObject('stats')
	##statbox.SetLabel("Data")
	##hist_sim.Draw("e")
	##statbox_sim = hist_sim.GetListOfFunctions().FindObject('stats')
	##statbox_sim.SetName("Simulation")
	#
	#statbox.SetX1NDC(0.2)
	#statbox.SetX2NDC(0.4)
	##statbox.SetY1NDC(y1)
	##statbox.SetY2NDC(y2)
	#
	#hist_sim.Draw("e")
	#can.Update()
	#statbox_sim = hist_sim.GetListOfFunctions().FindObject('stats')
	#statbox_sim.SetTextColor(kRed)
	#
	#
	#
	#
	#if hist.GetMaximum() > hist_sim.GetMaximum():
	#	hist.Draw("e")
	#	gPad.Modified()
	#	hist_sim.Draw("esames")
	#else:
	#	hist_sim.Draw("e")
	#	gPad.Modified()
	#	hist.Draw("esames")
	#
	#
	#
	#
	#can.SaveAs("clusterPlots_Pos" + cfg.position + ".pdf")
	can.SaveAs(name + ".pdf]")

lhcbStyle.setLHCbStyle()
gStyle.SetOptStat(1)

#measurements = [
#		 ["1431786652", "a_at_0deg", "pos_a_0_deg"]
#		,["1432091294", "a_at_10deg", "pos_a_10_deg"]
#		,["1432264510", "a_at_20deg", "pos_a_20_deg"]
#		,["1432350729", "a_at_30deg", "pos_a_30_deg"]
#		,["1432169457", "c_at_0deg", "pos_c_0_deg"]
#		,["1432089102", "c_at_10deg", "pos_c_10_deg"]
#		,["1432187205", "c_at_20deg", "pos_c_20_deg"]
#		,["1432358809", "c_at_30deg", "pos_c_30_deg"]
#		]


#for measurement in measurements:
#
#	f = TFile("/home/tobi/SciFi/results/clusters/btsoftware_" + measurement[0] + "_datarun_ntuple_corrected_clusterAnalyis.root", "READ")
#	if not f.IsOpen():
#		print("Could not open data file for run number " + measurement[0])
#		continue
#	t = f.Get("clusterAnalysis")
#
#	f_sim = TFile("/home/tobi/SciFi/results/clusters/testbeam_simulation_position_" + measurement[1] + "_clusterAnalyisnoise.root", "READ")
#	if not f_sim.IsOpen():
#		print("Could not open sim file for configuration " + measurement[1])
#		continue
#	t_sim = f_sim.Get("clusterAnalysis")
#
#	doPlots(f, f_sim, measurement[2])

f = TFile("/home/tobi/SciFi/results/clusters/btsoftware_1432169457_datarun_ntuple_corrected_clusterAnalyis.root", "READ")
if not f.IsOpen():
	print("Could not open data file for run number")

t = f.Get("clusterAnalysis")

if not t:
	print("tree not found!")

f_sim = TFile("/home/tobi/SciFi/results/clusters/testbeam_simulation_position_c_at_0degnewLight_clusterAnalyis.root", "READ")
if not f_sim.IsOpen():
	print("Could not open sim file for configuration")

t_sim = f_sim.Get("clusterAnalysis")

if not t_sim:
	print("sim tree not found!")

doPlots(f, f_sim, "clusters_pos_C_0deg_newLight")


f = TFile("/home/tobi/SciFi/results/clusters/btsoftware_1431786652_datarun_ntuple_corrected_clusterAnalyis.root", "READ")
if not f.IsOpen():
	print("Could not open data file for run number")

t = f.Get("clusterAnalysis")
if not t:
	print("tree not found!")

f_sim = TFile("/home/tobi/SciFi/results/clusters/testbeam_simulation_position_a_at_0degnewLight_clusterAnalyis.root", "READ")
if not f_sim.IsOpen():
	print("Could not open sim file for configuration")

t_sim = f_sim.Get("clusterAnalysis")
if not t_sim:
	print("sim tree not found!")

doPlots(f, f_sim, "clusters_pos_a_0deg_newLight")



#f = TFile("/home/tobi/SciFi/results/clusters/btsoftware_1432169457_datarun_ntuple_corrected_clusterAnalyis.root", "READ")
#if not f.IsOpen():
#	print("Could not open data file for run number")
#
#t = f.Get("clusterAnalysis")
#
#f_sim = TFile("/home/tobi/SciFi/results/clusters/simulationResponse_fibMatVolCor_newTags_PosC_clusterAnalyis.root", "READ")
#if not f_sim.IsOpen():
#	print("Could not open sim file for configuration")
#
#t_sim = f_sim.Get("clusterAnalysis")
#
#doPlots(f, f_sim, "pos_C_0deg_defAtt")
