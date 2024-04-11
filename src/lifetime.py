import DataFormats.FWLite as fwlite
import ROOT
import math
import numpy as np

events =		fwlite.Events("file:output.root")
secondaryVertices =	fwlite.Handle("std::vector<reco::VertexCompositeCandidate>")
primaryVertices = 	fwlite.Handle("std::vector<reco::Vertex>")

#make historgram of \beta\gamma values
ksbetagamma = ROOT.TH1D("ksbetagamma", "K-shorts ; #beta#gamma", 416, 0, 52)

events.toBegin()
for event in events:
	event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
	for vertex in secondaryVertices.product():
		px = vertex.px()
		three_momentum = (   ( vertex.px() )**2 + ( vertex.py() )**2 + ( vertex.pz() )**2   )**(0.5)
		mass = vertex.mass()
		betagamma = three_momentum / mass
		ksbetagamma.Fill(betagamma)


# # # # # # # # # # # # # # # # # # # # # #

#===> Fit to a piecewise function: exponential decay past the peak and gaussian below it.

fitCan = ROOT.TCanvas("fitCan","Fitting #beta#gamma",800,800)
##Start of function abscissa
def find_first_nonzero_bin(hist):
	n_bins = hist.GetNbinsX()
	for bin in range(1, n_bins + 1):
		if hist.GetBinContent(bin) > 0:
			return bin
	return None

leftBin = find_first_nonzero_bin(ksbetagamma)
xMin = ksbetagamma.GetXaxis().GetBinCenter(leftBin)
print "x coordinate of leftmost histogram bin  :  ",xMin

##End of function abscissa
def find_last_nonzero_bin(hist):
	n_bins = hist.GetNbinsX()
	for bin in range(n_bins, 0, -1):
		if hist.GetBinContent(bin) > 0:
			return bin
	return None

rightBin = find_last_nonzero_bin(ksbetagamma)
xMax = ksbetagamma.GetXaxis().GetBinCenter(rightBin)
print "x coordinate of rightmost histogram bin :  ",xMax

##Piecewise changeover abscissa
binMax = ksbetagamma.GetMaximumBin()
xMid = ksbetagamma.GetXaxis().GetBinCenter(binMax)
print "x coordinate of histogram max value     :  ",xMid


# # # Piecewise by two fits

##Gaussian part of function : [xMin, xMid]
peak = ksbetagamma.GetBinContent(binMax)
def gauss_func(x,params):
	A = peak#params[0]
	B = params[0]
	k = xMid#params[2] #x-coord of center
	res = A*np.exp(-1.0*B*(x[0] - k)**2)
	return res

gauss = ROOT.TF1("gauss", gauss_func, xMin, xMid, 1)
gausfit = ksbetagamma.Fit(gauss, "SR", "")
gauss.Draw("same")

##Exponential part of function : [xMid+0.1, xMax]
def exp_dec_func(x, params):
	lifetime = params[0]
	xshift = params[1]
	res = np.exp(-(x[0] + xshift)/lifetime)
	return res

exp_dec = ROOT.TF1("exp_dec", exp_dec_func, xMid+0.1, xMax, 2)
exp_dec.SetParameters(1.0, 0.0)
expfit = ksbetagamma.Fit(exp_dec, "SR", "")
exp_dec.Draw("same")

fitCan.Draw()

fitted_exp = ksbetagamma.GetFunction("exp_dec")
lifetime = fitted_exp.GetParameter(0)
print "lifetime : ",lifetime


# # # Piecewise by fit to piecewise function
pieceCan = ROOT.TCanvas("pieceCan", "title", 800, 800) #add title

def piecewise_func(x, params):
	

piecewise = ROOT.TF1("piecewise", piecewis_func, xMin, xMas, ) #add no. params
piecefit = ksbetagamma.Fit(piecewise, "SR", "")
piecewise.Draw()




# # # # # # # # # # # # # # # # # # # # # #


c1 = ROOT.TCanvas( "c1", "c1", 800, 800)
ksbetagamma.Draw()
c1.SaveAs("ksbetagammas.png")
