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

#===> Fit to a piecewise function: exponential decay past the peak and x^3 below it.

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

##Cubic part of function : [xMin, xMid]
def cube_func(x, params):
	res = params[0]*pow(x[0], 3)
	return res

cube = ROOT.TF1("cube", cube_func, xMin, (xMid + xMin)/2.0, 1)
#cube.SetParameters(1.0)
#ksbetagamma.Fit(cube, "SR", "", xMin, (xMid + xMin)/2.0)


##Planck's law-type function
def planck_func(x, params):
        a = params[0]
        b = params[1]
        c = 1.0#params[2]
        d = 0.0#params[3]
        k = 0.5

        res = (   a*pow(x[0]+k,3) / (  np.exp( b*(x[0] + k) ) - c  )   ) + d
        return res

planck = ROOT.TF1("planck", planck_func, xMin, xMid, 2)
#planck.SetParameters(67,1.0/5.35)
#ksbetagamma.Fit(planck, "SR", "")


##gaussian
peak = ksbetagamma.GetBinContent(binMax)
def gauss_func(x,params):
	A = peak#params[0]
	B = params[0]
	k = xMid#params[2] #x-coord of center
	res = A*np.exp(-1.0*B*(x[0] - k)**2)
	return res

gauss = ROOT.TF1("gauss", gauss_func, xMin, xMid, 1)
#ksbetagamma.Fit(gauss, "SR", "")
#gauss.Draw("same")

##Exponential part of function : (xMid, xMax]
def exp_dec_func(x, params):
	lifetime = params[0]
	xshift = params[1]
	res = np.exp(-(x[0] + xshift)/lifetime)
	return res

expMin = 3.0
exp_dec = ROOT.TF1("exp_dec", exp_dec_func, expMin, xMax, 2)
exp_dec.SetParameters(1.0, 0.0)
ksbetagamma.Fit(exp_dec, "SR", "", expMin, xMax)

fitted_exp = ksbetagamma.GetFunction("exp_dec")
lifetime = fitted_exp.GetParameter(0)
print "lifetime : ",lifetime

##product of gaussian and exponential decay
peak = ksbetagamma.GetBinContent(binMax)
def gauxp_func(x,params):
	H = params[0] #peak height
	a = params[1] #1/(2*sigma^2)
	b = params[2] #reciprocal of (expoential's) mean
	k = params[3] #center of peak
	expArg = ( -1.0*a*(x[0] - k)**2 ) - ( b*x[0] )
	res = H*np.exp( expArg )
	return res

gauxp = ROOT.TF1("gauxp", gauxp_func, xMin, xMax, 4)
gauxp.SetParameters(peak, 1999.8, lifetime, xMid) #2500.0, 2.5
ksbetagamma.Fit(gauxp, "SR", "")
gauxp.Draw("same")




# # # # # # # # # # # # # # # # # # # # # #


c1 = ROOT.TCanvas( "c1", "c1", 800, 800)
ksbetagamma.Draw()
c1.SaveAs("ksbetagammas.png")
