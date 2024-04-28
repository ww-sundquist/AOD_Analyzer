import DataFormats.FWLite as fwlite
import ROOT
import math
import numpy as np

events =		fwlite.Events("file:output.root")
secondaryVertices =	fwlite.Handle("std::vector<reco::VertexCompositeCandidate>")
primaryVertices = 	fwlite.Handle("std::vector<reco::Vertex>")

#make historgram of \beta\gamma values
ksbetagamma = ROOT.TH1D("ksbetagamma", "K-shorts ; #beta#gamma", 416, 0, 52)

ksNum = 0 #number of k-shorts
events.toBegin()
for event in events:
	event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
	for vertex in secondaryVertices.product():
		px = vertex.px()
		three_momentum = (   ( vertex.px() )**2 + ( vertex.py() )**2 + ( vertex.pz() )**2   )**(0.5)
		mass = vertex.mass()
		betagamma = three_momentum / mass
		#betagamma = betagamma * 2.998*10.**8. #unitless
		ksbetagamma.Fill(betagamma)
		ksNum = ksNum + 1


print "ksNum :  ",ksNum

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
#print "x coordinate of leftmost histogram bin  :  ",xMin

##End of function abscissa
def find_last_nonzero_bin(hist):
	n_bins = hist.GetNbinsX()
	for bin in range(n_bins, 0, -1):
		if hist.GetBinContent(bin) > 0:
			return bin
	return None

rightBin = find_last_nonzero_bin(ksbetagamma)
xMax = ksbetagamma.GetXaxis().GetBinCenter(rightBin)
#print "x coordinate of rightmost histogram bin :  ",xMax

##Piecewise changeover abscissa
binMax = ksbetagamma.GetMaximumBin()
xMid = ksbetagamma.GetXaxis().GetBinCenter(binMax)
#print "x coordinate of histogram max value     :  ",xMid


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
#gausfit = ksbetagamma.Fit(gauss, "SR", "")
#gauss.Draw("same")

##Exponential part of function : [xMid+0.1, xMax]
def exp_dec_func(x, params):
	lifetime = params[0]
	xshift = params[1]
	res = np.exp(-(x[0] + xshift)/lifetime)
	return res

exp_dec = ROOT.TF1("exp_dec", exp_dec_func, xMid+0.1, xMax, 2)
exp_dec.SetParameters(1.0, 0.0)
#expfit = ksbetagamma.Fit(exp_dec, "SR", "")
#exp_dec.Draw("same")

#fitCan.Draw()

#fitted_exp = ksbetagamma.GetFunction("exp_dec")
#lifetime = fitted_exp.GetParameter(0)
#print "lifetime : ",lifetime


# # # Piecewise by fit to piecewise function
pieceCan = ROOT.TCanvas("pieceCan", "title", 800, 800) #add title

def piecewise_func(x, params):
	x = x[0]
	H = params[0] #height of central peak at xMid
	# Gausian parameters
	w = params[1] #narrowness of gaussian : w = 1/(2sigma^2)
	k = xMid #center of gaussian
	# Exponential parameters
	L = params[2] #mean val of exp decay
	if x < k:
		y = H*np.exp(  -1.0*w*( x - k )**2  )
	if x > k:
		y = ( H*np.exp(xMid / L) )*np.exp(-1.0*x/L)
	if x == k:
		y = H
	return y

piecewise = ROOT.TF1("piecewise", piecewise_func, xMin, xMax, 3) #add no. params
piecewise.SetParameters(1.0, 1.0, 1.0)
piecefit = ksbetagamma.Fit(piecewise, "SR", "")
piecewise.Draw()

pieceCan.Draw()
pieceCan.SaveAs("ksbgfit.png")

# # # # # # # # # # # # # # # # # # # # # #

# Now that we have our pdf, we can generate a collection of random numbers using it

#start by generating 1,000 random numbers on the range startX to endX, quantized every stepsize
startX = 0.0
endX = 2.0 #40.0
stepsize = 0.01 #1.0
samplesize = 1000#ksNum

fitted_piece = ksbetagamma.GetFunction("piecewise")
param0 = fitted_piece.GetParameter(0)
param1 = fitted_piece.GetParameter(1)
param2 = fitted_piece.GetParameter(2)
param3 = fitted_piece.GetParameter(3)
params = [param0, param1, param2, param3]


#xVals = np.arange(startX+1,endX+1, stepsize) #list of x-values
xVals = np.linspace(startX, endX, samplesize) #list of x-values
xProbs = [] #list of assoc. probabilities
xProbs_unnormed = []
for i in xVals:
	res = piecewise_func([i], params)
	xProbs_unnormed.append(res)

#normalize probabilies
for i in xProbs_unnormed:
	res = i / sum(xProbs_unnormed)
	xProbs.append(res)

genNums	= [] #generated numbers
nBins = (endX-startX)/stepsize
#genNums_hist = ROOT.TH1D("genNums_hist", "#beta#gamma", int(round(nBins)), startX, endX)
#for i in range(samplesize):
#	res = np.random.choice(xVals, p=xProbs)
#	genNums.append(res)
#	genNums_hist.Fill(res)

#print "generated numbers  : ",genNums


# # # # # # # # # # # # # # # # # # # # # #

# Test case:  I want to generate 3 distributions
##	1. For samplesize number of instances, pick a number on [0,1] with equal probability, take the natural log and multiply by c\tau
##	2. For samplesize number of instances, pick a number on [startX,endX] with probability density function given by my piecewise_fun()
##	3. For samplesize number of instances, (A) pick a number on [startX,endX] with pdf given by my piecewise_fun(); (B) pick a number on [0,1] with equal probability, take the natural log and multiply by c\tau; and (C) multiply the results of A and B together

bg_hist = ROOT.TH1D("bg_hist", "#beta#gamma distribution", int(round(nBins)), startX, endX) #beta gamma
ln_hist = ROOT.TH1D("ln_hist", "Exponential decay distribution", int(round(nBins)), startX, endX) #exponential decay
genDist_hist = ROOT.TH1D("genDist_hist", "Generated K-short lifetime disribution", int(round(nBins)), startX, endX) #generated k-short lifetimes

#take a particular tau
tau = 10*8.95**(-11.)#12.0 #seconds

for i in range(samplesize):

	bg = np.random.choice(xVals, p=xProbs) #\beta\gamma value generated from pdf of best fit
	bg_hist.Fill(bg)

	c = 2.998*(10.0**8.0) #meters/second
	ln = -1.0 * c * tau * np.log( np.random.choice(np.linspace( 0.0, 1.0, samplesize, endpoint=True )) )
	ln_hist.Fill(ln)
	#print "ln : ",ln

	gen = bg*ln
	genDist_hist.Fill(gen)
	#print "Gen : ",gen

#get 'experimentat' \bega\gamma c\tau
ksdist_hist = ROOT.TH1D("ksdist_hist","3D Kshort SV-PV separations; Distance [cm]; N_{events}", int(round(nBins)), startX, endX)

events.toBegin()
counter = 0
for event in events:
	if counter < 1000:
		event.getByLabel("offlinePrimaryVertices", primaryVertices)
		pv = primaryVertices.product()[0] #get the first PV
		event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
		for vertex in secondaryVertices.product():
			dist = ( (pv.x() - vertex.vx())**2 + (pv.y() - vertex.vy())**2 + (pv.z() - vertex.vz())**2)**(0.5) #pythagoras
			#ksdist_hist.Fill(dist)
			ksdist_hist.Fill(dist/100.) #try converting cm to m
			counter = counter + 1
	else:
		break

# Print gen distribution alongside experiment
c_compare = ROOT.TCanvas( "c_compare", "c_compare", 1000, 1000 )
genDist_hist.SetLineColor(ROOT.kBlue)
genDist_hist.Draw("hist same")
ksdist_hist.SetLineColor(ROOT.kRed)
ksdist_hist.Draw("hist same")
c_compare.Draw()
c_compare.SaveAs("gen_vs_exp_lifetimedist.png")

chi2 = genDist_hist.Chi2Test(ksdist_hist , "UU" )
print "-------\nChi-squared value:  ",chi2

# # # # # # # # # # # # # # # # # # # # # #

# reproduce above test case in functions that can be called iteratively to generate chi2-dist'n

def generateDist(tau, samplesize, startX, endX):
	genDist_hist = ROOT.TH1D("genDist_hist", "Generated K-short lifetime disribution", int(round(nBins)), startX, endX) #generated k-short lifetimes
	
	for count in range(samplesize):

		bg = np.random.choice(xVals, p=xProbs) #\beta\gamma value generated from pdf of best fit

		c = 2.998*(10.0**8.0) #meters/second
		ln = -1.0 * c * tau * np.log( np.random.choice(np.linspace( 0.0, 1.0, samplesize, endpoint=True )) )
	
		gen = bg*ln
		genDist_hist.Fill(gen)

	return genDist_hist

def getChi2(tau, ksdist_hist, samplesize, startX, endX):
	genDist_hist = generateDist(tau, samplesize, startX, endX)
	chi2 = genDist_hist.Chi2Test( ksdist_hist , "UU" )
	return chi2

# - - - - - - - - - - - - - - - - - - - - #

#establish the range of lifetime values over which we are searching
n = 1000 #number of points
#tau_range = np.linspace( 1.*(10.**(-11.)), 1.*(10.**(-9)), n, endpoint=True )
tau_range = np.linspace( 0.0, 10.0*(10.**(-11)), n, endpoint=True )

chi2_vals = [] #array of y-coords
for lifetime in tau_range:
	y = getChi2(tau, ksdist_hist, samplesize, startX, endX)
	chi2_vals.append(y)

chi2_dist = ROOT.TGraph( n )
for i in range(n):
	chi2_dist.SetPoint( i , tau_range[i], chi2_vals[i] )

c_chi2dist = ROOT.TCanvas( "c_chi2dist", "c_chi2dist", 1200, 800 )
chi2_dist.SetMarkerStyle(21)
chi2_dist.SetMarkerSize(1)
chi2_dist.Draw("AP")
chi2_dist.SetTitle("Comparing generated and experimental distributions;Lifetime (s);#chi-square")
c_chi2dist.Draw()
c_chi2dist.SaveAs("chi2dist.png")

# Questions:
#	-- Why does the gen dist lifetime by inspection seem to be around a factor of 10 larger than expected?
#	-- What is throwing all the "There is a bing in h1 with less than 1 event." errors in the Chi2Test?
#	-- Why does the chi2 distribution look like a mess rather than following a trend?

# # # # # # # # # # # # # # # # # # # # # #

# Take a look at the exp dec dist
#c_ed = ROOT.TCanvas( "c_ed", "c_ed", 800, 800 )
#ln_hist.SetFillStyle(3005)
#ln_hist.Draw("hist")
#c_ed.Draw()

# # # # # # # # # # # # # # # # # # # # # #


# Do the same with our exp. dec. pdf
#
#intermediate = []
#for i in range(samplesize):
#	res = np.random.choice(np.linspace( 0.0, 1.0, samplesize, endpoint=True ))
#	intermediate.append(res)

#def getExpNums(tau, dist, startX, endX, nBins):
#	expNums_hist = ROOT.TH1D("expNums_hist", "Exponential decay", int(round(nBins)), startX, endX)
#	c = 2.998*(10**8) #speed of light
#	for i in dist:
#		res = -1.0 * tau * c * np.log(i)
#		#print "i : ",i,"  | res : ",res
#		expNums_hist.Fill(res)
#	return expNums_hist

#take a particular tau
#tau = 10.0**(-10.)#12.0
#expNums_hist = getExpNums(tau, intermediate, startX, endX, nBins)

#genLifetimes = genNums_hist * expNums_hist

#scale genLifetimes
#genLifetimes = genLifetimes.Scale(0.001)

#draw hists
#c2 = ROOT.TCanvas("c2","c2", 1000, 800)

#genNums_hist.SetLineColor(ROOT.kRed)
#genNums_hist.SetFillColor(ROOT.kRed)
#genNums_hist.SetFillStyle(3003)
#genNums_hist.Draw("same")

#expNums_hist.SetLineColor(ROOT.kBlue)
#expNums_hist.SetFillColor(ROOT.kBlue)
#expNums_hist.SetFillStyle(3005)
#expNums_hist.Draw("same")

#genLifetimes.SetLineColor(ROOT.kBlack)
#genLifetimes.SetFillColor(ROOT.kBlack)
#genLifetimes.SetFillStyle(3003)
#genLifetimes.SetTitle("Product (scaled by 1/1000)")
#genLifetimes.Draw("hist same")

#legend2 = c2.BuildLegend(0.7,0.7,0.9,0.9,"Distributions")
#genNums_hist.SetTitle("Generating lifetime distributions; #beta#gamma c#tau (1/MeV); N_{events}")
#c2.SetTitle("Generating lifetime distributions")
#c2.Draw()
#c2.SaveAs("genLifetimes.png")

#unscale genLifetimes
#genLifetimes.Scale(1000)

#get 'experimentat' \bega\gamma c\tau
#ksdist_hist = ROOT.TH1D("ksdist_hist","3D Kshort SV-PV separations; Distance [cm]; N_{events}", 100, 0, 50)

#events.toBegin()
#for event in events:
#	event.getByLabel("offlinePrimaryVertices", primaryVertices)
#	pv = primaryVertices.product()[0] #get the first PV
#	event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
#	for vertex in secondaryVertices.product():
#		dist = ( (pv.x() - vertex.vx())**2 + (pv.y() - vertex.vy())**2 + (pv.z() - vertex.vz())**2)**(0.5) #pythagoras
#		#ksdist_hist.Fill(dist)
#		ksdist_hist.Fill(dist/100.) #try converting cm to m

#c_troubleshoot = ROOT.TCanvas( "c_troubleshoot", "c_troubleshoot", 900, 900 )
#ksdist_hist.SetLineColor(ROOT.kGreen)
#ksdist_hist.Draw("hist same")
#ksbetagamma.SetLineColor(ROOT.kPink)
#ksbetagamma.Draw("hist same")
#c_troubleshoot.Draw()
#c_troubleshoot.SaveAs("troubleshoot.png")

#def normalize(hist):
#	hist.Scale(1./hist.Integral(0,100,"width"))
#	return hist

#exp_norm = normalize(ksdist_hist)

#gen_norm = normalize(genLifetimes)


#c = ROOT.TCanvas( "c", "c", 800, 800 )

#ksdist_hist.Draw("hist same")
#ksdist_hist.SetLineColor(ROOT.kRed)
#ksdist_hist.SetFillStyle(3003)

#exp_norm.SetLineColor(ROOT.kRed)
#exp_norm.Draw("hist same")

#genLifetimes.Draw("hist same")
#genLifetimes.SetLineColor(ROOT.kBlue)
#genLifetimes.SetFillStyle(3003)

#gen_norm.SetLineColor(ROOT.kBlue)
#gen_norm.Draw("hist same")

#ksdist_hist.SetTitle("Generated and experimental lifetimes; #beta #gamma c #tau (1/MeV); N_{events}")
#c.SaveAs("exp_and_gen.png")


#chi2 = genLifetimes.Chi2Test(ksdist_hist , "UU" )
#rint "-------\nChi-squared value:  ",chi2
#print "-------\nexp_norm area:  ",exp_norm.Integral(),"    gen_norm area:  ",gen_norm.Integral()

#to do:
# - find chi-square comparison between getLifetimes and experiment
# - loop over different tau values and optimize







# # # # # # # # # # # # # # # # # # # # # #

c1 = ROOT.TCanvas( "c1", "c1", 800, 800)
ksbetagamma.Draw()
c1.SaveAs("ksbetagammas.png")
