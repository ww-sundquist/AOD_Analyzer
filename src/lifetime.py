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

xMin = 0.1#restrict domain of fit, values taken by eye
xMax = 50.0

###Fit beta gamma dist'n to f(x) = (A exp(-x) x) - (B exp(-x^2) x) - (C exp(-x) x^2) - (D exp(-x^2) x^2) - (E exp(-x^2) x^3) + (F exp(-x) x^4)

#Parameter definitions:
# A : [0]
# B : [1]
# C : [2]
# D : [3]
# E : [4]
# F : [5]

def banerjee_func(x, params):
	a = params[0]*math.exp(-1.0*x[0])*x[0]
	b = params[1]*math.exp(-1.0*(x[0])**2)*x[0]
	c = params[2]*math.exp(-1.0*x[0])*(x[0])**2
	d = params[3]*math.exp(-1.0*(x[0])**2)*(x[0])**2
	e = params[4]*math.exp(-1.0*(x[0])**2)*pow(x[0],3)
	f = params[5]*math.exp(-1.0*x[0])*pow(x[0],4)
	res = a - b - c + d - e + f
	return res

#Initialize parameters
#init_A = 1.0
#init_B = 1.0
#init_C = 1.0
#init_D = 1.0
#init_E = 1.0
#init_F = 1.0

banerjee = ROOT.TF1("banerjee", banerjee_func, xMin, xMax, 6)
#banerjee.SetParameters(init_A, init_B, init_C, init_D, init_E, init_F)

#ksbetagamma.Fit(banerjee, "SR", "", xMin, xMax)


###Fit beta gamms dist'n to skewed gaussian

def skewed_normal_func(x, params):
	loc = params[0]
	scale = params[1]
	skew = params[2]
	
	arg1 = (skew / (scale * math.sqrt(2)))*(x[0] - loc) #argument of erf
	arg2 = (-1.0 / (2 * scale**2))*(x[0] - loc)**2 #argument of exp
	coeff = 1.0 / (scale * math.sqrt(2*math.pi))
	res = coeff * (1 - math.erf(arg1)) * math.exp(arg2)
	return res

skewed_normal = ROOT.TF1("skewed_normal", skewed_normal_func, xMin, xMax, 4)
#skewed_normal.SetParameters(1.0, 1.0, 1.0, 1.0)

#ksbetagamma.Fit(skewed_normal, "SR", "", xMin, xMax)


###Fit to exponential decay past peak, then maybe use to initialize planck's law?

def exp_dec_func(x, params):
	lifetime = params[0]
	xshift = params[1]
	res = np.exp(-(x[0] + xshift)/lifetime)
	return res

expMin = 5.0
exp_dec = ROOT.TF1("exp_dec", exp_dec_func, expMin, xMax, 2)
exp_dec.SetParameters(1.0, 0.0)
#ksbetagamma.Fit(exp_dec, "SR", "", expMin, xMax)


###Fit to x^3 below peak, init'ze?
cubeMax = 3.0
def cube_func(x, params):
	res = params[0]*pow(x[0], 3)
	return res

cube = ROOT.TF1("cube", cube_func, xMin, cubeMax, 1)
#cube.SetParameters(1.0)
#ksbetagamma.Fit(cube, "SR", "", xMin, cubeMax)


###Planck's law-type function
def planck_func(x, params):
	a = params[0]
	b = params[1]
	c = 1.0#params[2]
	d = 0.0#params[3]
	k = 0.5

	res = (   a*pow(x[0]+k,3) / (  np.exp( b*(x[0] + k) ) - c  )   ) + d
	return res

planck = ROOT.TF1("planck", planck_func, xMin, xMax, 2)
planck.SetParameters(2.5*66.0,0.5/5.3)
ksbetagamma.Fit(planck, "SR", "", xMin, xMax)


# # # # # # # # # # # # # # # # # # # # # #


c1 = ROOT.TCanvas( "c1", "c1", 800, 800)
ksbetagamma.Draw()
c1.SaveAs("ksbetagammas.png")
