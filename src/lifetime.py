import DataFormats.FWLite as fwlite
import ROOT
import math

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

xMin = 0.0#restrict domain of fit, values taken by eye
xMax = 40.0

#Fit beta gamma dist'n to f(x) = (A exp(-x) x) - (B exp(-x^2) x) - (C exp(-x) x^2) - (D exp(-x^2) x^2) - (E exp(-x^2) x^3) + (F exp(-x) x^4)

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

ksbetagamma.Fit(banerjee, "SR", "", xMin, xMax)

# # # # # # # # # # # # # # # # # # # # # #


c1 = ROOT.TCanvas( "c1", "c1", 800, 800)
ksbetagamma.Draw()
c1.SaveAs("ksbetagammas.png")
