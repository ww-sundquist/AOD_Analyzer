import DataFormats.FWLite as fwlite
import ROOT

events =		fwlite.Events("file:output.root")
secondaryVertices =	fwlite.Handle("std::vector<reco::VertexCompositeCandidate>")
primaryVertices = 	fwlite.Handle("std::vector<reco::Vertex>")

#make histogram of Kshort masses
ksmass_hist = ROOT.TH1D("ksmass_hist", "Kshort vertex masses; Mass [GeV]", 100, 0.4, 0.6)

events.toBegin()
for i, event in enumerate(events):
#    print "Event:", i
    event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
    for j, vertex in enumerate(secondaryVertices.product()):
#   	 print "    Vertex:", j, vertex.vx(), vertex.vy(), vertex.vz()
#	 print "    Mass:", vertex.mass()
	 ksmass_hist.Fill(vertex.mass())
#    if i > 10: break

c = ROOT.TCanvas( "c", "c", 800, 800 )


# Fits # # # # # # # # # # # # # # # # # # # # # # # #

xMin = 0.48#restrict domain of fit, values taken by eye
xMax = 0.52
bottom = 0.0
top = 4300.0

# naive gaussian
fitMin = ROOT.TLine(xMin, bottom, xMin, top)
fitMax = ROOT.TLine(xMax, bottom, xMax, top)
fitMin.SetLineColor(ROOT.kTeal)
fitMin.SetLineStyle(2)
fitMax.SetLineColor(ROOT.kTeal)
fitMax.SetLineStyle(2)
#gausFit = ksmass_hist.Fit("gaus","SR", "", xMin, xMax)


#Product of two gaussian

# Parameter naming convention
# [0] : mu_1, mean of 1st gaussian
# [1] : sigma_1, std. dev. of 1st gaussian
# [2] : mu_2, mean of 2nd gaussian
# [3] : sigma_2, std. dev. of 2nd gaussian

def gausProd_func(x, params):
	pi = 3.1415926535897932385
	e = 2.7182818284590452356
	expArg = (-1.0/2)*( ( (x[0] - params[0])/(params[1]) )**2 + ( (x[0] - params[2])/(params[3]) )**2 )
	prefac = 1.0/(2*pi*params[1]*params[3])
	res = prefac * (e**expArg)
	return res

# Initialize parameters
init_mu_1 = 0.0
init_sigma_1 = 1.0
init_mu_2 = 0.0
init_sigma_2 = 1.0

gausProd = ROOT.TF1("gausProd", gausProd_func, xMin, xMax, 4)
gausProd.SetParameters(init_mu_1, init_sigma_1, init_mu_2, init_sigma_2)

ksmass_hist.Fit(gausProd, "SR", "", xMin, xMax)

# # # # # # # # # # # # # # # # # # # # # # # # # # #


ksmass_hist.Draw()
fitMin.Draw()
fitMax.Draw()

#legend = ROOT.TLegend()
#legend.Draw()

c.SaveAs("ksmasses.png")

#make histogram of Kshort SV-PV separation
#ksdist_hist = ROOT.TH1D("ksdist_hist","3D Kshort SV-PV separations; Distance [cm]; N_{events}", 100, 0, 50)
#
#events.toBegin()
#for event in events:
#	event.getByLabel("offlinePrimaryVertices", primaryVertices)
#	pv = primaryVertices.product()[0] #get the first PV
#	event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
#	for vertex in secondaryVertices.product():
#		dist = ( (pv.x() - vertex.vx())**2 + (pv.y() - vertex.vy())**2 + (pv.z() - vertex.vz())**2)**(0.5) #pythagoras
#		ksdist_hist.Fill(dist)
#
#c2 = ROOT.TCanvas( "c2", "c2", 1000, 800 )
#
#ksdist_hist.Draw()
#
#c2.SaveAs("ksdists.png")
