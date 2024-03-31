import DataFormats.FWLite as fwlite
import ROOT

events =		fwlite.Events("file:output.root")
secondaryVertices =	fwlite.Handle("std::vector<reco::VertexCompositeCandidate>")
primaryVertices = 	fwlite.Handle("std::vector<reco::Vertex>")

#make histogram of Kshort masses
ksmass_hist = ROOT.TH1D("ksmass_hist", "Kshort vertex masses; Mass [GeV]", 100, 0.4, 0.6)

events.toBegin()
for i, event in enumerate(events):
    print "Event:", i
    event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
    for j, vertex in enumerate(secondaryVertices.product()):
#   	 print "    Vertex:", j, vertex.vx(), vertex.vy(), vertex.vz()
#	 print "    Mass:", vertex.mass()
	 ksmass_hist.Fill(vertex.mass())
#    if i > 10: break

c = ROOT.TCanvas( "c", "c", 800, 800 )

ksmass_hist.Draw()

c.SaveAs("ksmasses.png")

#make histogram of Kshort SV-PV separation
ksdist_hist = ROOT.TH1D("ksdist_hist","3D Kshort SV-PV separations; Distance [cm]; N_{events}", 100, 0, 50)

for event in events:
	event.getByLabel("offlinePrimaryVertices", primaryVertices)
	pv = primaryVertices.product()[0] #get the first PV
	event.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondaryVertices)
	for vertex in secondaryVertices.product():
		dist = ( (pv.x() - vertex.vx())**2 + (pv.y() - vertex.vy())**2 + (pv.z() - vertex.vz())**2)**(0.5) #pythagoras
		ksdist_hist.Fill(dist)

c = ROOT.TCanvas( "c", "c", 800, 800 )

ksdist_hist.Draw()

c.SaveAs("ksdists.png")
