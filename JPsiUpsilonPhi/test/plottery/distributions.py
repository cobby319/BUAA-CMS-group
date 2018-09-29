import ROOT as r
import os
import plottery as ply

histos= {'pT_JPsi_':'J_pT','DeltaR_JpsiPi1':'JP1_deltaR_','DeltaR_JpsiPi2':'JP2_deltaR_','JPiPi_lxy_':'JPP_lxy'}
file = r.TFile("../Merged.root")
for hist in histos:
    h= r.TH1F()
    file.GetObject("histos/"+hist,h)
    ply.plot_hist(
        bgs=[h],
        options = {
            "do_stack": False,
            "yaxis_log": True,
            "output_name": "../distributions.pdf",
            "output_ic": True,
            }
        )