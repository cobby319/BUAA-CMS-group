import ROOT as r
import os
import plottery as ply

histos= {'pT_JPsi_tot':'pT_JPsi_final','DeltaR_JpsiPi1_':'DeltaR_JpsiPi1_final','DeltaR_JpsiPi2_':'DeltaR_JpsiPi2_final','JPiPi_lxy_':'JPiPi_lxy_final'}
file = r.TFile("Merged.root")
for hist in histos:
    h1= r.TH1F()
    h2= r.TH1F()
    file.GetObject("histos/"+hist,h1)
    file.GetObject("histos/"+histos[hist],h2)
    ply.plot_hist(
        bgs=[h1],
        legend_labels = ["before cut"],
        sigs = [h2],
        sig_labels = ["signal mass region"],
        options = {
            "legend_scalex": 0.7,
            "legend_scaley": 1.5,
            "cms_label": "Preliminary",
            "lumi_value": "-inf",
            "do_stack": False,
            "yaxis_log": True,
            "output_name": histos[hist]+".pdf",
            "yaxis_range": [0.1,1000000000],
            "yaxis_moreloglabels":False,
            #"yaxis_noexponents":True,
            #"output_ic": True,
            }
        )
