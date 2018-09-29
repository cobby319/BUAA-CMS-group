import ROOT as r
import os
import plottery as ply

histos= ['Pi_pt1','Pi_pt2','Pi_eta1','Pi_eta2','J_mass','pT_JPsi_','DeltaR_JpsiPi1_','DeltaR_JpsiPi2_','JPiPi_lxy_','Deta_JpsiPi1_','Deta_JpsiPi2_','Dphi_JpsiPi1_','Dphi_JpsiPi2_','mumNHits_','mumNPHits_','mupNHits_','mupNPHits_','mu1soft_','mu2soft_','mu1tight_','mu2tight_','mu1loose_','mu2loose_','J_vertexchi2_','Pi_nhits1_','Pi_npixelhits1_','Pi_nhits2_','Pi_npixelhits2_','Pi_vertexchisq1_','Pi_vertexchisq2_']
file = r.TFile("Merged.root")
if not os.path.exists('distribution_PLOTS'): 
    os.makedirs('distribution_PLOTS')
for hist in histos:
    h1= r.TH1F()
    h2= r.TH1F()
    file.GetObject("histos/"+hist,h1)
    file.GetObject("histos/"+hist+"final",h2)
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
            "output_name": "distribution_PLOTS/test0929_"+hist+".pdf",
            "yaxis_range": [0.1,1000000000],
            "yaxis_moreloglabels":False,
            #"yaxis_noexponents":True,
            #"output_ic": True,
            }
        )
