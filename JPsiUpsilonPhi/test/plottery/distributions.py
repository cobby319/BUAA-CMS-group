import ROOT as r
import os
import plottery as ply

histos= ['pT_JPsi_','DeltaR_JpsiPi1_','DeltaR_JpsiPi2_','JPiPi_lxy_','M_JPsi','"Deta_JpsiPi1','Deta_JpsiPi2','Dphi_JpsiPi1','"Dphi_JpsiPi2','mumNHits','mumNPHits','mupNHits','mupNPHits','mu1soft','mu2soft','mu1tight','mu2tight','mu1loose','mu2loose','J_vertexchi2','Pi_nhits1','Pi_npixelhits1','Pi_nhits2','Pi_npixelhits2','Pi_vertexchisq1','Pi_vertexchisq2']
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
