import ROOT as r
import os
import plottery as ply

histos= {'Pi_pt1_':'pT of track1','Pi_pt2_':'pT of track2','Pi_eta1_':'#eta of track1','Pi_eta2_':'#eta of track2','J_mass_':'dimuon invariant mass','pT_JPsi_':'dimuon pT','DeltaR_JpsiPi1_':'#DeltaR of J/#psi & track1','DeltaR_JpsiPi2_':'#DeltaR of J/#psi & track2','JPiPi_lxy_':'transverse path of the final reconstructed vertex','Deta_JpsiPi1_':'#Delta#eta of J/#psi & track1','Deta_JpsiPi2_':'#Delta#eta of J/#psi & track2','Dphi_JpsiPi1_':'#Delta#phi of J/#psi & track1','Dphi_JpsiPi2_':'#Delta#phi of J/#psi & track2','mumNHits_':'#mu^{-} total hits','mumNPHits_':'#mu^{-} pixel hits','mupNHits_':'#mu^{+} total hits','mupNPHits_':'#mu^{+} pixel hits','mu1soft_':'muon1 soft ID','mu2soft_':'muon2 soft ID','mu1tight_':'muon1 tight ID','mu2tight_':'muon2 tight ID','mu1loose_':'muon1 loose ID','mu2loose_':'muon2 loose ID','J_vertexchi2_':'#chi^{2} value of J/#psi vertex fit','Pi_nhits1_':'#pi1 total hits','Pi_npixelhits1_':'#pi1 pixel hits','Pi_nhits2_':'#pi2 total hits','Pi_npixelhits2_':'#pi2 pixel hits','Pi_vertexchisq1_':'#chi^{2} value of J/#psi+#pi vertex fit','Pi_vertexchisq2_':'#chi^{2} value of J/#psi+#pi#pi vertex fit'}
histos_xlabel= {'Pi_pt1_':'#pi1pT (GeV)','Pi_pt2_':'#pi2pT (GeV)','Pi_eta1_':'#eta of #pi1','Pi_eta2_':'#eta of #pi2','J_mass_':'M_{#mu^{+}#mu^{-}} (GeV)','pT_JPsi_':'pT_{#mu^{+}#mu^{-}} (GeV)','DeltaR_JpsiPi1_':'#DeltaR_{J/#psi#pi1}','DeltaR_JpsiPi2_':'#DeltaR_{J/#psi#pi2}','JPiPi_lxy_':'L_{xy,J/#psi#pi#pi}','Deta_JpsiPi1_':'#Delta#eta_{J/#psi#pi1}','Deta_JpsiPi2_':'#Delta#eta_{J/#psi#pi2}','Dphi_JpsiPi1_':'#Delta#phi_{J/#psi#pi1}','Dphi_JpsiPi2_':'#Delta#phi_{J/#psi#pi2}','mumNHits_':'N_{hits,#mu^{-}}','mumNPHits_':'N_{pixel hits,#mu^{-}}','mupNHits_':'N_{hits,#mu^{+}}','mupNPHits_':'N_{pixel hits,#mu^{+}}','mu1soft_':'muon1 soft ID','mu2soft_':'muon2 soft ID','mu1tight_':'muon1 tight ID','mu2tight_':'muon2 tight ID','mu1loose_':'muon1 loose ID','mu2loose_':'muon2 loose ID','J_vertexchi2_':'#chi^{2}_{J/#psi}','Pi_nhits1_':'N_{hits,#pi^{1}}','Pi_npixelhits1_':'N_{pixel hits,#pi^{1}}','Pi_nhits2_':'N_{hits,#pi^{2}}','Pi_npixelhits2_':'N_{pixel hits,#pi^{2}}','Pi_vertexchisq1_':'#chi^{2}_{J/#psi+#pi}','Pi_vertexchisq2_':'#chi^{2}_{J/#psi+#pi#pi}'}
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
        colors = [r.kGreen+3],
        legend_labels = ["Before Selection"],
        sigs = [h2],
        sig_labels = ["Signal Region"],
        options = {
            "legend_scalex": 0.7,
            "legend_scaley": 1.5,
            "cms_label": "",
            "lumi_value": "78.77",
            "do_stack": False,
            "yaxis_log": True,
            "output_name": "distribution_PLOTS/test0929_"+hist+".pdf",
            "yaxis_range": [0.1,1000000000],
            "yaxis_moreloglabels":False,
            "xaxis_label": histos_xlabel[hist],
            "title": histos[hist]
            #"yaxis_noexponents":True,
            #"output_ic": True,
            }
        )
