import ROOT as r
import os
import plottery as ply

histos= {
'PiPi_mass_':'#pi#pi invariant mass',
'Pion2_mass_':'pion2 invariant mass',
'Pion1_mass_':'pion1 invariant mass',
'DeltaR_PiPi_':'#DeltaR_{#pi#pi}',
'Dphi_PiPi_':'#Delta#phi_{#pi#pi}',
'Deta_PiPi_':'#Delta#eta_{#pi#pi}',
'Pi_pt1_':'pT of track1',
'Pi_pt2_':'pT of track2',
'Pi_eta1_':'#eta of track1',
'Pi_eta2_':'#eta of track2',
'J_mass_':'dimuon invariant mass',
'pT_JPsi_':'dimuon pT',
'DeltaR_JpsiPi1_':'#DeltaR of J/#psi & track1',
'DeltaR_JpsiPi2_':'#DeltaR of J/#psi & track2',
'JPiPi_lxy_':'Transverse L_{xy}',
'Deta_JpsiPi1_':'#Delta#eta of J/#psi & track1',
'Deta_JpsiPi2_':'#Delta#eta of J/#psi & track2',
'Dphi_JpsiPi1_':'#Delta#phi of J/#psi & track1',
'Dphi_JpsiPi2_':'#Delta#phi of J/#psi & track2',
'mumNHits_':'#mu^{-} total hits',
'mumNPHits_':'#mu^{-} pixel hits',
'mupNHits_':'#mu^{+} total hits',
'mupNPHits_':'#mu^{+} pixel hits',
'mu1soft_':'muon1 soft ID',
'mu2soft_':'muon2 soft ID',
'mu1tight_':'muon1 tight ID',
'mu2tight_':'muon2 tight ID',
'mu1loose_':'muon1 loose ID',
'mu2loose_':'muon2 loose ID',
'J_vertexchi2_':'J/#psi vertex fit#chi^{2}',
'Pi_nhits1_':'#pi1 total hits',
'Pi_npixelhits1_':'#pi1 pixel hits',
'Pi_nhits2_':'#pi2 total hits',
'Pi_npixelhits2_':'#pi2 pixel hits',
'Pi_vertexchisq1_':'J/#psi+#pi vertex fit #chi^{2}',
'Pi_vertexchisq2_':'J/#psi+#pi#pi vertex fit #chi^{2}',
'Pi1_hcalFraction_':'#pi1 hcalFraction',
'Pi2_hcalFraction_':'#pi2 hcalFraction',
'Pi1_vertexNdof_':'#pi1 vertexNdof',
'Pi2_vertexNdof_':'#pi2 vertexNdof',
'Pi1_vertexNchi2_':'#pi1 vertexNormalizedchi2',
'Pi2_vertexNchi2_':'#pi2 vertexNormalizedchi2',
'Pi1_lambda_':'#pi1 #lambda',
'Pi2_lambda_':'#pi2 #lambda',
'Pi1_lambdaError_':'#pi1 #lambda Error',
'Pi2_lambdaError_':'#pi2 #lambda Error',
'Pi1_qoverp_':'#pi1 q/p',
'Pi2_qoverp_':'#pi2 q/p',
'Pi1_qoverpError_':'#pi1 q/p Error',
'Pi2_qoverpError_':'#pi2 q/p Error',
'Pi1_validTkFraction_':'#pi1 validTkFraction',
'Pi2_validTkFraction_':'#pi2 validTkFraction',
'Pi1_numberOfMothers_':'#pi1 numberOfMothers',
'Pi2_numberOfMothers_':'#pi2 numberOfMothers',
'Pi1_numberOfSourceCandidatePtrs_':'#pi1 numberOfSourceCandidatePtrs',
'Pi2_numberOfSourceCandidatePtrs_':'#pi2 numberOfSourceCandidatePtrs',
'Pi1_pdgId_':'#pi1 pdgId',
'Pi2_pdgId_':'#pi2 pdgId',
'Pi1_numberOfValidHitsOnTrack_':'#pi1 numberOfValidHitsOnTrack',
'Pi2_numberOfValidHitsOnTrack_':'#pi2 numberOfValidHitsOnTrack',
'Pi1_innerDetId_':'#pi1 innerDetId',
'Pi2_innerDetId_':'#pi2 innerDetId',
'Pi1_innerOk_':'#pi1 innerOK',
'Pi2_innerOk_':'#pi2 innerOK',
'Pi1_isCaloMuon_':'#pi1 isCaloMuon',
'Pi2_isCaloMuon_':'#pi2 isCaloMuon',
'Pi1_isConvertedPhoton_':'#pi1 isConvertedPhoton',
'Pi2_isConvertedPhoton_':'#pi2 isConvertedPhoton',
'Pi1_isElectron_':'#pi1 isElectron',
'Pi2_isElectron_':'#pi2 isElectron',
#'Pi1_isMuon_':'#pi1 isMuon',
#'Pi2_isMuon_':'#pi2 isMuon',
#'Pi1_isPhoton_':'#pi1 isPhoton',
#'Pi2_isPhoton_':'#pi2 isPhoton',
'Pi1_isGlobalMuon_':'#pi1 isGlobalMuon',
'Pi2_isGlobalMuon_':'#pi2 isGlobalMuon',
'Pi1_isJet_':'#pi1 isJet',
'Pi2_isJet_':'#pi2 isJet',
'Pi1_isLonglived_':'#pi1 isLonglived',
'Pi2_isLonglived_':'#pi2 isLonglived',
'Pi1_massConstraint_':'#pi1 hasmassConstraint',
'Pi2_massConstraint_':'#pi2 hasmassConstraint',
'J_Prob_':'J/#psi vertexfit probability',
'JPi_Prob_':'J/#psi+#pi vertexfit probability',
'JPiPi_Prob_':'J/#psi+#pi#pi vertexfit probability'
}

histos_xlabel= {
'PiPi_mass_':'M_{#pi#pi} (GeV)',
'Pion2_mass_':'m_{#pi2}',
'Pion1_mass_':'m_{#pi1}',
'DeltaR_PiPi_':'#DeltaR_{#pi#pi}',
'Dphi_PiPi_':'#Delta#phi_{#pi#pi}',
'Deta_PiPi_':'#Delta#eta_{#pi#pi}',
'Pi_pt1_':'#pi1pT (GeV)',
'Pi_pt2_':'#pi2pT (GeV)',
'Pi_eta1_':'#eta of #pi1',
'Pi_eta2_':'#eta of #pi2',
'J_mass_':'M_{#mu^{+}#mu^{-}} (GeV)',
'pT_JPsi_':'pT_{#mu^{+}#mu^{-}} (GeV)',
'DeltaR_JpsiPi1_':'#DeltaR_{J/#psi#pi1}',
'DeltaR_JpsiPi2_':'#DeltaR_{J/#psi#pi2}',
'JPiPi_lxy_':'L_{xy,J/#psi#pi#pi}',
'Deta_JpsiPi1_':'#Delta#eta_{J/#psi#pi1}',
'Deta_JpsiPi2_':'#Delta#eta_{J/#psi#pi2}',
'Dphi_JpsiPi1_':'#Delta#phi_{J/#psi#pi1}',
'Dphi_JpsiPi2_':'#Delta#phi_{J/#psi#pi2}',
'mumNHits_':'N_{hits,#mu^{-}}',
'mumNPHits_':'N_{pixel hits,#mu^{-}}',
'mupNHits_':'N_{hits,#mu^{+}}',
'mupNPHits_':'N_{pixel hits,#mu^{+}}',
'mu1soft_':'muon1 soft ID',
'mu2soft_':'muon2 soft ID',
'mu1tight_':'muon1 tight ID',
'mu2tight_':'muon2 tight ID',
'mu1loose_':'muon1 loose ID',
'mu2loose_':'muon2 loose ID',
'J_vertexchi2_':'#chi^{2}_{J/#psi}',
'Pi_nhits1_':'N_{hits,#pi^{1}}',
'Pi_npixelhits1_':'N_{pixel hits,#pi^{1}}',
'Pi_nhits2_':'N_{hits,#pi^{2}}',
'Pi_npixelhits2_':'N_{pixel hits,#pi^{2}}',
'Pi_vertexchisq1_':'#chi^{2}_{J/#psi+#pi}',
'Pi_vertexchisq2_':'#chi^{2}_{J/#psi+#pi#pi}',
'Pi1_hcalFraction_':'',
'Pi2_hcalFraction_':'',
'Pi1_vertexNdof_':'Ndof',
'Pi2_vertexNdof_':'Ndof',
'Pi1_vertexNchi2_':'#chi^{2}',
'Pi2_vertexNchi2_':'#chi^{2}',
'Pi1_lambda_':'#lambda',
'Pi2_lambda_':'#lambda',
'Pi1_lambdaError_':'Error_{#lambda}',
'Pi2_lambdaError_':'Error_{#lambda}',
'Pi1_qoverp_':'',
'Pi2_qoverp_':'',
'Pi1_qoverpError_':'',
'Pi2_qoverpError_':'',
'Pi1_validTkFraction_':'',
'Pi2_validTkFraction_':'',
'Pi1_numberOfMothers_':'',
'Pi2_numberOfMothers_':'',
'Pi1_numberOfSourceCandidatePtrs_':'',
'Pi2_numberOfSourceCandidatePtrs_':'',
'Pi1_pdgId_':'',
'Pi2_pdgId_':'',
'Pi1_numberOfValidHitsOnTrack_':'',
'Pi2_numberOfValidHitsOnTrack_':'',
'Pi1_innerDetId_':'',
'Pi2_innerDetId_':'',
'Pi1_innerOk_':'',
'Pi2_innerOk_':'',
'Pi1_isCaloMuon_':'',
'Pi2_isCaloMuon_':'',
'Pi1_isConvertedPhoton_':'',
'Pi2_isConvertedPhoton_':'',
'Pi1_isElectron_':'',
'Pi2_isElectron_':'',
#'Pi1_isMuon_':'',
#'Pi2_isMuon_':'',
#'Pi1_isPhoton_':'',
#'Pi2_isPhoton_':'',
'Pi1_isGlobalMuon_':'',
'Pi2_isGlobalMuon_':'',
'Pi1_isJet_':'',
'Pi2_isJet_':'',
'Pi1_isLonglived_':'',
'Pi2_isLonglived_':'',
'Pi1_massConstraint_':'',
'Pi2_massConstraint_':'',
'J_Prob_':'',
'JPi_Prob_':'',
'JPiPi_Prob_':''
}
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
hm2 = r.TH2F()
file.GetObject("histos/M_JpsiPi1&M_JpsiPi2_total")
ply.plot_hist_2d(
        hm2,
        options = {
            #"zaxis_log": True,
            "bin_text_smart": True,
            "output_name": "distribution_PLOTS/dalitz.pdf",
            #"us_flag": True,
            "output_ic": True,
            #"zaxis_noexponents": True,
            }
        )
