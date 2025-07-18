NanoAOD versions
=================

CMS has different versions of the nanoAOD format, with each data-taking run typically introducing a new version that includes updates, bug fixes, or new features. 
For Run-2, the UL datasets were fully reprocessed using nanoAODv9. In Run-3, a new nanoAOD version is introduced each year, starting with nanoAODv12. 
The expected "final" version for both Run-2 and Run-3 will be nanoAODv15. 
However, these version updates also bring changes to branch data types in the nanoAODs. Since C++ code in CROWN requires explicit type definitions, 
these changes can cause compatibility issues. 
Below is a list of known differences between nanoAODv9 and nanoAODv12, and between nanoAODv12 and nanoAODv15.

.. list-table:: Changes from nanoAOD v9 to v12
   :widths: 75 100
   :header-rows: 1

   * - Branch 
     - Type or name change
   * - nCorrT1METJet
     - UInt_t -> Int_t
   * - Electron_cutBased
     - Int_t -> UChar_t
   * - Electron_genPartIdx
     - Int_t -> Short_t
   * - Electron_jetIdx
     - Int_t -> Short_t
   * - Electron_mvaFall17V2Iso
     - Electron_mvaIso
   * - Electron_mvaFall17V2noIso
     - Electron_mvaNoIso
   * - Electron_mvaFall17V2Iso_WP80
     - Electron_mvaIso_WP80
   * - Electron_mvaFall17V2Iso_WP90
     - Electron_mvaIso_WP90
   * - Electron_mvaFall17V2Iso_WPL
     - Electron_mvaIso_WPL
   * - Electron_mvaFall17V2noIso_WP80
     - Electron_mvaNoIso_WP80
   * - Electron_mvaFall17V2noIso_WP90
     - Electron_mvaNoIso_WP90
   * - Electron_mvaFall17V2noIso_WPL
     - Electron_mvaNoIso_WPL
   * - Electron_photonIdx
     - Int_t -> Short_t
   * - Electron_tightCharge
     - Int_t -> UChar_t
   * - nElectron
     - UInt_t -> Int_t
   * - FatJet_electronIdx3SJ
     - Int_t -> Short_t
   * - FatJet_genJetAK8Idx
     - Int_t -> Short_t
   * - FatJet_hadronFlavour
     - Int_t -> UChar_t
   * - FatJet_jetId
     - Int_t -> UChar_t
   * - FatJet_muonIdx3SJ
     - Int_t -> Short_t
   * - FatJet_particleNetMD_QCD
     - FatJet_particleNetLegacy_QCD
   * - FatJet_particleNetMD_Xbb
     - FatJet_particleNetLegacy_Xbb
   * - FatJet_particleNetMD_Xcc
     - FatJet_particleNetLegacy_Xcc
   * - FatJet_particleNetMD_Xqq
     - FatJet_particleNetLegacy_Xqq
   * - FatJet_particleNet_H4qvsQCD
     - FatJet_particleNetWithMass_H4qvsQCD
   * - FatJet_particleNet_HbbvsQCD
     - FatJet_particleNetWithMass_HbbvsQCD
   * - FatJet_particleNet_HccvsQCD
     - FatJet_particleNetWithMass_HccvsQCD
   * - FatJet_particleNet_TvsQCD
     - FatJet_particleNetWithMass_TvsQCD
   * - FatJet_particleNet_WvsQCD
     - FatJet_particleNetWithMass_WvsQCD
   * - FatJet_particleNet_ZvsQCD
     - FatJet_particleNetWithMass_ZvsQCD
   * - FatJet_particleNet_QCD
     - FatJet_particleNetWithMass_QCD
   * - FatJet_particleNet_mass
     - FatJet_particleNetLegacy_mass
   * - FatJet_subJetIdx1
     - Int_t -> Short_t
   * - FatJet_subJetIdx2
     - Int_t -> Short_t
   * - nFatJet
     - UInt_t -> Int_t
   * - FsrPhoton_muonIdx
     - Int_t -> Short_t
   * - nFsrPhoton
     - UInt_t -> Int_t
   * - nGenDressedLepton
     - UInt_t -> Int_t
   * - nGenIsolatedPhoton
     - UInt_t -> Int_t
   * - GenJet_partonFlavour
     - Int_t -> Short_t
   * - nGenJet
     - UInt_t -> Int_t
   * - GenJetAK8_partonFlavour
     - Int_t -> Short_t
   * - nGenJetAK8
     - UInt_t -> Int_t
   * - GenPart_genPartIdxMother
     - Int_t -> Short_t
   * - GenPart_statusFlags
     - Int_t -> UShort_t
   * - nGenPart
     - UInt_t -> Int_t
   * - GenVisTau_charge
     - Int_t -> Short_t
   * - GenVisTau_genPartIdxMother
     - Int_t -> Short_t
   * - GenVisTau_status
     - Int_t -> UChar_t
   * - nGenVisTau
     - UInt_t -> Int_t
   * - IsoTrack_charge
     - Int_t -> Short_t
   * - IsoTrack_fromPV
     - Int_t -> Short_t
   * - nIsoTrack
     - UInt_t -> Int_t
   * - Jet_electronIdx1
     - Int_t -> Short_t
   * - Jet_electronIdx2
     - Int_t -> Short_t
   * - Jet_genJetIdx
     - Int_t -> Short_t
   * - Jet_hadronFlavour
     - Int_t -> UChar_t
   * - Jet_jetId
     - Int_t -> UChar_t
   * - Jet_muonIdx1
     - Int_t -> Short_t
   * - Jet_muonIdx2
     - Int_t -> Short_t
   * - Jet_nElectrons
     - Int_t -> UChar_t
   * - Jet_nMuons
     - Int_t -> UChar_t
   * - Jet_puId
     - Int_t -> UChar_t
   * - Jet_partonFlavour
     - Int_t -> Short_t
   * - nJet
     - UInt_t -> Int_t
   * - LowPtElectron_convWP
     - Int_t -> UChar_t
   * - LowPtElectron_genPartIdx
     - Int_t -> Short_t
   * - nLowPtElectron
     - UInt_t -> Int_t
   * - Muon_fsrPhotonIdx
     - Int_t -> Short_t
   * - Muon_genPartIdx
     - Int_t -> Short_t
   * - Muon_jetIdx
     - Int_t -> Short_t
   * - Muon_nStations
     - Int_t -> UChar_t
   * - Muon_nTrackerLayers
     - Int_t -> UChar_t
   * - Muon_tightCharge
     - Int_t -> UChar_t
   * - nMuon
     - UInt_t -> Int_t
   * - nOtherPV
     - UInt_t -> Int_t
   * - nPSWeight
     - UInt_t -> Int_t
   * - PV_npvs
     - Int_t -> UChar_t
   * - PV_npvsGood
     - Int_t -> UChar_t
   * - Photon_cutBased
     - Int_t -> UChar_t
   * - Photon_electronIdx
     - Int_t -> Short_t
   * - Photon_jetIdx
     - Int_t -> Short_t
   * - Photon_genPartIdx
     - Int_t -> Short_t
   * - nPhoton
     - UInt_t -> Int_t
   * - SV_charge
     - Int_t -> Short_t
   * - nSV
     - UInt_t -> Int_t
   * - nSoftActivityJet
     - UInt_t -> Int_t
   * - nSubGenJetAK8
     - UInt_t -> Int_t
   * - SubJet_hadronFlavour
     - Int_t -> UChar_t
   * - nSubJet
     - UInt_t -> Int_t
   * - Tau_charge
     - Int_t -> Short_t
   * - Tau_decayMode
     - Int_t -> UChar_t
   * - Tau_idDeepTau2017v2p1VSe
     - same type, content changed
   * - Tau_idDeepTau2017v2p1VSjet
     - same type, content changed
   * - Tau_idDeepTau2017v2p1VSmu
     - same type, content changed
   * - Tau_jetIdx
     - Int_t -> Short_t
   * - Tau_genPartIdx
     - Int_t -> Short_t
   * - nTau
     - UInt_t -> Int_t
   * - TrigObj_id
     - Int_t -> UShort_t
   * - TrigObj_l1charge
     - Int_t -> Short_t
   * - nTrigObj
     - UInt_t -> Int_t
   * - boostedTau_idAntiEle2018
     - same type, content changed
   * - boostedTau_idAntiMu
     - same type, content changed
   * - boostedTau_idMVAnewDM2017v2
     - same type, content changed
   * - boostedTau_idMVAoldDM2017v2
     - same type, content changed
   * - boostedTau_idMVAoldDMdR032017v2
     - same type, content changed
   * - boostedTau_jetIdx
     - Int_t -> Short_t
   * - boostedTau_rawAntiEleCat2018
     - Int_t -> Short_t
   * - boostedTau_genPartIdx
     - Int_t -> Short_t
   * - nboostedTau
     - UInt_t -> Int_t
   * - fixedGridRhoFastjetAll
     - Rho_fixedGridRhoFastjetAll
   * - fixedGridRhoFastjetCentral
     - Rho_fixedGridRhoFastjetCentral
   * - fixedGridRhoFastjetCentralCalo
     - Rho_fixedGridRhoFastjetCentralCalo
   * - fixedGridRhoFastjetCentralChargedPileUp
     - Rho_fixedGridRhoFastjetCentralChargedPileUp
   * - fixedGridRhoFastjetCentralNeutral
     - Rho_fixedGridRhoFastjetCentralNeutral

.. list-table:: Changes from nanoAOD v12 to v15
   :widths: 75 100
   :header-rows: 1

   * - Branch
     - Type or name change
   * - BeamSpot_type
     - Char_t -> Short_t
   * - Electron_seediEtaOriX
     - Char_t -> Short_t
   * - Electron_seediPhiOriY
     - Int_t -> Short_t
   * - Photon_seediEtaOriX
     - Char_t -> Short_t
   * - Photon_seediPhiOriY
     - Int_t -> Short_t
   * - TrigObj_filterBits
     - Int_t -> ULong64_t
   * - HLT_AK8DiPFJet250_250_MassSD50
     - HLT_AK8DiPFJet250_250_SoftDropMass50
   * - HLT_AK8DiPFJet260_260_MassSD30
     - HLT_AK8DiPFJet260_260_SoftDropMass30
   * - HLT_AK8DiPFJet270_270_MassSD30
     - HLT_AK8DiPFJet270_270_SoftDropMass30
   * - HLT_AK8PFJet400_MassSD30
     - HLT_AK8PFJet400_SoftDropMass30
   * - HLT_AK8PFJet450_MassSD30
     - HLT_AK8PFJet450_SoftDropMass30
   * - MET_covXX
     - PFMET_covXX
   * - MET_covXY
     - PFMET_covXY
   * - MET_covYY
     - PFMET_covYY
   * - MET_phi
     - PFMET_phi
   * - MET_pt
     - PFMET_pt
   * - MET_significance
     - PFMET_significance
   * - MET_sumEt
     - PFMET_sumEt
   * - MET_sumPtUnclustered
     - PFMET_sumPtUnclustered
   * - RawMET_phi
     - RawPFMET_phi
   * - RawMET_pt
     - RawPFMET_pt
   * - RawMET_sumEt
     - RawPFMET_sumEt
   * - TkMET_phi
     - TrkMET_phi
   * - TkMET_pt
     - TrkMET_pt
   * - TkMET_sumEt
     - TrkMET_sumEt

.. list-table:: Changes from nanoAOD v12 to v15
   :widths: 75 100
   :header-rows: 1

   * - Dropped out
  - New variables
* - ChsMET_phi
  - boostedTau_rawBoostedDeepTauRunIIv2p0VSe
* - ChsMET_pt
  - boostedTau_rawBoostedDeepTauRunIIv2p0VSjet
* - ChsMET_sumEt
  - boostedTau_rawBoostedDeepTauRunIIv2p0VSmu
* - Electron_mvaTTH
  - CorrT1METJet_EmEF
* - FatJet_btagDDBvLV2
  - CorrT1METJet_muonSubtrDeltaEta
* - FatJet_btagDDCvBV2
  - CorrT1METJet_muonSubtrDeltaPhi
* - FatJet_btagDDCvLV2
  - CorrT1METJet_rawMass
* - FatJet_btagDeepB
  - Dataset_ScoutingPFMonitor
* - FatJet_btagHbb
  - Dataset_ScoutingPFRun3
* - FatJet_jetId
  - DST_PFScouting_AXOLoose
* - FatJet_nBHadrons
  - DST_PFScouting_AXONominal
* - FatJet_nCHadrons
  - DST_PFScouting_AXOTight
* - Flag_METFilters
  - DST_PFScouting_AXOVLoose
* - HLT_AK4CaloJet100
  - DST_PFScouting_AXOVTight
* - HLT_AK4CaloJet120
  - DST_PFScouting_CICADALoose
* - HLT_AK4CaloJet30
  - DST_PFScouting_CICADAMedium
* - HLT_AK4CaloJet40
  - DST_PFScouting_CICADATight
* - HLT_AK4CaloJet50
  - DST_PFScouting_CICADAVLoose
* - HLT_AK4CaloJet80
  - DST_PFScouting_CICADAVTight
* - HLT_AK4PFJet100
  - DST_PFScouting_DatasetMuon
* - HLT_AK4PFJet120
  - DST_PFScouting_DoubleEG
* - HLT_AK4PFJet30
  - DST_PFScouting_DoubleMuon
* - HLT_AK4PFJet50
  - DST_PFScouting_JetHT
* - HLT_AK4PFJet80
  - DST_PFScouting_SingleMuon
* - HLT_AK8DiPFJet250_250_MassSD30
  - DST_PFScouting_SinglePhotonEB
* - HLT_AK8PFHT750_TrimMass50
  - DST_PFScouting_ZeroBias
* - HLT_AK8PFHT800_TrimMass50
  - Electron_ecalEnergy
* - HLT_AK8PFHT850_TrimMass50
  - Electron_ecalEnergyError
* - HLT_AK8PFHT900_TrimMass50
  - Electron_fbrem
* - HLT_AK8PFJet15
  - Electron_gsfTrketaMode
* - HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35
  - Electron_gsfTrkphiMode
* - HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30
  - Electron_gsfTrkpMode
* - HLT_AK8PFJet25
  - Electron_gsfTrkpModeErr
* - HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35
  - Electron_ipLengthSig
* - HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetTauTau0p30
  - Electron_IPx
* - HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetBB0p35
  - Electron_IPy
* - HLT_AK8PFJet275_SoftDropMass40_PFAK8ParticleNetTauTau0p30
  - Electron_IPz
* - HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2
  - Electron_isEB
* - HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4
  - Electron_isEcalDriven
* - HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02
  - Electron_jetDF
* - HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1
  - Electron_mvaIso_WPHZZ
* - HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17
  - Electron_pfRelIso04_all
* - HLT_AK8PFJet360_TrimMass30
  - Electron_PreshowerEnergy
* - HLT_AK8PFJet380_TrimMass30
  - Electron_promptMVA
* - HLT_AK8PFJet400_SoftDropMass40
  - Electron_rawEnergy
* - HLT_AK8PFJet400_TrimMass30
  - Electron_superclusterEta
* - HLT_AK8PFJet420_MassSD30
  - FatJet_chEmEF
* - HLT_AK8PFJet420_TrimMass30
  - FatJet_chHEF
* - HLT_AK8PFJet425_SoftDropMass40
  - FatJet_chMultiplicity
* - HLT_AK8PFJet450_SoftDropMass40
  - FatJet_globalParT3_massCorrGeneric
* - HLT_AK8PFJetFwd15
  - FatJet_globalParT3_massCorrX2p
* - HLT_AK8PFJetFwd25
  - FatJet_globalParT3_QCD
* - HLT_CaloMET100_NotCleaned
  - FatJet_globalParT3_TopbWev
* - HLT_CaloMET110_NotCleaned
  - FatJet_globalParT3_TopbWmv
* - HLT_CaloMET250_NotCleaned
  - FatJet_globalParT3_TopbWq
* - HLT_CaloMET300_NotCleaned
  - FatJet_globalParT3_TopbWqq
* - HLT_CaloMET80_NotCleaned
  - FatJet_globalParT3_TopbWtauhv
* - HLT_DiJet110_35_Mjj650_PFMET110
  - FatJet_globalParT3_withMassTopvsQCD
* - HLT_DiJet110_35_Mjj650_PFMET120
  - FatJet_globalParT3_withMassWvsQCD
* - HLT_DiJet110_35_Mjj650_PFMET130
  - FatJet_globalParT3_withMassZvsQCD
* - HLT_Dimuon0_LowMass_L1_0er1p5R
  - FatJet_globalParT3_WvsQCD
* - HLT_Dimuon0_LowMass_L1_4R
  - FatJet_globalParT3_Xbb
* - HLT_Dimuon0_Upsilon_L1_4p5NoOS
  - FatJet_globalParT3_Xcc
* - HLT_Dimuon0_Upsilon_L1_5
  - FatJet_globalParT3_Xcs
* - HLT_Dimuon0_Upsilon_L1_5M
  - FatJet_globalParT3_Xqq
* - HLT_Dimuon0_Upsilon_Muon_L1_TM0
  - FatJet_globalParT3_Xtauhtaue
* - HLT_Dimuon10_PsiPrime_Barrel_Seagulls
  - FatJet_globalParT3_Xtauhtauh
* - HLT_Dimuon20_Jpsi_Barrel_Seagulls
  - FatJet_globalParT3_Xtauhtaum
* - HLT_DiPFJet15_FBEta3_NoCaloMatched
  - FatJet_globalParT3_XWW3q
* - HLT_DiPFJet15_NoCaloMatched
  - FatJet_globalParT3_XWW4q
* - HLT_DiPFJet25_FBEta3_NoCaloMatched
  - FatJet_globalParT3_XWWqqev
* - HLT_DiPFJet25_NoCaloMatched
  - FatJet_globalParT3_XWWqqmv
* - HLT_DiPFJetAve15_HFJEC
  - FatJet_hfEmEF
* - HLT_DiPFJetAve25_HFJEC
  - FatJet_hfHEF
* - HLT_DiPFJetAve35_HFJEC
  - FatJet_muEF
* - HLT_DiPhoton10sminlt0p1
  - FatJet_neEmEF
* - HLT_DiPhoton10sminlt0p12
  - FatJet_neHEF
* - HLT_DiPhoton10sminlt0p14
  - FatJet_neMultiplicity
* - HLT_DiPhoton10sminlt0p16
  - FatJet_particleNet_WVsQCD
* - HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55
  - FatJet_particleNetLegacy_mass
* - HLT_DoubleEle4_eta1p22_mMax6
  - FatJet_particleNetLegacy_QCD
* - HLT_DoubleEle4p5_eta1p22_mMax6
  - FatJet_particleNetLegacy_Xbb
* - HLT_DoubleEle5_eta1p22_mMax6
  - FatJet_particleNetLegacy_Xcc
* - HLT_DoubleEle5p5_eta1p22_mMax6
  - FatJet_particleNetLegacy_Xqq
* - HLT_DoubleEle6_eta1p22_mMax6
  - FatJetPFCand_jetIdx
* - HLT_DoubleEle7_eta1p22_mMax6
  - FatJetPFCand_pfCandIdx
* - HLT_DoubleEle7p5_eta1p22_mMax6
  - FiducialMET_phi
* - HLT_DoubleEle8p5_eta1p22_mMax6
  - FiducialMET_pt
* - HLT_DoubleEle9_eta1p22_mMax6
  - GenJet_nBHadrons
* - HLT_DoubleEle9p5_eta1p22_mMax6
  - GenJet_nCHadrons
* - HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1
  - GenJetAK8_nBHadrons
* - HLT_DoubleMu20_7_Mass0to30_L1_DM4
  - GenJetAK8_nCHadrons
* - HLT_DoubleMu20_7_Mass0to30_L1_DM4EG
  - GenPart_iso
* - HLT_DoubleMu20_7_Mass0to30_Photon23
  - HLT_AK8DiPFJet250_250_SoftDropMass40
* - HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi
  - HLT_AK8DiPFJet260_260_SoftDropMass40
* - HLT_DoubleMu40NoFiltersNoVtxDisplaced
  - HLT_AK8DiPFJet280_280_SoftDropMass30
* - HLT_DoublePFJets100_PFBTagDeepCSV_p71
  - HLT_AK8DiPFJet290_290_SoftDropMass30
* - HLT_DoublePFJets100_PFBTagDeepJet_p71
  - HLT_AK8PFJet220_SoftDropMass40
* - HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepCSV_p71
  - HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50
* - HLT_DoublePFJets116MaxDeta1p6_DoublePFBTagDeepJet_p71
  - HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53
* - HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepCSV_p71
  - HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55
* - HLT_DoublePFJets128MaxDeta1p6_DoublePFBTagDeepJet_p71
  - HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60
* - HLT_DoublePFJets200_PFBTagDeepCSV_p71
  - HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06
* - HLT_DoublePFJets200_PFBTagDeepJet_p71
  - HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10
* - HLT_DoublePFJets350_PFBTagDeepCSV_p71
  - HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03
* - HLT_DoublePFJets350_PFBTagDeepJet_p71
  - HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05
* - HLT_DoublePFJets40_Mass500_MediumDeepTauPFTauHPS45_L2NN_MediumDeepTauPFTauHPS20_eta2p1
  - HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06
* - HLT_DoublePFJets40_PFBTagDeepCSV_p71
  - HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10
* - HLT_DoublePFJets40_PFBTagDeepJet_p71
  - HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03
* - HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1
  - HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05
* - HLT_DoubleTrkMu_16_6_NoFiltersNoVtx
  - HLT_AK8PFJet275_Nch40
* - HLT_Ele145_CaloIdVT_GsfTrkIdT
  - HLT_AK8PFJet275_Nch45
* - HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30
  - HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06
* - HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL
  - HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10
* - HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5
  - HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03
* - HLT_Ele15_WPLoose_Gsf
  - HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05
* - HLT_Ele17_WPLoose_Gsf
  - HLT_AK8PFJet380_SoftDropMass30
* - HLT_Ele200_CaloIdVT_GsfTrkIdT
  - HLT_AK8PFJet425_SoftDropMass30
* - HLT_Ele20_eta2p1_WPLoose_Gsf
  - HLT_CscCluster100_Ele5
* - HLT_Ele20_WPLoose_Gsf
  - HLT_CscCluster100_Mu5
* - HLT_Ele20_WPTight_Gsf
  - HLT_CscCluster100_PNetTauhPFJet10_Loose
* - HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1
  - HLT_CscCluster50_Photon20Unseeded
* - HLT_Ele250_CaloIdVT_GsfTrkIdT
  - HLT_CscCluster50_Photon30Unseeded
* - HLT_Ele27_Ele37_CaloIdL_MW
  - HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP
* - HLT_Ele27_WPTight_Gsf
  - HLT_DiPFJetAve260_HFJEC
* - HLT_Ele28_WPTight_Gsf
  - HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT
* - HLT_Ele300_CaloIdVT_GsfTrkIdT
  - HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT
* - HLT_Ele35_WPTight_Gsf_L1EGMT
  - HLT_DiphotonMVA14p25_Mass90
* - HLT_ExpressMuons
  - HLT_DiphotonMVA14p25_Tight_Mass90
* - HLT_HcalIsolatedbunch
  - HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24
* - HLT_HT200_L1SingleLLPJet_DisplacedDijet30_Inclusive1PtrkShortSig5
  - HLT_DoubleCscCluster100
* - HLT_HT200_L1SingleLLPJet_DisplacedDijet35_Inclusive1PtrkShortSig5
  - HLT_DoubleCscCluster75
* - HLT_HT430_DisplacedDijet30_Inclusive1PtrkShortSig5
  - HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm
* - HLT_HT430_DisplacedDijet35_Inclusive1PtrkShortSig5
  - HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm
* - HLT_HT430_DisplacedDijet60_DisplacedTrack
  - HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm
* - HLT_HT450_Beamspot
  - HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm
* - HLT_HT500_DisplacedDijet40_DisplacedTrack
  - HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy
* - HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1
  - HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1
* - HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1
  - HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng
* - HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1
  - HLT_DoubleMu2_Jpsi_LowPt
* - HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1
  - HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0
* - HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1
  - HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA
* - HLT_IsoMu27_MET90
  - HLT_DoubleMu4_3_LowMass_SS
* - HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1
  - HLT_DoublePFJets100_PNetBTag_0p11
* - HLT_IsoMu30
  - HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11
* - HLT_L1NotBptxOR
  - HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11
* - HLT_L1SingleMu18
  - HLT_DoublePFJets200_PNetBTag_0p11
* - HLT_L1SingleMu25
  - HLT_DoublePFJets350_PNetBTag_0p11
* - HLT_L1UnpairedBunchBptxMinus
  - HLT_DoublePFJets40_PNetBTag_0p11
* - HLT_L1UnpairedBunchBptxPlus
  - HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60
* - HLT_L2Mu10
  - HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75
* - HLT_L2Mu50
  - HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3
* - HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1
  - HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3
* - HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1
  - HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53
* - HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1
  - HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1
* - HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight
  - HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1
* - HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight
  - HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1
* - HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight
  - HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40
* - HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight
  - HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06
* - HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60
  - HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40
* - HLT_Mu12
  - HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06
* - HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71
  - HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10
* - HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71
  - HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7
* - HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71
  - HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8
* - HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71
  - HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5
* - HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71
  - HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5
* - HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71
  - HLT_HT350
* - HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71
  - HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive
* - HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71
  - HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive
* - HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71
  - HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive
* - HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71
  - HLT_HT350_DelayedJet40_SingleDelay3nsInclusive
* - HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71
  - HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive
* - HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71
  - HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive
* - HLT_Mu12_DoublePhoton20
  - HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5
* - HLT_Mu12_IP6_part0
  - HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5
* - HLT_Mu12_IP6_part1
  - HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5
* - HLT_Mu12_IP6_part2
  - HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5
* - HLT_Mu12_IP6_part3
  - HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5
* - HLT_Mu12_IP6_part4
  - HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive
* - HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5
  - HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless
* - HLT_Mu18_Mu9
  - HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless
* - HLT_Mu18_Mu9_DZ
  - HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive
* - HLT_Mu18_Mu9_SameSign_DZ
  - HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive
* - HLT_Mu20_Mu10
  - HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive
* - HLT_Mu20_Mu10_DZ
  - HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless
* - HLT_Mu20_Mu10_SameSign
  - HLT_HT430_DelayedJet40_SingleDelay1nsInclusive
* - HLT_Mu20_Mu10_SameSign_DZ
  - HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive
* - HLT_Mu20_TkMu0_Phi
  - HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless
* - HLT_Mu23_Mu12
  - HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive
* - HLT_Mu23_Mu12_DZ
  - HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive
* - HLT_Mu23_Mu12_SameSign
  - HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless
* - HLT_Mu23_Mu12_SameSign_DZ
  - HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive
* - HLT_Mu25_TkMu0_Onia
  - HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive
* - HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight
  - HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive
* - HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight
  - HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1
* - HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60
  - HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1
* - HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5
  - HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1
* - HLT_Mu7_IP4_part0
  - HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1
* - HLT_Mu7_IP4_part1
  - HLT_IsoMu24_eta2p1_PFHT250
* - HLT_Mu7_IP4_part2
  - HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25
* - HLT_Mu7_IP4_part3
  - HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50
* - HLT_Mu7_IP4_part4
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP3_part0
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP3_part1
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP3_part2
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1
* - HLT_Mu8_IP3_part3
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP3_part4
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60
* - HLT_Mu8_IP5_part0
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75
* - HLT_Mu8_IP5_part1
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring
* - HLT_Mu8_IP5_part2
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring
* - HLT_Mu8_IP5_part3
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP5_part4
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring
* - HLT_Mu8_IP6_part0
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP6_part1
  - HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1
* - HLT_Mu8_IP6_part2
  - HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50
* - HLT_Mu8_IP6_part3
  - HLT_IsoMu24_OneProng32
* - HLT_Mu8_IP6_part4
  - HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1
* - HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60
  - HLT_IsoMu50_AK8PFJet220_SoftDropMass40
* - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5
  - HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06
* - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5
  - HLT_IsoMu50_AK8PFJet230_SoftDropMass40
* - HLT_Mu9_IP0_part0
  - HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06
* - HLT_Mu9_IP3_part0
  - HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10
* - HLT_Mu9_IP4_part0
  - HLT_IsoTrk200_L1SingleMuShower
* - HLT_Mu9_IP4_part1
  - HLT_IsoTrk400_L1SingleMuShower
* - HLT_Mu9_IP4_part2
  - HLT_L1AXOVTight
* - HLT_Mu9_IP4_part3
  - HLT_L1SingleLLPJet
* - HLT_Mu9_IP4_part4
  - HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless
* - HLT_Mu9_IP5_part0
  - HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive
* - HLT_Mu9_IP5_part1
  - HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless
* - HLT_Mu9_IP5_part2
  - HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive
* - HLT_Mu9_IP5_part3
  - HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless
* - HLT_Mu9_IP5_part4
  - HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive
* - HLT_Mu9_IP6_part0
  - HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless
* - HLT_Mu9_IP6_part1
  - HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive
* - HLT_Mu9_IP6_part2
  - HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless
* - HLT_Mu9_IP6_part3
  - HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive
* - HLT_Mu9_IP6_part4
  - HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive
* - HLT_OnlineMonitorGroup
  - HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless
* - HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5
  - HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive
* - HLT_PFHT350MinPFJet15
  - HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless
* - HLT_PFHT400_FivePFJet_100_100_60_30_30
  - HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive
* - HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5
  - HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive
* - HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepJet_4p5
  - HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive
* - HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5
  - HLT_L2Mu10NoVtx_2Cha_CosmicSeed
* - HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepJet_4p5
  - HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm
* - HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94
  - HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm
* - HLT_PFHT400_SixPFJet32_DoublePFBTagDeepJet_2p94
  - HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm
* - HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59
  - HLT_L3Mu30NoVtx_DxyMin0p01cm
* - HLT_PFHT450_SixPFJet36_PFBTagDeepJet_1p59
  - HLT_L3Mu50NoVtx_DxyMin0p01cm
* - HLT_PFHT700_PFMET95_PFMHT95_IDTight
  - HLT_Mu0_Barrel
* - HLT_PFHT800_PFMET85_PFMHT85_IDTight
  - HLT_Mu0_Barrel_L1HP10
* - HLT_PFJet15
  - HLT_Mu0_Barrel_L1HP11
* - HLT_PFJet25
  - HLT_Mu0_Barrel_L1HP6
* - HLT_PFJetFwd15
  - HLT_Mu0_Barrel_L1HP6_IP6
* - HLT_PFJetFwd25
  - HLT_Mu0_Barrel_L1HP7
* - HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1
  - HLT_Mu0_Barrel_L1HP8
* - HLT_PFMET100_PFMHT100_IDTight_PFHT60
  - HLT_Mu0_Barrel_L1HP9
* - HLT_PFMET105_PFJet100_looseRecoiling
  - HLT_Mu10_Barrel_L1HP11_IP6
* - HLT_PFMET110_PFJet100
  - HLT_Mu12_DoublePFJets100_PNetBTag_0p11
* - HLT_PFMET110_PFJet100_looseRecoiling
  - HLT_Mu12_DoublePFJets200_PNetBTag_0p11
* - HLT_PFMET110_PFMHT110_IDTight
  - HLT_Mu12_DoublePFJets350_PNetBTag_0p11
* - HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1
  - HLT_Mu12_DoublePFJets40_PNetBTag_0p11
* - HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1
  - HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11
* - HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1
  - HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11
* - HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1
  - HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53
* - HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60
  - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30
* - HLT_PFMETNoMu110_PFMHTNoMu110_IDTight
  - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30
* - HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60
  - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30
* - HLT_PFMETTypeOne110_PFMHT110_IDTight
  - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30
* - HLT_PFMETTypeOne120_PFMHT120_IDTight
  - HLT_Mu50_L1SingleMuShower
* - HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60
  - HLT_Mu6_Barrel_L1HP7_IP6
* - HLT_PFMETTypeOne130_PFMHT130_IDTight
  - HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5
* - HLT_Photon100EB_TightID_TightIso
  - HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5
* - HLT_Photon100EE_TightID_TightIso
  - HLT_Mu7_Barrel_L1HP8_IP6
* - HLT_Photon100EEHE10
  - HLT_Mu8_Barrel_L1HP9_IP6
* - HLT_Photon120EB_TightID_TightIso
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30
* - HLT_Photon20
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50
* - HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50
* - HLT_Photon60_R9Id90_CaloIdL_IsoL
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250
* - HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25
* - HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20
* - HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55
* - HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280
* - HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30
* - HLT_Photon90_CaloIdL_PFHT700
  - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55
* - HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
  - HLT_Mu9_Barrel_L1HP10_IP6
* - HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1
  - HLT_PFHT250_QuadPFJet25
* - HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2
  - HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50
* - HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2
  - HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55
* - HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
  - HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50
* - HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1
  - HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55
* - HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2
  - HLT_PFHT280_QuadPFJet30
* - HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2
  - HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50
* - HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
  - HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55
* - HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1
  - HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60
* - HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2
  - HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60
* - HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2
  - HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0
* - HLT_QuadPFJet70_50_40_30
  - HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3
* - HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65
  - HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70
* - HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65
  - HLT_PFHT400_FivePFJet_120_120_60_30_30
* - HLT_QuadPFJet70_50_45_35_PFBTagParticleNet_2BTagSum0p65
  - HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3
* - HLT_QuadPFJet98_83_71_15
  - HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6
* - HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
  - HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50
* - HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepJet_1p3_7p7_VBF1
  - HLT_PFHT450_SixPFJet36_PNetBTag0p35
* - HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2
  - HLT_PFJet110
* - HLT_QuadPFJet98_83_71_15_PFBTagDeepJet_1p3_VBF2
  - HLT_PFJet200_TimeGt2p5ns
* - HLT_Rsq0p35
  - HLT_PFJet200_TimeLtNeg2p5ns
* - HLT_Rsq0p40
  - HLT_PFJet40_GPUvsCPU
* - HLT_RsqMR300_Rsq0p09_MR200
  - HLT_Photon110EB_TightID_TightIso_AK8CaloJet30
* - HLT_RsqMR300_Rsq0p09_MR200_4jet
  - HLT_Photon110EB_TightID_TightIso_AK8PFJet30
* - HLT_RsqMR320_Rsq0p09_MR200
  - HLT_Photon110EB_TightID_TightIso_CaloJet30
* - HLT_RsqMR320_Rsq0p09_MR200_4jet
  - HLT_Photon110EB_TightID_TightIso_PFJet30
* - HLT_SingleJet30_Mu12_SinglePFJet40
  - HLT_Photon32_OneProng32_M50To105
* - HLT_SinglePhoton10_Eta3p1ForPPRef
  - HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34
* - HLT_SinglePhoton20_Eta3p1ForPPRef
  - HLT_Photon40EB
* - HLT_SinglePhoton30_Eta3p1ForPPRef
  - HLT_Photon40EB_TightID_TightIso
* - HLT_TripleJet110_35_35_Mjj650_PFMET110
  - HLT_Photon45EB
* - HLT_TripleJet110_35_35_Mjj650_PFMET120
  - HLT_Photon45EB_TightID_TightIso
* - HLT_TripleJet110_35_35_Mjj650_PFMET130
  - HLT_Photon50_TimeGt2p5ns
* - HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx
  - HLT_Photon50_TimeLtNeg2p5ns
* - HLT_TrkMu16NoFiltersNoVtx
  - HLT_Photon50EB
* - HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx
  - HLT_Photon50EB_TightID_TightIso
* - HLT_TrkMu6NoFiltersNoVtx
  - HLT_Photon50EB_TightID_TightIso_AK8CaloJet30
* - HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1
  - HLT_Photon50EB_TightID_TightIso_AK8PFJet30
* - HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1
  - HLT_Photon50EB_TightID_TightIso_CaloJet30
* - HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1
  - HLT_Photon50EB_TightID_TightIso_PFJet30
* - Jet_btagRobustParTAK4B
  - HLT_Photon55EB_TightID_TightIso
* - Jet_btagRobustParTAK4CvB
  - HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350
* - Jet_btagRobustParTAK4CvL
  - HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380
* - Jet_btagRobustParTAK4QG
  - HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400
* - Jet_jetId
  - HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3
* - L1_DoubleEG10_er1p2_dR_Max0p6
  - HLT_Photon75EB_TightID_TightIso
* - L1_DoubleEG10p5_er1p2_dR_Max0p6
  - HLT_Photon90EB_TightID_TightIso
* - L1_DoubleEG4_er1p2_dR_Max0p9
  - HLT_PPSRandom
* - L1_DoubleEG4p5_er1p2_dR_Max0p9
  - HLT_QuadPFJet100_88_70_30
* - L1_DoubleEG5_er1p2_dR_Max0p9
  - HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight
* - L1_DoubleEG5p5_er1p2_dR_Max0p8
  - HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1
* - L1_DoubleEG6_er1p2_dR_Max0p8
  - HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2
* - L1_DoubleEG6p5_er1p2_dR_Max0p8
  - HLT_QuadPFJet105_88_75_30
* - L1_DoubleEG7_er1p2_dR_Max0p8
  - HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight
* - L1_DoubleEG7p5_er1p2_dR_Max0p7
  - HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1
* - L1_DoubleEG8_er1p2_dR_Max0p7
  - HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2
* - L1_DoubleEG8er2p5_HTT260er
  - HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1
* - L1_DoubleEG8er2p5_HTT340er
  - HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2
* - L1_DoubleEG8p5_er1p2_dR_Max0p7
  - HLT_QuadPFJet111_90_80_30
* - L1_DoubleEG9_er1p2_dR_Max0p7
  - HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight
* - L1_DoubleEG9p5_er1p2_dR_Max0p6
  - HLT_SingleEle8
* - L1_DoubleEG_LooseIso20_10_er2p5
  - HLT_SingleEle8_SingleEGL1
* - L1_DoubleEG_LooseIso22_10_er2p5
  - HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3
* - L1_DoubleIsoTau28er2p1_Mass_Max80
  - HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3
* - L1_DoubleIsoTau28er2p1_Mass_Max90
  - HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3
* - L1_DoubleIsoTau30er2p1_Mass_Max80
  - HLT_VBF_DiPFJet125_45_Mjj1050
* - L1_DoubleIsoTau30er2p1_Mass_Max90
  - HLT_VBF_DiPFJet125_45_Mjj1200
* - L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5
  - HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1
* - L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5
  - HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3
* - L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5
  - HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1
* - L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp
  - HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3
* - L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5
  - HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf
* - L1_DoubleJet_100_30_DoubleJet30_Mass_Min620
  - HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf
* - L1_DoubleJet_110_35_DoubleJet35_Mass_Min620
  - HLT_VBF_DiPFJet50_Mjj650_Photon22
* - L1_DoubleJet_115_40_DoubleJet40_Mass_Min620
  - HLT_VBF_DiPFJet50_Mjj750_Photon22
* - L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28
  - HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60
* - L1_DoubleJet_120_45_DoubleJet45_Mass_Min620
  - HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60
* - L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28
  - HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85
* - L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ
  - HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85
* - L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp
  - HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL
* - L1_DoubleJet_80_30_Mass_Min420_Mu8
  - HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL
* - L1_DoubleJet_90_30_DoubleJet30_Mass_Min620
  - HLT_VBF_DoublePNetTauhPFJet20_eta2p2
* - L1_DoubleMu0er2p0_SQ_dR_Max1p4
  - HTXS_dPhijj
* - L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4
  - HTXS_Mjj
* - L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8
  - HTXS_ptHjj
* - L1_DoubleMu3_SQ_HTT240er
  - HTXS_V_pt
* - L1_DoubleMu3_SQ_HTT260er
  - Jet_btagPNetCvNotB
* - L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4
  - Jet_btagUParTAK4B
* - L1_ETMHF110_HTT60er_NotSecondBunchInTrain
  - Jet_btagUParTAK4CvB
* - L1_ETMHF120_NotSecondBunchInTrain
  - Jet_btagUParTAK4CvL
* - L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1
  - Jet_btagUParTAK4CvNotB
* - L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6
  - Jet_btagUParTAK4Ele
* - L1_ETT1200
  - Jet_btagUParTAK4Mu
* - L1_ETT1600
  - Jet_btagUParTAK4probb
* - L1_LooseIsoEG30er2p1_HTT100er
  - Jet_btagUParTAK4probbb
* - L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6
  - Jet_btagUParTAK4QvG
* - L1_Mu22er2p1_IsoTau28er2p1
  - Jet_btagUParTAK4SvCB
* - L1_Mu22er2p1_IsoTau36er2p1
  - Jet_btagUParTAK4SvUDG
* - L1_Mu3_Jet120er2p5_dR_Max0p8
  - Jet_btagUParTAK4TauVJet
* - L1_Mu3_Jet35er2p5_dR_Max0p4
  - Jet_btagUParTAK4UDG
* - L1_Mu3_Jet80er2p5_dR_Max0p4
  - Jet_chMultiplicity
* - L1_SingleIsoEG24er1p5
  - Jet_hfEmEF
* - L1_SingleIsoEG26er1p5
  - Jet_hfHEF
* - L1_SingleIsoTau32er2p1
  - Jet_muonSubtrDeltaEta
* - L1_SingleJet140er2p5_ETMHF70
  - Jet_muonSubtrDeltaPhi
* - L1_SingleJet140er2p5_ETMHF80
  - Jet_neMultiplicity
* - L1_SingleJet140er2p5_ETMHF90
  - Jet_puIdDisc
* - L1_SingleJet60_FWD3p0
  - Jet_UParTAK4RegPtRawCorr
* - L1_SingleJet60er2p5
  - Jet_UParTAK4RegPtRawCorrNeutrino
* - L1_SingleJet90_FWD3p0
  - Jet_UParTAK4RegPtRawRes
* - L1_SingleJet90er2p5
  - Jet_UParTAK4V1RegPtRawCorr
* - L1_SingleMu10er1p5
  - Jet_UParTAK4V1RegPtRawCorrNeutrino
* - L1_SingleMu12er1p5
  - Jet_UParTAK4V1RegPtRawRes
* - L1_SingleMu14er1p5
  - L1_AXO_Loose
* - L1_SingleMu16er1p5
  - L1_AXO_Nominal
* - L1_SingleMu18er1p5
  - L1_AXO_Tight
* - L1_SingleMu6er1p5
  - L1_AXO_VLoose
* - L1_SingleMu7er1p5
  - L1_AXO_VTight
* - L1_SingleMu8er1p5
  - L1_CICADA_Loose
* - L1_SingleMu9er1p5
  - L1_CICADA_Medium
* - L1_SingleTau70er2p1
  - L1_CICADA_Tight
* - L1_TripleEG16er2p5
  - L1_CICADA_VLoose
* - L1_TripleEG_16_12_8_er2p5
  - L1_CICADA_VTight
* - L1_TripleEG_16_15_8_er2p5
  - L1_DoubleIsoTau32er2p1_Mass_Max80
* - L1_TripleMu_2SQ_1p5SQ_0OQ
  - L1_DoubleJet120er2p5_Mu3_dR_Max0p8
* - L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12
  - L1_DoubleJet16er2p5_Mu3_dR_Max0p4
* - L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12
  - L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5
* - L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17
  - L1_DoubleJet35er2p5_Mu3_dR_Max0p4
* - L1_UnprefireableEvent
  - L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5
* - MET_fiducialGenPhi
  - L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2
* - MET_fiducialGenPt
  - L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5
* - MET_sumPtUnclustered
  - L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2
* - Muon_mvaTTH
  - L1_DoubleJet60er2p5_Mu3_dR_Max0p4
* - PuppiMET_phiJERDown
  - L1_DoubleJet80er2p5_Mu3_dR_Max0p4
* - PuppiMET_phiJERUp
  - L1_DoubleJet_110_35_DoubleJet35_Mass_Min800
* - PuppiMET_phiJESDown
  - L1_DoubleJet_110_35_DoubleJet35_Mass_Min850
* - PuppiMET_phiJESUp
  - L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50
* - PuppiMET_ptJERDown
  - L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50
* - PuppiMET_ptJERUp
  - L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65
* - PuppiMET_ptJESDown
  - L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65
* - PuppiMET_ptJESUp
  - L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ
* - SubJet_btagDeepB
  - L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ
* - Tau_idDeepTau2017v2p1VSe
  - L1_DoubleMu0_Upt6_SQ_er2p0
* - Tau_idDeepTau2017v2p1VSjet
  - L1_DoubleMu0_Upt7_SQ_er2p0
* - Tau_idDeepTau2017v2p1VSmu
  - L1_DoubleMu0_Upt8_SQ_er2p0
* - Tau_rawDeepTau2017v2p1VSe
  - L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2
* - Tau_rawDeepTau2017v2p1VSjet
  - L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2
* - Tau_rawDeepTau2017v2p1VSmu
  - L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2
* - 
  - L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6
* - 
  - L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6
* - 
  - L1_DoubleMu5_SQ_OS_dR_Max1p6
* - 
  - L1_DoubleMu6_Upt6_SQ_er2p0
* - 
  - L1_DoubleMu7_Upt7_SQ_er2p0
* - 
  - L1_DoubleMu8_Upt8_SQ_er2p0
* - 
  - L1_FinalOR_BXmin1
* - 
  - L1_FinalOR_BXmin2
* - 
  - L1_HTMHF100
* - 
  - L1_HTMHF120
* - 
  - L1_HTMHF125
* - 
  - L1_HTMHF130
* - 
  - L1_HTMHF150
* - 
  - L1_LooseIsoEG14er2p5_HTT200er
* - 
  - L1_LooseIsoEG16er2p5_HTT200er
* - 
  - L1_Mu12_HTT150er
* - 
  - L1_Mu14_HTT150er
* - 
  - L1_SingleJet120_FWD2p5
* - 
  - L1_SingleJet120er1p3
* - 
  - L1_SingleJet35_FWD2p5
* - 
  - L1_SingleJet35er1p3
* - 
  - L1_SingleJet60_FWD2p5
* - 
  - L1_SingleJet90_FWD2p5
* - 
  - L1_SingleMu0_SQ13_BMTF
* - 
  - L1_SingleMu0_SQ14_BMTF
* - 
  - L1_SingleMu0_SQ15_BMTF
* - 
  - L1_SingleMu0_Upt10
* - 
  - L1_SingleMu0_Upt10_BMTF
* - 
  - L1_SingleMu0_Upt10_EMTF
* - 
  - L1_SingleMu0_Upt10_OMTF
* - 
  - L1_SingleMu0_Upt10_SQ14_BMTF
* - 
  - L1_SingleMu0_Upt15_SQ14_BMTF
* - 
  - L1_SingleMu0_Upt20_SQ14_BMTF
* - 
  - L1_SingleMu0_Upt25_SQ14_BMTF
* - 
  - L1_SingleMu10_SQ14_BMTF
* - 
  - L1_SingleMu11_SQ14_BMTF
* - 
  - L1_SingleMu22_BMTF_NEG
* - 
  - L1_SingleMu22_BMTF_POS
* - 
  - L1_SingleMu22_EMTF_NEG
* - 
  - L1_SingleMu22_EMTF_POS
* - 
  - L1_SingleMu22_OMTF_NEG
* - 
  - L1_SingleMu22_OMTF_POS
* - 
  - L1_SingleMu5_SQ14_BMTF
* - 
  - L1_SingleMu6_SQ14_BMTF
* - 
  - L1_SingleMu7_SQ14_BMTF
* - 
  - L1_SingleMu8_SQ14_BMTF
* - 
  - L1_SingleMu9_SQ14_BMTF
* - 
  - L1_SingleMuOpen_BMTF
* - 
  - L1_SingleMuOpen_EMTF
* - 
  - L1_SingleMuOpen_OMTF
* - 
  - L1_TripleMu_3SQ_2p5SQ_0
* - 
  - L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12
* - 
  - L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12
* - 
  - L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12
* - 
  - L1_TwoMuShower_Loose
* - 
  - L1_UnprefireableEvent_FirstBxInTrain
* - 
  - L1_UnprefireableEvent_TriggerRules
* - 
  - LHEPart_firstMotherIdx
* - 
  - LHEPart_lastMotherIdx
* - 
  - MC_PFScouting
* - 
  - Muon_bestTrackType
* - 
  - Muon_dxybsErr
* - 
  - Muon_ipLengthSig
* - 
  - Muon_IPx
* - 
  - Muon_IPy
* - 
  - Muon_IPz
* - 
  - Muon_jetDF
* - 
  - Muon_pnScore_heavy
* - 
  - Muon_pnScore_light
* - 
  - Muon_pnScore_prompt
* - 
  - Muon_pnScore_tau
* - 
  - Muon_promptMVA
* - 
  - Muon_softMvaRun3
* - 
  - Muon_tuneP_charge
* - 
  - Muon_tuneP_pterr
* - 
  - Muon_VXBS_Cov00
* - 
  - Muon_VXBS_Cov03
* - 
  - Muon_VXBS_Cov33
* - 
  - nFatJetPFCand
* - 
  - nPFCand
* - 
  - nPVBS
* - 
  - nTauProd
* - 
  - nTrackGenJetAK4
* - 
  - orbitNumber
* - 
  - PFCand_eta
* - 
  - PFCand_mass
* - 
  - PFCand_pdgId
* - 
  - PFCand_phi
* - 
  - PFCand_pt
* - 
  - PFMET_phiUnclusteredDown
* - 
  - PFMET_phiUnclusteredUp
* - 
  - PFMET_ptUnclusteredDown
* - 
  - PFMET_ptUnclusteredUp
* - 
  - Photon_hoe_Tower
* - 
  - Photon_superclusterEta
* - 
  - Pileup_pthatmax
* - 
  - PuppiMET_covXX
* - 
  - PuppiMET_covXY
* - 
  - PuppiMET_covYY
* - 
  - PuppiMET_significance
* - 
  - PuppiMET_sumPtUnclustered
* - 
  - PV_sumpt2
* - 
  - PV_sumpx
* - 
  - PV_sumpy
* - 
  - PVBS_chi2
* - 
  - PVBS_cov00
* - 
  - PVBS_cov10
* - 
  - PVBS_cov11
* - 
  - PVBS_cov20
* - 
  - PVBS_cov21
* - 
  - PVBS_cov22
* - 
  - PVBS_x
* - 
  - PVBS_y
* - 
  - PVBS_z
* - 
  - SubJet_area
* - 
  - SubJet_btagDeepFlavB
* - 
  - SubJet_btagUParTAK4B
* - 
  - SubJet_subGenJetAK8Idx
* - 
  - SubJet_UParTAK4RegPtRawCorr
* - 
  - SubJet_UParTAK4RegPtRawCorrNeutrino
* - 
  - SubJet_UParTAK4RegPtRawRes
* - 
  - SubJet_UParTAK4V1RegPtRawCorr
* - 
  - SubJet_UParTAK4V1RegPtRawCorrNeutrino
* - 
  - SubJet_UParTAK4V1RegPtRawRes
* - 
  - Tau_decayModeUParT
* - 
  - Tau_hasRefitSV
* - 
  - Tau_ipLengthSig
* - 
  - Tau_IPx
* - 
  - Tau_IPy
* - 
  - Tau_IPz
* - 
  - Tau_probDM0UParT
* - 
  - Tau_probDM10UParT
* - 
  - Tau_probDM11UParT
* - 
  - Tau_probDM1UParT
* - 
  - Tau_probDM2UParT
* - 
  - Tau_ptCorrUParT
* - 
  - Tau_qConfUParT
* - 
  - Tau_rawUParTVSe
* - 
  - Tau_rawUParTVSjet
* - 
  - Tau_rawUParTVSmu
* - 
  - Tau_refitSVchi2
* - 
  - Tau_refitSVcov00
* - 
  - Tau_refitSVcov10
* - 
  - Tau_refitSVcov11
* - 
  - Tau_refitSVcov20
* - 
  - Tau_refitSVcov21
* - 
  - Tau_refitSVcov22
* - 
  - Tau_refitSVx
* - 
  - Tau_refitSVy
* - 
  - Tau_refitSVz
* - 
  - TauProd_eta
* - 
  - TauProd_pdgId
* - 
  - TauProd_phi
* - 
  - TauProd_pt
* - 
  - TauProd_tauIdx
* - 
  - TauSpinner_weight_cp_0
* - 
  - TauSpinner_weight_cp_0_alt
* - 
  - TauSpinner_weight_cp_0p25
* - 
  - TauSpinner_weight_cp_0p25_alt
* - 
  - TauSpinner_weight_cp_0p375
* - 
  - TauSpinner_weight_cp_0p375_alt
* - 
  - TauSpinner_weight_cp_0p5
* - 
  - TauSpinner_weight_cp_0p5_alt
* - 
  - TauSpinner_weight_cp_minus0p25
* - 
  - TauSpinner_weight_cp_minus0p25_alt
* - 
  - TrackGenJetAK4_eta
* - 
  - TrackGenJetAK4_phi
* - 
  - TrackGenJetAK4_pt
