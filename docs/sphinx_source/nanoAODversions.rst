NanoAOD versions
=================

CMS has different versions of the nanoAOD format. Usually in each data taking run an new version is introduced that can have new features or bug fixes compared to the previous version. The UL version of Run-2 data was fully processed in nanoAODv9, for Run-3 each year a new version is used starting with nanoAODv12. A preliminary final version for both Run-2 and Run-3 will be nanoAODv15. The changes between version also affect the types of the branches in the nanoAOD which can lead to problems in the C++ code of CROWN because the types need to be specified. In the following a list of known changes between nanoAOD v9 and v12 is given.

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