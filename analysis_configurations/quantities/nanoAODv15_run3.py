from code_generation.quantity import NanoAODQuantity

BeamSpot_sigmaZ = NanoAODQuantity("BeamSpot_sigmaZ")
"""dtype: Float_t; description: Width of BeamSpot in z (cm) """
BeamSpot_sigmaZError = NanoAODQuantity("BeamSpot_sigmaZError")
"""dtype: Float_t; description: Error on width of BeamSpot in z (cm) """
BeamSpot_type = NanoAODQuantity("BeamSpot_type")
"""dtype: Short_t; description: BeamSpot type (Unknown = -1, Fake = 0, LHC = 1, Tracker = 2) """
BeamSpot_z = NanoAODQuantity("BeamSpot_z")
"""dtype: Float_t; description: BeamSpot center, z coordinate (cm) """
BeamSpot_zError = NanoAODQuantity("BeamSpot_zError")
"""dtype: Float_t; description: Error on BeamSpot center, z coordinate (cm) """

CaloMET_phi = NanoAODQuantity("CaloMET_phi")
"""dtype: Float_t; description: phi """
CaloMET_pt = NanoAODQuantity("CaloMET_pt")
"""dtype: Float_t; description: pt """
CaloMET_sumEt = NanoAODQuantity("CaloMET_sumEt")
"""dtype: Float_t; description: scalar sum of Et """

nCorrT1METJet = NanoAODQuantity("nCorrT1METJet")
"""dtype: Int_t; description: Additional low-pt ak4 Puppi jets for Type-1 MET re-correction """
CorrT1METJet_EmEF = NanoAODQuantity("CorrT1METJet_EmEF")
"""dtype: Float_t; description: charged+neutral Electromagnetic Energy Fraction """
CorrT1METJet_area = NanoAODQuantity("CorrT1METJet_area")
"""dtype: Float_t; description: jet catchment area, for JECs """
CorrT1METJet_eta = NanoAODQuantity("CorrT1METJet_eta")
"""dtype: Float_t; description: eta """
CorrT1METJet_muonSubtrDeltaEta = NanoAODQuantity("CorrT1METJet_muonSubtrDeltaEta")
"""dtype: Float_t; description: muon-subtracted raw eta - eta """
CorrT1METJet_muonSubtrDeltaPhi = NanoAODQuantity("CorrT1METJet_muonSubtrDeltaPhi")
"""dtype: Float_t; description: muon-subtracted raw phi - phi """
CorrT1METJet_muonSubtrFactor = NanoAODQuantity("CorrT1METJet_muonSubtrFactor")
"""dtype: Float_t; description: 1-(muon-subtracted raw pt)/(raw pt) """
CorrT1METJet_phi = NanoAODQuantity("CorrT1METJet_phi")
"""dtype: Float_t; description: phi """
CorrT1METJet_rawMass = NanoAODQuantity("CorrT1METJet_rawMass")
"""dtype: Float_t; description: mass()*jecFactor('Uncorrected') """
CorrT1METJet_rawPt = NanoAODQuantity("CorrT1METJet_rawPt")
"""dtype: Float_t; description: pt()*jecFactor('Uncorrected') """

DST_PFScouting_AXOLoose = NanoAODQuantity("DST_PFScouting_AXOLoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_AXONominal = NanoAODQuantity("DST_PFScouting_AXONominal")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_AXOTight = NanoAODQuantity("DST_PFScouting_AXOTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_AXOVLoose = NanoAODQuantity("DST_PFScouting_AXOVLoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_AXOVTight = NanoAODQuantity("DST_PFScouting_AXOVTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_CICADALoose = NanoAODQuantity("DST_PFScouting_CICADALoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_CICADAMedium = NanoAODQuantity("DST_PFScouting_CICADAMedium")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_CICADATight = NanoAODQuantity("DST_PFScouting_CICADATight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_CICADAVLoose = NanoAODQuantity("DST_PFScouting_CICADAVLoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_CICADAVTight = NanoAODQuantity("DST_PFScouting_CICADAVTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_DatasetMuon = NanoAODQuantity("DST_PFScouting_DatasetMuon")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_DoubleEG = NanoAODQuantity("DST_PFScouting_DoubleEG")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_DoubleMuon = NanoAODQuantity("DST_PFScouting_DoubleMuon")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_JetHT = NanoAODQuantity("DST_PFScouting_JetHT")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_SingleMuon = NanoAODQuantity("DST_PFScouting_SingleMuon")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_SinglePhotonEB = NanoAODQuantity("DST_PFScouting_SinglePhotonEB")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
DST_PFScouting_ZeroBias = NanoAODQuantity("DST_PFScouting_ZeroBias")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

Dataset_ScoutingPFMonitor = NanoAODQuantity("Dataset_ScoutingPFMonitor")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
Dataset_ScoutingPFRun3 = NanoAODQuantity("Dataset_ScoutingPFRun3")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

DeepMETResolutionTune_phi = NanoAODQuantity("DeepMETResolutionTune_phi")
"""dtype: Float_t; description: DeepmET ResolutionTune phi """
DeepMETResolutionTune_pt = NanoAODQuantity("DeepMETResolutionTune_pt")
"""dtype: Float_t; description: DeepMET ResolutionTune pt """

DeepMETResponseTune_phi = NanoAODQuantity("DeepMETResponseTune_phi")
"""dtype: Float_t; description: DeepMET ResponseTune phi """
DeepMETResponseTune_pt = NanoAODQuantity("DeepMETResponseTune_pt")
"""dtype: Float_t; description: DeepMET ResponseTune pt """

nElectron = NanoAODQuantity("nElectron")
"""dtype: Int_t; description: slimmedElectrons after basic selection (pt > 5 ) """
Electron_IPx = NanoAODQuantity("Electron_IPx")
"""dtype: Float_t; description: x coordinate of impact parameter vector """
Electron_IPy = NanoAODQuantity("Electron_IPy")
"""dtype: Float_t; description: y coordinate of impact parameter vector """
Electron_IPz = NanoAODQuantity("Electron_IPz")
"""dtype: Float_t; description: z coordinate of impact parameter vector """
Electron_PreshowerEnergy = NanoAODQuantity("Electron_PreshowerEnergy")
"""dtype: Float_t; description: energy deposited in preshower """
Electron_charge = NanoAODQuantity("Electron_charge")
"""dtype: Int_t; description: electric charge """
Electron_convVeto = NanoAODQuantity("Electron_convVeto")
"""dtype: Bool_t; description: pass conversion veto """
Electron_cutBased = NanoAODQuantity("Electron_cutBased")
"""dtype: UChar_t; description: cut-based ID RunIII Winter22: fail ==0, veto >=1 (to veto, ask for <1), loose >=2, medium >=3, tight >=4 """
Electron_cutBased_HEEP = NanoAODQuantity("Electron_cutBased_HEEP")
"""dtype: Bool_t; description: cut-based HEEP ID """
Electron_deltaEtaSC = NanoAODQuantity("Electron_deltaEtaSC")
"""dtype: Float_t; description: delta eta (SC,ele) with sign """
Electron_dr03EcalRecHitSumEt = NanoAODQuantity("Electron_dr03EcalRecHitSumEt")
"""dtype: Float_t; description: Non-PF Ecal isolation within a delta R cone of 0.3 with electron pt > 35 GeV """
Electron_dr03HcalDepth1TowerSumEt = NanoAODQuantity("Electron_dr03HcalDepth1TowerSumEt")
"""dtype: Float_t; description: Non-PF Hcal isolation within a delta R cone of 0.3 with electron pt > 35 GeV """
Electron_dr03TkSumPt = NanoAODQuantity("Electron_dr03TkSumPt")
"""dtype: Float_t; description: Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV """
Electron_dr03TkSumPtHEEP = NanoAODQuantity("Electron_dr03TkSumPtHEEP")
"""dtype: Float_t; description: Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV used in HEEP ID """
Electron_dxy = NanoAODQuantity("Electron_dxy")
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm """
Electron_dxyErr = NanoAODQuantity("Electron_dxyErr")
"""dtype: Float_t; description: dxy uncertainty, in cm """
Electron_dz = NanoAODQuantity("Electron_dz")
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm """
Electron_dzErr = NanoAODQuantity("Electron_dzErr")
"""dtype: Float_t; description: dz uncertainty, in cm """
Electron_eInvMinusPInv = NanoAODQuantity("Electron_eInvMinusPInv")
"""dtype: Float_t; description: 1/E_SC - 1/p_trk """
Electron_ecalEnergy = NanoAODQuantity("Electron_ecalEnergy")
"""dtype: Float_t; description: energy after ECAL-only regression applied """
Electron_ecalEnergyError = NanoAODQuantity("Electron_ecalEnergyError")
"""dtype: Float_t; description: ecalEnergy error """
Electron_energyErr = NanoAODQuantity("Electron_energyErr")
"""dtype: Float_t; description: energy error of the cluster-track combination """
Electron_eta = NanoAODQuantity("Electron_eta")
"""dtype: Float_t; description: eta """
Electron_fbrem = NanoAODQuantity("Electron_fbrem")
"""dtype: Float_t; description: Fraction of brem """
Electron_fsrPhotonIdx = NanoAODQuantity("Electron_fsrPhotonIdx")
"""dtype: Short_t; description: Index of the lowest-dR/ET2 among associated FSR photons """
Electron_genPartFlav = NanoAODQuantity("Electron_genPartFlav")
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 electrons or photons: 1 = prompt electron (including gamma*->mu mu), 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched """
Electron_genPartIdx = NanoAODQuantity("Electron_genPartIdx")
"""dtype: Short_t; description: Index into genParticle list for MC matching to status==1 electrons or photons """
Electron_gsfTrketaMode = NanoAODQuantity("Electron_gsfTrketaMode")
"""dtype: Float_t; description: GSF track etaMode """
Electron_gsfTrkpMode = NanoAODQuantity("Electron_gsfTrkpMode")
"""dtype: Float_t; description: GSF track pMode """
Electron_gsfTrkpModeErr = NanoAODQuantity("Electron_gsfTrkpModeErr")
"""dtype: Float_t; description: GSF track pMode error """
Electron_gsfTrkphiMode = NanoAODQuantity("Electron_gsfTrkphiMode")
"""dtype: Float_t; description: GSF track phiMode """
Electron_hoe = NanoAODQuantity("Electron_hoe")
"""dtype: Float_t; description: H over E """
Electron_ip3d = NanoAODQuantity("Electron_ip3d")
"""dtype: Float_t; description: 3D impact parameter wrt first PV, in cm """
Electron_ipLengthSig = NanoAODQuantity("Electron_ipLengthSig")
"""dtype: Float_t; description: significance of impact parameter """
Electron_isEB = NanoAODQuantity("Electron_isEB")
"""dtype: Bool_t; description: object in barrel if true derived from the seedCrystal and detID information """
Electron_isEcalDriven = NanoAODQuantity("Electron_isEcalDriven")
"""dtype: Bool_t; description: is ECAL driven if true """
Electron_isPFcand = NanoAODQuantity("Electron_isPFcand")
"""dtype: Bool_t; description: electron is PF candidate """
Electron_jetDF = NanoAODQuantity("Electron_jetDF")
"""dtype: Float_t; description: value of the DEEPJET b tagging algorithm discriminator of the associated jet (0 if none) """
Electron_jetIdx = NanoAODQuantity("Electron_jetIdx")
"""dtype: Short_t; description: index of the associated jet (-1 if none) """
Electron_jetNDauCharged = NanoAODQuantity("Electron_jetNDauCharged")
"""dtype: UChar_t; description: number of charged daughters of the closest jet """
Electron_jetPtRelv2 = NanoAODQuantity("Electron_jetPtRelv2")
"""dtype: Float_t; description: Relative momentum of the lepton with respect to the closest jet after subtracting the lepton """
Electron_jetRelIso = NanoAODQuantity("Electron_jetRelIso")
"""dtype: Float_t; description: Relative isolation in matched jet (1/ptRatio-1), -1 if none """
Electron_lostHits = NanoAODQuantity("Electron_lostHits")
"""dtype: UChar_t; description: number of missing inner hits """
Electron_mass = NanoAODQuantity("Electron_mass")
"""dtype: Float_t; description: mass """
Electron_mvaHZZIso = NanoAODQuantity("Electron_mvaHZZIso")
"""dtype: Float_t; description: HZZ MVA Iso ID score """
Electron_pdgId = NanoAODQuantity("Electron_pdgId")
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth) """
Electron_pfRelIso04_all = NanoAODQuantity("Electron_pfRelIso04_all")
"""dtype: Float_t; description: PF relative isolation dR=0.4, total (with rho*EA PU Winter22V1 corrections) """
Electron_phi = NanoAODQuantity("Electron_phi")
"""dtype: Float_t; description: phi """
Electron_photonIdx = NanoAODQuantity("Electron_photonIdx")
"""dtype: Short_t; description: index of the first associated photon (-1 if none) """
Electron_promptMVA = NanoAODQuantity("Electron_promptMVA")
"""dtype: Float_t; description: Prompt MVA lepton ID score. Corresponds to the previous mvaTTH """
Electron_pt = NanoAODQuantity("Electron_pt")
"""dtype: Float_t; description: pt """
Electron_r9 = NanoAODQuantity("Electron_r9")
"""dtype: Float_t; description: R9 of the supercluster, calculated with full 5x5 region """
Electron_rawEnergy = NanoAODQuantity("Electron_rawEnergy")
"""dtype: Float_t; description: raw energy of Supercluster """
Electron_scEtOverPt = NanoAODQuantity("Electron_scEtOverPt")
"""dtype: Float_t; description: (supercluster transverse energy)/pt-1 """
Electron_seedGain = NanoAODQuantity("Electron_seedGain")
"""dtype: UChar_t; description: Gain of the seed crystal """
Electron_seediEtaOriX = NanoAODQuantity("Electron_seediEtaOriX")
"""dtype: Short_t; description: iEta or iX of seed crystal. iEta is barrel-only, iX is endcap-only. iEta runs from -85 to +85, with no crystal at iEta=0. iX runs from 1 to 100. """
Electron_seediPhiOriY = NanoAODQuantity("Electron_seediPhiOriY")
"""dtype: Short_t; description: iPhi or iY of seed crystal. iPhi is barrel-only, iY is endcap-only. iPhi runs from 1 to 360. iY runs from 1 to 100. """
Electron_sieie = NanoAODQuantity("Electron_sieie")
"""dtype: Float_t; description: sigma_IetaIeta of the supercluster, calculated with full 5x5 region """
Electron_sip3d = NanoAODQuantity("Electron_sip3d")
"""dtype: Float_t; description: 3D impact parameter significance wrt first PV, in cm """
Electron_superclusterEta = NanoAODQuantity("Electron_superclusterEta")
"""dtype: Float_t; description: supercluster eta """
Electron_svIdx = NanoAODQuantity("Electron_svIdx")
"""dtype: Short_t; description: index of matching secondary vertex """
Electron_tightCharge = NanoAODQuantity("Electron_tightCharge")
"""dtype: UChar_t; description: Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent) """
Electron_vidNestedWPBitmap = NanoAODQuantity("Electron_vidNestedWPBitmap")
"""dtype: Int_t; description: VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleEBEECut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEBEECut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 3 bits per cut """
Electron_vidNestedWPBitmapHEEP = NanoAODQuantity("Electron_vidNestedWPBitmapHEEP")
"""dtype: Int_t; description: VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleFull5x5SigmaIEtaIEtaWithSatCut,GsfEleFull5x5E2x5OverE5x5WithSatCut,GsfEleHadronicOverEMLinearCut,GsfEleTrkPtIsoCut,GsfEleEmHadD1IsoRhoCut,GsfEleDxyCut,GsfEleMissingHitsCut,GsfEleEcalDrivenCut), 1 bits per cut """

Electron_miniPFRelIso_all = NanoAODQuantity("Electron_miniPFRelIso_all")
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU Winter22V1 corrections) """
Electron_miniPFRelIso_chg = NanoAODQuantity("Electron_miniPFRelIso_chg")
"""dtype: Float_t; description: mini PF relative isolation, charged component """

Electron_mvaIso = NanoAODQuantity("Electron_mvaIso")
"""dtype: Float_t; description: MVA Iso ID score, Winter22V1 """
Electron_mvaIso_WP80 = NanoAODQuantity("Electron_mvaIso_WP80")
"""dtype: Bool_t; description: MVA Iso ID WP80, Winter22V1 """
Electron_mvaIso_WP90 = NanoAODQuantity("Electron_mvaIso_WP90")
"""dtype: Bool_t; description: MVA Iso ID WP90, Winter22V1 """
Electron_mvaIso_WPHZZ = NanoAODQuantity("Electron_mvaIso_WPHZZ")
"""dtype: Bool_t; description: MVA Iso ID WPHZZ, Winter22V1 """

Electron_mvaNoIso = NanoAODQuantity("Electron_mvaNoIso")
"""dtype: Float_t; description: MVA noIso ID score, Winter22V1 """
Electron_mvaNoIso_WP80 = NanoAODQuantity("Electron_mvaNoIso_WP80")
"""dtype: Bool_t; description: MVA noIso ID WP80, Winter22V1 """
Electron_mvaNoIso_WP90 = NanoAODQuantity("Electron_mvaNoIso_WP90")
"""dtype: Bool_t; description: MVA noIso ID WP90, Winter22V1 """

Electron_pfRelIso03_all = NanoAODQuantity("Electron_pfRelIso03_all")
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (with rho*EA PU Winter22V1 corrections) """
Electron_pfRelIso03_chg = NanoAODQuantity("Electron_pfRelIso03_chg")
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component """

nFatJet = NanoAODQuantity("nFatJet")
"""dtype: Int_t; description: slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis """
FatJet_area = NanoAODQuantity("FatJet_area")
"""dtype: Float_t; description: jet catchment area, for JECs """
FatJet_chEmEF = NanoAODQuantity("FatJet_chEmEF")
"""dtype: Float_t; description: charged Electromagnetic Energy Fraction """
FatJet_chHEF = NanoAODQuantity("FatJet_chHEF")
"""dtype: Float_t; description: charged Hadron Energy Fraction """
FatJet_chMultiplicity = NanoAODQuantity("FatJet_chMultiplicity")
"""dtype: Short_t; description: (Puppi-weighted) Number of charged particles in the jet """
FatJet_electronIdx3SJ = NanoAODQuantity("FatJet_electronIdx3SJ")
"""dtype: Short_t; description: index of electron matched to jet """
FatJet_eta = NanoAODQuantity("FatJet_eta")
"""dtype: Float_t; description: eta """
FatJet_genJetAK8Idx = NanoAODQuantity("FatJet_genJetAK8Idx")
"""dtype: Short_t; description: index of matched gen AK8 jet """
FatJet_hadronFlavour = NanoAODQuantity("FatJet_hadronFlavour")
"""dtype: UChar_t; description: flavour from hadron ghost clustering """
FatJet_hfEmEF = NanoAODQuantity("FatJet_hfEmEF")
"""dtype: Float_t; description: electromagnetic Energy Fraction in HF """
FatJet_hfHEF = NanoAODQuantity("FatJet_hfHEF")
"""dtype: Float_t; description: hadronic Energy Fraction in HF """
FatJet_lsf3 = NanoAODQuantity("FatJet_lsf3")
"""dtype: Float_t; description: Lepton Subjet Fraction (3 subjets) """
FatJet_mass = NanoAODQuantity("FatJet_mass")
"""dtype: Float_t; description: mass """
FatJet_msoftdrop = NanoAODQuantity("FatJet_msoftdrop")
"""dtype: Float_t; description: Corrected soft drop mass with PUPPI """
FatJet_muEF = NanoAODQuantity("FatJet_muEF")
"""dtype: Float_t; description: muon Energy Fraction """
FatJet_muonIdx3SJ = NanoAODQuantity("FatJet_muonIdx3SJ")
"""dtype: Short_t; description: index of muon matched to jet """
FatJet_n2b1 = NanoAODQuantity("FatJet_n2b1")
"""dtype: Float_t; description: N2 with beta=1 (for jets with raw pT>250 GeV) """
FatJet_n3b1 = NanoAODQuantity("FatJet_n3b1")
"""dtype: Float_t; description: N3 with beta=1 (for jets with raw pT>250 GeV) """
FatJet_nConstituents = NanoAODQuantity("FatJet_nConstituents")
"""dtype: UChar_t; description: Number of particles in the jet """
FatJet_neEmEF = NanoAODQuantity("FatJet_neEmEF")
"""dtype: Float_t; description: neutral Electromagnetic Energy Fraction """
FatJet_neHEF = NanoAODQuantity("FatJet_neHEF")
"""dtype: Float_t; description: neutral Hadron Energy Fraction """
FatJet_neMultiplicity = NanoAODQuantity("FatJet_neMultiplicity")
"""dtype: Short_t; description: (Puppi-weighted) Number of neutral particles in the jet """
FatJet_phi = NanoAODQuantity("FatJet_phi")
"""dtype: Float_t; description: phi """
FatJet_pt = NanoAODQuantity("FatJet_pt")
"""dtype: Float_t; description: pt """
FatJet_rawFactor = NanoAODQuantity("FatJet_rawFactor")
"""dtype: Float_t; description: 1 - Factor to get back to raw pT """
FatJet_subJetIdx1 = NanoAODQuantity("FatJet_subJetIdx1")
"""dtype: Short_t; description: index of first subjet """
FatJet_subJetIdx2 = NanoAODQuantity("FatJet_subJetIdx2")
"""dtype: Short_t; description: index of second subjet """
FatJet_tau1 = NanoAODQuantity("FatJet_tau1")
"""dtype: Float_t; description: Nsubjettiness (1 axis) """
FatJet_tau2 = NanoAODQuantity("FatJet_tau2")
"""dtype: Float_t; description: Nsubjettiness (2 axis) """
FatJet_tau3 = NanoAODQuantity("FatJet_tau3")
"""dtype: Float_t; description: Nsubjettiness (3 axis) """
FatJet_tau4 = NanoAODQuantity("FatJet_tau4")
"""dtype: Float_t; description: Nsubjettiness (4 axis) """

nFatJetPFCand = NanoAODQuantity("nFatJetPFCand")
"""dtype: Int_t; description:  """
FatJetPFCand_jetIdx = NanoAODQuantity("FatJetPFCand_jetIdx")
"""dtype: Int_t; description: Index of the parent jet """
FatJetPFCand_pfCandIdx = NanoAODQuantity("FatJetPFCand_pfCandIdx")
"""dtype: Int_t; description: Index in the PFCand table """

FatJet_globalParT3_QCD = NanoAODQuantity("FatJet_globalParT3_QCD")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 QCD score. """
FatJet_globalParT3_TopbWev = NanoAODQuantity("FatJet_globalParT3_TopbWev")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 Top->bWev score """
FatJet_globalParT3_TopbWmv = NanoAODQuantity("FatJet_globalParT3_TopbWmv")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 Top->bWmv score """
FatJet_globalParT3_TopbWq = NanoAODQuantity("FatJet_globalParT3_TopbWq")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 Top->bWq score """
FatJet_globalParT3_TopbWqq = NanoAODQuantity("FatJet_globalParT3_TopbWqq")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 Top->bWqq score """
FatJet_globalParT3_TopbWtauhv = NanoAODQuantity("FatJet_globalParT3_TopbWtauhv")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 Top->bWtauhv score """
FatJet_globalParT3_WvsQCD = NanoAODQuantity("FatJet_globalParT3_WvsQCD")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 (Xqq+Xcs/Xqq+Xcs+QCD) binarized score. """
FatJet_globalParT3_XWW3q = NanoAODQuantity("FatJet_globalParT3_XWW3q")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->WW3q score """
FatJet_globalParT3_XWW4q = NanoAODQuantity("FatJet_globalParT3_XWW4q")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->WW4q score """
FatJet_globalParT3_XWWqqev = NanoAODQuantity("FatJet_globalParT3_XWWqqev")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->WWqqev score """
FatJet_globalParT3_XWWqqmv = NanoAODQuantity("FatJet_globalParT3_XWWqqmv")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->WWqqmv score """
FatJet_globalParT3_Xbb = NanoAODQuantity("FatJet_globalParT3_Xbb")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->bb score. Note: For sig vs bkg (e.g. bkg=QCD) tagging, use sig/(sig+bkg) to construct the discriminator """
FatJet_globalParT3_Xcc = NanoAODQuantity("FatJet_globalParT3_Xcc")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->cc score """
FatJet_globalParT3_Xcs = NanoAODQuantity("FatJet_globalParT3_Xcs")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->cs score """
FatJet_globalParT3_Xqq = NanoAODQuantity("FatJet_globalParT3_Xqq")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->qq (ss/dd/uu) score """
FatJet_globalParT3_Xtauhtaue = NanoAODQuantity("FatJet_globalParT3_Xtauhtaue")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->tauhtaue score """
FatJet_globalParT3_Xtauhtauh = NanoAODQuantity("FatJet_globalParT3_Xtauhtauh")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->tauhtauh score """
FatJet_globalParT3_Xtauhtaum = NanoAODQuantity("FatJet_globalParT3_Xtauhtaum")
"""dtype: Float_t; description: Mass-decorrelated GlobalParT-3 X->tauhtaum score """
FatJet_globalParT3_massCorrGeneric = NanoAODQuantity(
    "FatJet_globalParT3_massCorrGeneric"
)
"""dtype: Float_t; description: GlobalParT-3 mass regression corrector with respect to the original jet mass, optimised for generic jet cases. Use (massCorrGeneric * mass * (1 - rawFactor)) to get the regressed mass """
FatJet_globalParT3_massCorrX2p = NanoAODQuantity("FatJet_globalParT3_massCorrX2p")
"""dtype: Float_t; description: GlobalParT-3 mass regression corrector with respect to the original jet mass, optimised for resonance 2-prong (bb/cc/cs/ss/qq) jets. Use (massCorrX2p * mass * (1 - rawFactor)) to get the regressed mass """
FatJet_globalParT3_withMassTopvsQCD = NanoAODQuantity(
    "FatJet_globalParT3_withMassTopvsQCD"
)
"""dtype: Float_t; description: GlobalParT-3 tagger (w/mass) Top vs QCD discriminator """
FatJet_globalParT3_withMassWvsQCD = NanoAODQuantity("FatJet_globalParT3_withMassWvsQCD")
"""dtype: Float_t; description: GlobalParT-3 tagger (w/mass) W vs QCD discriminator """
FatJet_globalParT3_withMassZvsQCD = NanoAODQuantity("FatJet_globalParT3_withMassZvsQCD")
"""dtype: Float_t; description: GlobalParT-3 tagger (w/mass) Z vs QCD discriminator """

FatJet_particleNet_QCD = NanoAODQuantity("FatJet_particleNet_QCD")
"""dtype: Float_t; description: ParticleNet tagger QCD(0+1+2HF) sum """
FatJet_particleNet_QCD0HF = NanoAODQuantity("FatJet_particleNet_QCD0HF")
"""dtype: Float_t; description: ParticleNet tagger QCD 0 HF (b/c) score """
FatJet_particleNet_QCD1HF = NanoAODQuantity("FatJet_particleNet_QCD1HF")
"""dtype: Float_t; description: ParticleNet tagger QCD 1 HF (b/c) score """
FatJet_particleNet_QCD2HF = NanoAODQuantity("FatJet_particleNet_QCD2HF")
"""dtype: Float_t; description: ParticleNet tagger QCD 2 HF (b/c) score """
FatJet_particleNet_WVsQCD = NanoAODQuantity("FatJet_particleNet_WVsQCD")
"""dtype: Float_t; description: ParticleNet W->qq vs. QCD score: Xqq+Xcc/(Xqq+Xcc+QCD) """
FatJet_particleNet_XbbVsQCD = NanoAODQuantity("FatJet_particleNet_XbbVsQCD")
"""dtype: Float_t; description: ParticleNet X->bb vs. QCD score: Xbb/(Xbb+QCD) """
FatJet_particleNet_XccVsQCD = NanoAODQuantity("FatJet_particleNet_XccVsQCD")
"""dtype: Float_t; description: ParticleNet X->cc vs. QCD score: Xcc/(Xcc+QCD) """
FatJet_particleNet_XggVsQCD = NanoAODQuantity("FatJet_particleNet_XggVsQCD")
"""dtype: Float_t; description: ParticleNet X->gg vs. QCD score: Xgg/(Xgg+QCD) """
FatJet_particleNet_XqqVsQCD = NanoAODQuantity("FatJet_particleNet_XqqVsQCD")
"""dtype: Float_t; description: ParticleNet X->qq (uds) vs. QCD score: Xqq/(Xqq+QCD) """
FatJet_particleNet_XteVsQCD = NanoAODQuantity("FatJet_particleNet_XteVsQCD")
"""dtype: Float_t; description: ParticleNet X->e tau_h vs. QCD score: Xte/(Xte+QCD) """
FatJet_particleNet_XtmVsQCD = NanoAODQuantity("FatJet_particleNet_XtmVsQCD")
"""dtype: Float_t; description: ParticleNet X->mu tau_h vs. QCD score: Xtm/(Xtm+QCD) """
FatJet_particleNet_XttVsQCD = NanoAODQuantity("FatJet_particleNet_XttVsQCD")
"""dtype: Float_t; description: ParticleNet X->tau_h tau_h vs. QCD score: Xtt/(Xtt+QCD) """
FatJet_particleNet_massCorr = NanoAODQuantity("FatJet_particleNet_massCorr")
"""dtype: Float_t; description: ParticleNet mass regression, relative correction to JEC-corrected jet mass (no softdrop) """

FatJet_particleNetLegacy_QCD = NanoAODQuantity("FatJet_particleNetLegacy_QCD")
"""dtype: Float_t; description: Mass-decorrelated ParticleNet Legacy Run-2 tagger raw QCD score """
FatJet_particleNetLegacy_Xbb = NanoAODQuantity("FatJet_particleNetLegacy_Xbb")
"""dtype: Float_t; description: Mass-decorrelated ParticleNet Legacy Run-2 tagger raw X->bb score. For X->bb vs QCD tagging, use Xbb/(Xbb+QCD) """
FatJet_particleNetLegacy_Xcc = NanoAODQuantity("FatJet_particleNetLegacy_Xcc")
"""dtype: Float_t; description: Mass-decorrelated ParticleNet Legacy Run-2 tagger raw X->cc score. For X->cc vs QCD tagging, use Xcc/(Xcc+QCD) """
FatJet_particleNetLegacy_Xqq = NanoAODQuantity("FatJet_particleNetLegacy_Xqq")
"""dtype: Float_t; description: Mass-decorrelated ParticleNet Legacy Run-2 tagger raw X->qq (uds) score. For X->qq vs QCD tagging, use Xqq/(Xqq+QCD). For W vs QCD tagging, use (Xcc+Xqq)/(Xcc+Xqq+QCD) """
FatJet_particleNetLegacy_mass = NanoAODQuantity("FatJet_particleNetLegacy_mass")
"""dtype: Float_t; description: ParticleNet Legacy Run-2 mass regression """

FatJet_particleNetWithMass_H4qvsQCD = NanoAODQuantity(
    "FatJet_particleNetWithMass_H4qvsQCD"
)
"""dtype: Float_t; description: ParticleNet tagger (w/ mass) H(->VV->qqqq) vs QCD discriminator """
FatJet_particleNetWithMass_HbbvsQCD = NanoAODQuantity(
    "FatJet_particleNetWithMass_HbbvsQCD"
)
"""dtype: Float_t; description: ParticleNet tagger (w/mass) H(->bb) vs QCD discriminator """
FatJet_particleNetWithMass_HccvsQCD = NanoAODQuantity(
    "FatJet_particleNetWithMass_HccvsQCD"
)
"""dtype: Float_t; description: ParticleNet tagger (w/mass) H(->cc) vs QCD discriminator """
FatJet_particleNetWithMass_QCD = NanoAODQuantity("FatJet_particleNetWithMass_QCD")
"""dtype: Float_t; description: ParticleNet tagger (w/ mass) QCD(bb,cc,b,c,others) sum """
FatJet_particleNetWithMass_TvsQCD = NanoAODQuantity("FatJet_particleNetWithMass_TvsQCD")
"""dtype: Float_t; description: ParticleNet tagger (w/ mass) top vs QCD discriminator """
FatJet_particleNetWithMass_WvsQCD = NanoAODQuantity("FatJet_particleNetWithMass_WvsQCD")
"""dtype: Float_t; description: ParticleNet tagger (w/ mass) W vs QCD discriminator """
FatJet_particleNetWithMass_ZvsQCD = NanoAODQuantity("FatJet_particleNetWithMass_ZvsQCD")
"""dtype: Float_t; description: ParticleNet tagger (w/ mass) Z vs QCD discriminator """

FiducialMET_phi = NanoAODQuantity("FiducialMET_phi")
"""dtype: Float_t; description: phi """
FiducialMET_pt = NanoAODQuantity("FiducialMET_pt")
"""dtype: Float_t; description: pt """

Flag_BadChargedCandidateFilter = NanoAODQuantity("Flag_BadChargedCandidateFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_BadChargedCandidateSummer16Filter = NanoAODQuantity(
    "Flag_BadChargedCandidateSummer16Filter"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_BadPFMuonDzFilter = NanoAODQuantity("Flag_BadPFMuonDzFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_BadPFMuonFilter = NanoAODQuantity("Flag_BadPFMuonFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_BadPFMuonSummer16Filter = NanoAODQuantity("Flag_BadPFMuonSummer16Filter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_CSCTightHalo2015Filter = NanoAODQuantity("Flag_CSCTightHalo2015Filter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_CSCTightHaloFilter = NanoAODQuantity("Flag_CSCTightHaloFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_CSCTightHaloTrkMuUnvetoFilter = NanoAODQuantity(
    "Flag_CSCTightHaloTrkMuUnvetoFilter"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_EcalDeadCellBoundaryEnergyFilter = NanoAODQuantity(
    "Flag_EcalDeadCellBoundaryEnergyFilter"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_EcalDeadCellTriggerPrimitiveFilter = NanoAODQuantity(
    "Flag_EcalDeadCellTriggerPrimitiveFilter"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_HBHENoiseFilter = NanoAODQuantity("Flag_HBHENoiseFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_HBHENoiseIsoFilter = NanoAODQuantity("Flag_HBHENoiseIsoFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_HcalStripHaloFilter = NanoAODQuantity("Flag_HcalStripHaloFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_chargedHadronTrackResolutionFilter = NanoAODQuantity(
    "Flag_chargedHadronTrackResolutionFilter"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_ecalBadCalibFilter = NanoAODQuantity("Flag_ecalBadCalibFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_ecalLaserCorrFilter = NanoAODQuantity("Flag_ecalLaserCorrFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_eeBadScFilter = NanoAODQuantity("Flag_eeBadScFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_globalSuperTightHalo2016Filter = NanoAODQuantity(
    "Flag_globalSuperTightHalo2016Filter"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_globalTightHalo2016Filter = NanoAODQuantity("Flag_globalTightHalo2016Filter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_goodVertices = NanoAODQuantity("Flag_goodVertices")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_hcalLaserEventFilter = NanoAODQuantity("Flag_hcalLaserEventFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_hfNoisyHitsFilter = NanoAODQuantity("Flag_hfNoisyHitsFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_muonBadTrackFilter = NanoAODQuantity("Flag_muonBadTrackFilter")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_trkPOGFilters = NanoAODQuantity("Flag_trkPOGFilters")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """

Flag_trkPOG_logErrorTooManyClusters = NanoAODQuantity(
    "Flag_trkPOG_logErrorTooManyClusters"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_trkPOG_manystripclus53X = NanoAODQuantity("Flag_trkPOG_manystripclus53X")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """
Flag_trkPOG_toomanystripclus53X = NanoAODQuantity("Flag_trkPOG_toomanystripclus53X")
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT) """

nFsrPhoton = NanoAODQuantity("nFsrPhoton")
"""dtype: Int_t; description: Final state radiation photons emitted by muons or electrons """
FsrPhoton_dROverEt2 = NanoAODQuantity("FsrPhoton_dROverEt2")
"""dtype: Float_t; description: deltaR to associated muon divided by photon et2 """
FsrPhoton_electronIdx = NanoAODQuantity("FsrPhoton_electronIdx")
"""dtype: Short_t; description: index of associated electron """
FsrPhoton_eta = NanoAODQuantity("FsrPhoton_eta")
"""dtype: Float_t; description: eta """
FsrPhoton_muonIdx = NanoAODQuantity("FsrPhoton_muonIdx")
"""dtype: Short_t; description: index of associated muon """
FsrPhoton_phi = NanoAODQuantity("FsrPhoton_phi")
"""dtype: Float_t; description: phi """
FsrPhoton_pt = NanoAODQuantity("FsrPhoton_pt")
"""dtype: Float_t; description: pt """
FsrPhoton_relIso03 = NanoAODQuantity("FsrPhoton_relIso03")
"""dtype: Float_t; description: relative isolation in a 0.3 cone without CHS """

nGenDressedLepton = NanoAODQuantity("nGenDressedLepton")
"""dtype: Int_t; description: Dressed leptons from Rivet-based ParticleLevelProducer """
GenDressedLepton_eta = NanoAODQuantity("GenDressedLepton_eta")
"""dtype: Float_t; description: eta """
GenDressedLepton_hasTauAnc = NanoAODQuantity("GenDressedLepton_hasTauAnc")
"""dtype: Bool_t; description: true if Dressed lepton has a tau as ancestor """
GenDressedLepton_mass = NanoAODQuantity("GenDressedLepton_mass")
"""dtype: Float_t; description: mass """
GenDressedLepton_pdgId = NanoAODQuantity("GenDressedLepton_pdgId")
"""dtype: Int_t; description: PDG id """
GenDressedLepton_phi = NanoAODQuantity("GenDressedLepton_phi")
"""dtype: Float_t; description: phi """
GenDressedLepton_pt = NanoAODQuantity("GenDressedLepton_pt")
"""dtype: Float_t; description: pt """

nGenIsolatedPhoton = NanoAODQuantity("nGenIsolatedPhoton")
"""dtype: Int_t; description: Isolated photons from Rivet-based ParticleLevelProducer """
GenIsolatedPhoton_eta = NanoAODQuantity("GenIsolatedPhoton_eta")
"""dtype: Float_t; description: eta """
GenIsolatedPhoton_mass = NanoAODQuantity("GenIsolatedPhoton_mass")
"""dtype: Float_t; description: mass """
GenIsolatedPhoton_phi = NanoAODQuantity("GenIsolatedPhoton_phi")
"""dtype: Float_t; description: phi """
GenIsolatedPhoton_pt = NanoAODQuantity("GenIsolatedPhoton_pt")
"""dtype: Float_t; description: pt """

nGenJet = NanoAODQuantity("nGenJet")
"""dtype: Int_t; description: slimmedGenJets, i.e. ak4 Jets made with visible genparticles """
GenJet_eta = NanoAODQuantity("GenJet_eta")
"""dtype: Float_t; description: eta """
GenJet_hadronFlavour = NanoAODQuantity("GenJet_hadronFlavour")
"""dtype: UChar_t; description: flavour from hadron ghost clustering """
GenJet_mass = NanoAODQuantity("GenJet_mass")
"""dtype: Float_t; description: mass """
GenJet_nBHadrons = NanoAODQuantity("GenJet_nBHadrons")
"""dtype: UChar_t; description: number of b-hadrons """
GenJet_nCHadrons = NanoAODQuantity("GenJet_nCHadrons")
"""dtype: UChar_t; description: number of c-hadrons """
GenJet_partonFlavour = NanoAODQuantity("GenJet_partonFlavour")
"""dtype: Short_t; description: flavour from parton matching """
GenJet_phi = NanoAODQuantity("GenJet_phi")
"""dtype: Float_t; description: phi """
GenJet_pt = NanoAODQuantity("GenJet_pt")
"""dtype: Float_t; description: pt """

nGenJetAK8 = NanoAODQuantity("nGenJetAK8")
"""dtype: Int_t; description: slimmedGenJetsAK8, i.e. ak8 Jets made with visible genparticles """
GenJetAK8_eta = NanoAODQuantity("GenJetAK8_eta")
"""dtype: Float_t; description: eta """
GenJetAK8_hadronFlavour = NanoAODQuantity("GenJetAK8_hadronFlavour")
"""dtype: UChar_t; description: flavour from hadron ghost clustering """
GenJetAK8_mass = NanoAODQuantity("GenJetAK8_mass")
"""dtype: Float_t; description: mass """
GenJetAK8_nBHadrons = NanoAODQuantity("GenJetAK8_nBHadrons")
"""dtype: UChar_t; description: number of b-hadrons """
GenJetAK8_nCHadrons = NanoAODQuantity("GenJetAK8_nCHadrons")
"""dtype: UChar_t; description: number of c-hadrons """
GenJetAK8_partonFlavour = NanoAODQuantity("GenJetAK8_partonFlavour")
"""dtype: Short_t; description: flavour from parton matching """
GenJetAK8_phi = NanoAODQuantity("GenJetAK8_phi")
"""dtype: Float_t; description: phi """
GenJetAK8_pt = NanoAODQuantity("GenJetAK8_pt")
"""dtype: Float_t; description: pt """

GenMET_phi = NanoAODQuantity("GenMET_phi")
"""dtype: Float_t; description: phi """
GenMET_pt = NanoAODQuantity("GenMET_pt")
"""dtype: Float_t; description: pt """

nGenPart = NanoAODQuantity("nGenPart")
"""dtype: Int_t; description: interesting gen particles  """
GenPart_eta = NanoAODQuantity("GenPart_eta")
"""dtype: Float_t; description: eta """
GenPart_genPartIdxMother = NanoAODQuantity("GenPart_genPartIdxMother")
"""dtype: Short_t; description: index of the mother particle """
GenPart_iso = NanoAODQuantity("GenPart_iso")
"""dtype: Float_t; description: Isolation for leptons """
GenPart_mass = NanoAODQuantity("GenPart_mass")
"""dtype: Float_t; description: Mass stored for all particles with the exception of quarks (except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411). For these particles, you can lookup the value from PDG. """
GenPart_pdgId = NanoAODQuantity("GenPart_pdgId")
"""dtype: Int_t; description: PDG id """
GenPart_phi = NanoAODQuantity("GenPart_phi")
"""dtype: Float_t; description: phi """
GenPart_pt = NanoAODQuantity("GenPart_pt")
"""dtype: Float_t; description: pt """
GenPart_status = NanoAODQuantity("GenPart_status")
"""dtype: Int_t; description: Particle status. 1=stable """
GenPart_statusFlags = NanoAODQuantity("GenPart_statusFlags")
"""dtype: UShort_t; description: gen status flags stored bitwise, bits are: 0 : isPrompt, 1 : isDecayedLeptonHadron, 2 : isTauDecayProduct, 3 : isPromptTauDecayProduct, 4 : isDirectTauDecayProduct, 5 : isDirectPromptTauDecayProduct, 6 : isDirectHadronDecayProduct, 7 : isHardProcess, 8 : fromHardProcess, 9 : isHardProcessTauDecayProduct, 10 : isDirectHardProcessTauDecayProduct, 11 : fromHardProcessBeforeFSR, 12 : isFirstCopy, 13 : isLastCopy, 14 : isLastCopyBeforeFSR,  """

nGenProton = NanoAODQuantity("nGenProton")
"""dtype: Int_t; description:  """
GenProton_isPU = NanoAODQuantity("GenProton_isPU")
"""dtype: Bool_t; description: pileup proton? """
GenProton_px = NanoAODQuantity("GenProton_px")
"""dtype: Float_t; description: proton horizontal momentum """
GenProton_py = NanoAODQuantity("GenProton_py")
"""dtype: Float_t; description: proton vertical momentum """
GenProton_pz = NanoAODQuantity("GenProton_pz")
"""dtype: Float_t; description: proton longitudinal momentum """
GenProton_vz = NanoAODQuantity("GenProton_vz")
"""dtype: Float_t; description: proton vertex longitudinal coordinate """

nGenVisTau = NanoAODQuantity("nGenVisTau")
"""dtype: Int_t; description: gen hadronic taus  """
GenVisTau_charge = NanoAODQuantity("GenVisTau_charge")
"""dtype: Short_t; description: charge """
GenVisTau_eta = NanoAODQuantity("GenVisTau_eta")
"""dtype: Float_t; description: eta """
GenVisTau_genPartIdxMother = NanoAODQuantity("GenVisTau_genPartIdxMother")
"""dtype: Short_t; description: index of the mother particle """
GenVisTau_mass = NanoAODQuantity("GenVisTau_mass")
"""dtype: Float_t; description: mass """
GenVisTau_phi = NanoAODQuantity("GenVisTau_phi")
"""dtype: Float_t; description: phi """
GenVisTau_pt = NanoAODQuantity("GenVisTau_pt")
"""dtype: Float_t; description: pt """
GenVisTau_status = NanoAODQuantity("GenVisTau_status")
"""dtype: UChar_t; description: Hadronic tau decay mode. 0=OneProng0PiZero, 1=OneProng1PiZero, 2=OneProng2PiZero, 10=ThreeProng0PiZero, 11=ThreeProng1PiZero, 15=Other """

GenVtx_t0 = NanoAODQuantity("GenVtx_t0")
"""dtype: Float_t; description: gen vertex t0 """
GenVtx_x = NanoAODQuantity("GenVtx_x")
"""dtype: Float_t; description: gen vertex x """
GenVtx_y = NanoAODQuantity("GenVtx_y")
"""dtype: Float_t; description: gen vertex y """
GenVtx_z = NanoAODQuantity("GenVtx_z")
"""dtype: Float_t; description: gen vertex z """

Generator_binvar = NanoAODQuantity("Generator_binvar")
"""dtype: Float_t; description: MC generation binning value """
Generator_id1 = NanoAODQuantity("Generator_id1")
"""dtype: Int_t; description: id of first parton """
Generator_id2 = NanoAODQuantity("Generator_id2")
"""dtype: Int_t; description: id of second parton """
Generator_scalePDF = NanoAODQuantity("Generator_scalePDF")
"""dtype: Float_t; description: Q2 scale for PDF """
Generator_weight = NanoAODQuantity("Generator_weight")
"""dtype: Float_t; description: MC generator weight """
Generator_x1 = NanoAODQuantity("Generator_x1")
"""dtype: Float_t; description: x1 fraction of proton momentum carried by the first parton """
Generator_x2 = NanoAODQuantity("Generator_x2")
"""dtype: Float_t; description: x2 fraction of proton momentum carried by the second parton """
Generator_xpdf1 = NanoAODQuantity("Generator_xpdf1")
"""dtype: Float_t; description: x*pdf(x) for the first parton """
Generator_xpdf2 = NanoAODQuantity("Generator_xpdf2")
"""dtype: Float_t; description: x*pdf(x) for the second parton """

HLT_AK8DiPFJet270_270_SoftDropMass30 = NanoAODQuantity(
    "HLT_AK8DiPFJet270_270_SoftDropMass30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8DiPFJet280_280_SoftDropMass30 = NanoAODQuantity(
    "HLT_AK8DiPFJet280_280_SoftDropMass30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8DiPFJet290_290_SoftDropMass30 = NanoAODQuantity(
    "HLT_AK8DiPFJet290_290_SoftDropMass30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet140 = NanoAODQuantity("HLT_AK8PFJet140")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet200 = NanoAODQuantity("HLT_AK8PFJet200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet260 = NanoAODQuantity("HLT_AK8PFJet260")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet320 = NanoAODQuantity("HLT_AK8PFJet320")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet380_SoftDropMass30 = NanoAODQuantity("HLT_AK8PFJet380_SoftDropMass30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet40 = NanoAODQuantity("HLT_AK8PFJet40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet400 = NanoAODQuantity("HLT_AK8PFJet400")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet400_SoftDropMass30 = NanoAODQuantity("HLT_AK8PFJet400_SoftDropMass30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet425_SoftDropMass30 = NanoAODQuantity("HLT_AK8PFJet425_SoftDropMass30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet450 = NanoAODQuantity("HLT_AK8PFJet450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet450_SoftDropMass30 = NanoAODQuantity("HLT_AK8PFJet450_SoftDropMass30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet500 = NanoAODQuantity("HLT_AK8PFJet500")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet550 = NanoAODQuantity("HLT_AK8PFJet550")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet60 = NanoAODQuantity("HLT_AK8PFJet60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet80 = NanoAODQuantity("HLT_AK8PFJet80")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd140 = NanoAODQuantity("HLT_AK8PFJetFwd140")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd200 = NanoAODQuantity("HLT_AK8PFJetFwd200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd260 = NanoAODQuantity("HLT_AK8PFJetFwd260")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd320 = NanoAODQuantity("HLT_AK8PFJetFwd320")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd40 = NanoAODQuantity("HLT_AK8PFJetFwd40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd400 = NanoAODQuantity("HLT_AK8PFJetFwd400")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd450 = NanoAODQuantity("HLT_AK8PFJetFwd450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd500 = NanoAODQuantity("HLT_AK8PFJetFwd500")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd60 = NanoAODQuantity("HLT_AK8PFJetFwd60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJetFwd80 = NanoAODQuantity("HLT_AK8PFJetFwd80")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CaloJet500_NoJetID = NanoAODQuantity("HLT_CaloJet500_NoJetID")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CaloJet550_NoJetID = NanoAODQuantity("HLT_CaloJet550_NoJetID")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CaloMET350_NotCleaned = NanoAODQuantity("HLT_CaloMET350_NotCleaned")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CaloMET90_NotCleaned = NanoAODQuantity("HLT_CaloMET90_NotCleaned")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CaloMHT90 = NanoAODQuantity("HLT_CaloMHT90")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CascadeMu100 = NanoAODQuantity("HLT_CascadeMu100")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiEle27_WPTightCaloOnly_L1DoubleEG = NanoAODQuantity(
    "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8 = NanoAODQuantity(
    "HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve100_HFJEC = NanoAODQuantity("HLT_DiPFJetAve100_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve140 = NanoAODQuantity("HLT_DiPFJetAve140")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve160_HFJEC = NanoAODQuantity("HLT_DiPFJetAve160_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP = NanoAODQuantity(
    "HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve200 = NanoAODQuantity("HLT_DiPFJetAve200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve220_HFJEC = NanoAODQuantity("HLT_DiPFJetAve220_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve260 = NanoAODQuantity("HLT_DiPFJetAve260")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve260_HFJEC = NanoAODQuantity("HLT_DiPFJetAve260_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve300_HFJEC = NanoAODQuantity("HLT_DiPFJetAve300_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve320 = NanoAODQuantity("HLT_DiPFJetAve320")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve40 = NanoAODQuantity("HLT_DiPFJetAve40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve400 = NanoAODQuantity("HLT_DiPFJetAve400")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve500 = NanoAODQuantity("HLT_DiPFJetAve500")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve60 = NanoAODQuantity("HLT_DiPFJetAve60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve60_HFJEC = NanoAODQuantity("HLT_DiPFJetAve60_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve80 = NanoAODQuantity("HLT_DiPFJetAve80")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPFJetAve80_HFJEC = NanoAODQuantity("HLT_DiPFJetAve80_HFJEC")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10Time1ns = NanoAODQuantity("HLT_DiPhoton10Time1ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10Time1p2ns = NanoAODQuantity("HLT_DiPhoton10Time1p2ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10Time1p4ns = NanoAODQuantity("HLT_DiPhoton10Time1p4ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10Time1p6ns = NanoAODQuantity("HLT_DiPhoton10Time1p6ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10Time1p8ns = NanoAODQuantity("HLT_DiPhoton10Time1p8ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10Time2ns = NanoAODQuantity("HLT_DiPhoton10Time2ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiPhoton10_CaloIdL = NanoAODQuantity("HLT_DiPhoton10_CaloIdL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiSC30_18_EIso_AND_HE_Mass70 = NanoAODQuantity("HLT_DiSC30_18_EIso_AND_HE_Mass70")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon10_Upsilon_y1p4 = NanoAODQuantity("HLT_Dimuon10_Upsilon_y1p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon12_Upsilon_y1p4 = NanoAODQuantity("HLT_Dimuon12_Upsilon_y1p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24 = NanoAODQuantity(
    "HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleCscCluster100 = NanoAODQuantity("HLT_DoubleCscCluster100")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleCscCluster75 = NanoAODQuantity("HLT_DoubleCscCluster75")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle10_eta1p22_mMax6 = NanoAODQuantity("HLT_DoubleEle10_eta1p22_mMax6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle24_eta2p1_WPTight_Gsf = NanoAODQuantity(
    "HLT_DoubleEle24_eta2p1_WPTight_Gsf"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle25_CaloIdL_MW = NanoAODQuantity("HLT_DoubleEle25_CaloIdL_MW")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle27_CaloIdL_MW = NanoAODQuantity("HLT_DoubleEle27_CaloIdL_MW")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle33_CaloIdL_MW = NanoAODQuantity("HLT_DoubleEle33_CaloIdL_MW")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle6p5_eta1p22_mMax6 = NanoAODQuantity("HLT_DoubleEle6p5_eta1p22_mMax6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleIsoMu20_eta2p1 = NanoAODQuantity("HLT_DoubleIsoMu20_eta2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu50 = NanoAODQuantity("HLT_DoubleL2Mu50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm = NanoAODQuantity(
    "HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm = NanoAODQuantity(
    "HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm = NanoAODQuantity(
    "HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm = NanoAODQuantity(
    "HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1 = NanoAODQuantity(
    "HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1 = NanoAODQuantity(
    "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu43NoFiltersNoVtx = NanoAODQuantity("HLT_DoubleMu43NoFiltersNoVtx")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu48NoFiltersNoVtx = NanoAODQuantity("HLT_DoubleMu48NoFiltersNoVtx")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL = NanoAODQuantity(
    "HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePFJets100_PNetBTag_0p11 = NanoAODQuantity("HLT_DoublePFJets100_PNetBTag_0p11")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11 = NanoAODQuantity(
    "HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11 = NanoAODQuantity(
    "HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePFJets200_PNetBTag_0p11 = NanoAODQuantity("HLT_DoublePFJets200_PNetBTag_0p11")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePFJets350_PNetBTag_0p11 = NanoAODQuantity("HLT_DoublePFJets350_PNetBTag_0p11")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePFJets40_PNetBTag_0p11 = NanoAODQuantity("HLT_DoublePFJets40_PNetBTag_0p11")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePhoton33_CaloIdL = NanoAODQuantity("HLT_DoublePhoton33_CaloIdL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePhoton70 = NanoAODQuantity("HLT_DoublePhoton70")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePhoton85 = NanoAODQuantity("HLT_DoublePhoton85")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ECALHT800 = NanoAODQuantity("HLT_ECALHT800")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_EcalCalibration = NanoAODQuantity("HLT_EcalCalibration")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele115_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele115_CaloIdVT_GsfTrkIdT")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity(
    "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele135_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele135_CaloIdVT_GsfTrkIdT")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53 = NanoAODQuantity(
    "HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL = NanoAODQuantity(
    "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele17_CaloIdM_TrackIdM_PFJet30 = NanoAODQuantity(
    "HLT_Ele17_CaloIdM_TrackIdM_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele35_WPTight_Gsf = NanoAODQuantity("HLT_Ele35_WPTight_Gsf")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele38_WPTight_Gsf = NanoAODQuantity("HLT_Ele38_WPTight_Gsf")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele40_WPTight_Gsf = NanoAODQuantity("HLT_Ele40_WPTight_Gsf")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_EphemeralPhysics = NanoAODQuantity("HLT_EphemeralPhysics")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_EphemeralZeroBias = NanoAODQuantity("HLT_EphemeralZeroBias")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack = NanoAODQuantity(
    "HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack = NanoAODQuantity(
    "HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT300_Beamspot = NanoAODQuantity("HLT_HT300_Beamspot")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive = NanoAODQuantity(
    "HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT350 = NanoAODQuantity("HLT_HT350")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT400_DisplacedDijet40_DisplacedTrack = NanoAODQuantity(
    "HLT_HT400_DisplacedDijet40_DisplacedTrack"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive = NanoAODQuantity(
    "HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT425 = NanoAODQuantity("HLT_HT425")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT550_DisplacedDijet60_Inclusive = NanoAODQuantity(
    "HLT_HT550_DisplacedDijet60_Inclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT650_DisplacedDijet60_Inclusive = NanoAODQuantity(
    "HLT_HT650_DisplacedDijet60_Inclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HcalCalibration = NanoAODQuantity("HLT_HcalCalibration")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HcalNZS = NanoAODQuantity("HLT_HcalNZS")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HcalPhiSym = NanoAODQuantity("HLT_HcalPhiSym")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HighPtTkMu100 = NanoAODQuantity("HLT_HighPtTkMu100")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu20 = NanoAODQuantity("HLT_IsoMu20")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu27 = NanoAODQuantity("HLT_IsoMu27")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1 = NanoAODQuantity(
    "HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoTrackHB = NanoAODQuantity("HLT_IsoTrackHB")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoTrackHE = NanoAODQuantity("HLT_IsoTrackHE")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoTrk200_L1SingleMuShower = NanoAODQuantity("HLT_IsoTrk200_L1SingleMuShower")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoTrk400_L1SingleMuShower = NanoAODQuantity("HLT_IsoTrk400_L1SingleMuShower")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1AXOVTight = NanoAODQuantity("HLT_L1AXOVTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1ETMHadSeeds = NanoAODQuantity("HLT_L1ETMHadSeeds")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Mu6HT240 = NanoAODQuantity("HLT_L1Mu6HT240")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1SingleLLPJet = NanoAODQuantity("HLT_L1SingleLLPJet")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1SingleMuCosmics = NanoAODQuantity("HLT_L1SingleMuCosmics")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 = NanoAODQuantity(
    "HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX = NanoAODQuantity(
    "HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX = NanoAODQuantity(
    "HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L3Mu10NoVtx = NanoAODQuantity("HLT_L3Mu10NoVtx")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L3Mu10NoVtx_DxyMin0p01cm = NanoAODQuantity("HLT_L3Mu10NoVtx_DxyMin0p01cm")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L3Mu30NoVtx_DxyMin0p01cm = NanoAODQuantity("HLT_L3Mu30NoVtx_DxyMin0p01cm")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L3Mu50NoVtx_DxyMin0p01cm = NanoAODQuantity("HLT_L3Mu50NoVtx_DxyMin0p01cm")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L3dTksMu10_NoVtx_DxyMin0p01cm = NanoAODQuantity("HLT_L3dTksMu10_NoVtx_DxyMin0p01cm")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1 = NanoAODQuantity(
    "HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_MET105_IsoTrk50 = NanoAODQuantity("HLT_MET105_IsoTrk50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_MET120_IsoTrk50 = NanoAODQuantity("HLT_MET120_IsoTrk50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu10_Barrel_L1HP11_IP6 = NanoAODQuantity("HLT_Mu10_Barrel_L1HP11_IP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12eta2p3 = NanoAODQuantity("HLT_Mu12eta2p3")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12eta2p3_PFJet40 = NanoAODQuantity("HLT_Mu12eta2p3_PFJet40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu15 = NanoAODQuantity("HLT_Mu15")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu18_Mu9_SameSign = NanoAODQuantity("HLT_Mu18_Mu9_SameSign")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu20 = NanoAODQuantity("HLT_Mu20")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId = NanoAODQuantity(
    "HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu25_TkMu0_Phi = NanoAODQuantity("HLT_Mu25_TkMu0_Phi")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu27 = NanoAODQuantity("HLT_Mu27")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu27_Ele37_CaloIdL_MW = NanoAODQuantity("HLT_Mu27_Ele37_CaloIdL_MW")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL = NanoAODQuantity(
    "HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL = NanoAODQuantity(
    "HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL = NanoAODQuantity(
    "HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL = NanoAODQuantity(
    "HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu4_L1DoubleMu = NanoAODQuantity("HLT_Mu4_L1DoubleMu")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu55 = NanoAODQuantity("HLT_Mu55")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu6_Barrel_L1HP7_IP6 = NanoAODQuantity("HLT_Mu6_Barrel_L1HP7_IP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu7_Barrel_L1HP8_IP6 = NanoAODQuantity("HLT_Mu7_Barrel_L1HP8_IP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu9_Barrel_L1HP10_IP6 = NanoAODQuantity("HLT_Mu9_Barrel_L1HP10_IP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT1050 = NanoAODQuantity("HLT_PFHT1050")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT180 = NanoAODQuantity("HLT_PFHT180")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70 = NanoAODQuantity(
    "HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT350 = NanoAODQuantity("HLT_PFHT350")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT370 = NanoAODQuantity("HLT_PFHT370")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT430 = NanoAODQuantity("HLT_PFHT430")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT510 = NanoAODQuantity("HLT_PFHT510")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT590 = NanoAODQuantity("HLT_PFHT590")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT680 = NanoAODQuantity("HLT_PFHT680")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT700_PFMET85_PFMHT85_IDTight = NanoAODQuantity(
    "HLT_PFHT700_PFMET85_PFMHT85_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT780 = NanoAODQuantity("HLT_PFHT780")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT800_PFMET75_PFMHT75_IDTight = NanoAODQuantity(
    "HLT_PFHT800_PFMET75_PFMHT75_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT890 = NanoAODQuantity("HLT_PFHT890")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet110 = NanoAODQuantity("HLT_PFJet110")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet140 = NanoAODQuantity("HLT_PFJet140")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet260 = NanoAODQuantity("HLT_PFJet260")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet320 = NanoAODQuantity("HLT_PFJet320")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet40 = NanoAODQuantity("HLT_PFJet40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet400 = NanoAODQuantity("HLT_PFJet400")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet40_GPUvsCPU = NanoAODQuantity("HLT_PFJet40_GPUvsCPU")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet450 = NanoAODQuantity("HLT_PFJet450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet500 = NanoAODQuantity("HLT_PFJet500")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet550 = NanoAODQuantity("HLT_PFJet550")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet60 = NanoAODQuantity("HLT_PFJet60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet80 = NanoAODQuantity("HLT_PFJet80")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd140 = NanoAODQuantity("HLT_PFJetFwd140")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd200 = NanoAODQuantity("HLT_PFJetFwd200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd260 = NanoAODQuantity("HLT_PFJetFwd260")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd320 = NanoAODQuantity("HLT_PFJetFwd320")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd40 = NanoAODQuantity("HLT_PFJetFwd40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd400 = NanoAODQuantity("HLT_PFJetFwd400")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd450 = NanoAODQuantity("HLT_PFJetFwd450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd500 = NanoAODQuantity("HLT_PFJetFwd500")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd60 = NanoAODQuantity("HLT_PFJetFwd60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJetFwd80 = NanoAODQuantity("HLT_PFJetFwd80")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET105_IsoTrk50 = NanoAODQuantity("HLT_PFMET105_IsoTrk50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET130_PFMHT130_IDTight = NanoAODQuantity("HLT_PFMET130_PFMHT130_IDTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET140_PFMHT140_IDTight = NanoAODQuantity("HLT_PFMET140_PFMHT140_IDTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET250_NotCleaned = NanoAODQuantity("HLT_PFMET250_NotCleaned")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET300_NotCleaned = NanoAODQuantity("HLT_PFMET300_NotCleaned")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF = NanoAODQuantity(
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETTypeOne140_PFMHT140_IDTight = NanoAODQuantity(
    "HLT_PFMETTypeOne140_PFMHT140_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETTypeOne200_BeamHaloCleaned = NanoAODQuantity(
    "HLT_PFMETTypeOne200_BeamHaloCleaned"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PPSMaxTracksPerArm1 = NanoAODQuantity("HLT_PPSMaxTracksPerArm1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PPSMaxTracksPerRP4 = NanoAODQuantity("HLT_PPSMaxTracksPerRP4")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PPSRandom = NanoAODQuantity("HLT_PPSRandom")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon100EBHE10 = NanoAODQuantity("HLT_Photon100EBHE10")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon120 = NanoAODQuantity("HLT_Photon120")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon120_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon120_R9Id90_HE10_IsoM")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon150 = NanoAODQuantity("HLT_Photon150")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon165_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon165_R9Id90_HE10_IsoM")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon175 = NanoAODQuantity("HLT_Photon175")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon200 = NanoAODQuantity("HLT_Photon200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon20_HoverELoose = NanoAODQuantity("HLT_Photon20_HoverELoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon300_NoHE = NanoAODQuantity("HLT_Photon300_NoHE")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon30EB_TightID_TightIso = NanoAODQuantity("HLT_Photon30EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon30_HoverELoose = NanoAODQuantity("HLT_Photon30_HoverELoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon32_OneProng32_M50To105 = NanoAODQuantity("HLT_Photon32_OneProng32_M50To105")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon33 = NanoAODQuantity("HLT_Photon33")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34 = NanoAODQuantity(
    "HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon35_TwoProngs35 = NanoAODQuantity("HLT_Photon35_TwoProngs35")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon40EB = NanoAODQuantity("HLT_Photon40EB")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon40EB_TightID_TightIso = NanoAODQuantity("HLT_Photon40EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon45EB = NanoAODQuantity("HLT_Photon45EB")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon45EB_TightID_TightIso = NanoAODQuantity("HLT_Photon45EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50EB = NanoAODQuantity("HLT_Photon50EB")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon55EB_TightID_TightIso = NanoAODQuantity("HLT_Photon55EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon75 = NanoAODQuantity("HLT_Photon75")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon75EB_TightID_TightIso = NanoAODQuantity("HLT_Photon75EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon90 = NanoAODQuantity("HLT_Photon90")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon90EB_TightID_TightIso = NanoAODQuantity("HLT_Photon90EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon90_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon90_R9Id90_HE10_IsoM")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Physics = NanoAODQuantity("HLT_Physics")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Random = NanoAODQuantity("HLT_Random")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_SingleEle8 = NanoAODQuantity("HLT_SingleEle8")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_SingleEle8_SingleEGL1 = NanoAODQuantity("HLT_SingleEle8_SingleEGL1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Trimuon5_3p5_2_Upsilon_Muon = NanoAODQuantity("HLT_Trimuon5_3p5_2_Upsilon_Muon")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon = NanoAODQuantity(
    "HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx = NanoAODQuantity(
    "HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_UncorrectedJetE60_NoBPTX3BX = NanoAODQuantity("HLT_UncorrectedJetE60_NoBPTX3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_UncorrectedJetE70_NoBPTX3BX = NanoAODQuantity("HLT_UncorrectedJetE70_NoBPTX3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8DiPFJet250_250_SoftDropMass40 = NanoAODQuantity(
    "HLT_AK8DiPFJet250_250_SoftDropMass40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8DiPFJet250_250_SoftDropMass50 = NanoAODQuantity(
    "HLT_AK8DiPFJet250_250_SoftDropMass50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8DiPFJet260_260_SoftDropMass30 = NanoAODQuantity(
    "HLT_AK8DiPFJet260_260_SoftDropMass30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8DiPFJet260_260_SoftDropMass40 = NanoAODQuantity(
    "HLT_AK8DiPFJet260_260_SoftDropMass40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8PFJet220_SoftDropMass40 = NanoAODQuantity("HLT_AK8PFJet220_SoftDropMass40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50 = NanoAODQuantity(
    "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53 = NanoAODQuantity(
    "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55 = NanoAODQuantity(
    "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60 = NanoAODQuantity(
    "HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8PFJet230_SoftDropMass40 = NanoAODQuantity("HLT_AK8PFJet230_SoftDropMass40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10 = NanoAODQuantity(
    "HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03 = NanoAODQuantity(
    "HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05 = NanoAODQuantity(
    "HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10 = NanoAODQuantity(
    "HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03 = NanoAODQuantity(
    "HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05 = NanoAODQuantity(
    "HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8PFJet275_Nch40 = NanoAODQuantity("HLT_AK8PFJet275_Nch40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet275_Nch45 = NanoAODQuantity("HLT_AK8PFJet275_Nch45")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10 = NanoAODQuantity(
    "HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03 = NanoAODQuantity(
    "HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05 = NanoAODQuantity(
    "HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_BTagMu_AK4DiJet110_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet110_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK4DiJet170_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet170_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK4DiJet20_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet20_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK4DiJet40_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet40_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK4DiJet70_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet70_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK4Jet300_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4Jet300_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK8DiJet170_Mu5 = NanoAODQuantity("HLT_BTagMu_AK8DiJet170_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK8Jet170_DoubleMu5 = NanoAODQuantity("HLT_BTagMu_AK8Jet170_DoubleMu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_BTagMu_AK8Jet300_Mu5 = NanoAODQuantity("HLT_BTagMu_AK8Jet300_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_CDC_L2cosmic_10_er1p0 = NanoAODQuantity("HLT_CDC_L2cosmic_10_er1p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CDC_L2cosmic_5p5_er1p0 = NanoAODQuantity("HLT_CDC_L2cosmic_5p5_er1p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_CaloMET60_DTCluster50 = NanoAODQuantity("HLT_CaloMET60_DTCluster50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CaloMET60_DTClusterNoMB1S50 = NanoAODQuantity("HLT_CaloMET60_DTClusterNoMB1S50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_CscCluster_Loose = NanoAODQuantity("HLT_CscCluster_Loose")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CscCluster_Medium = NanoAODQuantity("HLT_CscCluster_Medium")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CscCluster_Tight = NanoAODQuantity("HLT_CscCluster_Tight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_CscCluster100_Ele5 = NanoAODQuantity("HLT_CscCluster100_Ele5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CscCluster100_Mu5 = NanoAODQuantity("HLT_CscCluster100_Mu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CscCluster100_PNetTauhPFJet10_Loose = NanoAODQuantity(
    "HLT_CscCluster100_PNetTauhPFJet10_Loose"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_CscCluster50_Photon20Unseeded = NanoAODQuantity("HLT_CscCluster50_Photon20Unseeded")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_CscCluster50_Photon30Unseeded = NanoAODQuantity("HLT_CscCluster50_Photon30Unseeded")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DiMu9_Ele9_CaloIdL_TrackIdL = NanoAODQuantity("HLT_DiMu9_Ele9_CaloIdL_TrackIdL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ = NanoAODQuantity(
    "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_Jpsi3p5_Muon2 = NanoAODQuantity("HLT_Dimuon0_Jpsi3p5_Muon2")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_LowMass = NanoAODQuantity("HLT_Dimuon0_LowMass")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_Jpsi = NanoAODQuantity("HLT_Dimuon0_Jpsi")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_Jpsi_L1_4R_0er1p5R = NanoAODQuantity("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_Jpsi_L1_NoOS = NanoAODQuantity("HLT_Dimuon0_Jpsi_L1_NoOS")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_Jpsi_NoVertexing = NanoAODQuantity("HLT_Dimuon0_Jpsi_NoVertexing")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R = NanoAODQuantity(
    "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_Jpsi_NoVertexing_NoOS = NanoAODQuantity("HLT_Dimuon0_Jpsi_NoVertexing_NoOS")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_LowMass_L1_0er1p5 = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_0er1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_LowMass_L1_4 = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_4")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_LowMass_L1_TM530 = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_TM530")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_Upsilon_Muon_NoL1Mass = NanoAODQuantity("HLT_Dimuon0_Upsilon_Muon_NoL1Mass")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_Upsilon_NoVertexing = NanoAODQuantity("HLT_Dimuon0_Upsilon_NoVertexing")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon0_Upsilon_L1_4p5 = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_Upsilon_L1_4p5er2p0 = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon0_Upsilon_L1_4p5er2p0M = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5er2p0M")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon14_Phi_Barrel_Seagulls = NanoAODQuantity("HLT_Dimuon14_Phi_Barrel_Seagulls")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon14_PsiPrime = NanoAODQuantity("HLT_Dimuon14_PsiPrime")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon14_PsiPrime_noCorrL1 = NanoAODQuantity("HLT_Dimuon14_PsiPrime_noCorrL1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon18_PsiPrime = NanoAODQuantity("HLT_Dimuon18_PsiPrime")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon18_PsiPrime_noCorrL1 = NanoAODQuantity("HLT_Dimuon18_PsiPrime_noCorrL1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon24_Phi_noCorrL1 = NanoAODQuantity("HLT_Dimuon24_Phi_noCorrL1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon24_Upsilon_noCorrL1 = NanoAODQuantity("HLT_Dimuon24_Upsilon_noCorrL1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Dimuon25_Jpsi = NanoAODQuantity("HLT_Dimuon25_Jpsi")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Dimuon25_Jpsi_noCorrL1 = NanoAODQuantity("HLT_Dimuon25_Jpsi_noCorrL1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT = NanoAODQuantity(
    "HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT = NanoAODQuantity(
    "HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT = NanoAODQuantity(
    "HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT = NanoAODQuantity(
    "HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT = NanoAODQuantity(
    "HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT = NanoAODQuantity(
    "HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId = NanoAODQuantity(
    "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55 = NanoAODQuantity(
    "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 = NanoAODQuantity(
    "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95 = NanoAODQuantity(
    "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DiphotonMVA14p25_Mass90 = NanoAODQuantity("HLT_DiphotonMVA14p25_Mass90")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DiphotonMVA14p25_Tight_Mass90 = NanoAODQuantity("HLT_DiphotonMVA14p25_Tight_Mass90")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleEle8_eta1p22_mMax6 = NanoAODQuantity("HLT_DoubleEle8_eta1p22_mMax6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350 = NanoAODQuantity(
    "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350 = NanoAODQuantity(
    "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleL2Mu23NoVtx_2Cha = NanoAODQuantity("HLT_DoubleL2Mu23NoVtx_2Cha")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed = NanoAODQuantity(
    "HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleL2Mu25NoVtx_2Cha = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed = NanoAODQuantity(
    "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4 = NanoAODQuantity(
    "HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4 = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4 = NanoAODQuantity(
    "HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4 = NanoAODQuantity("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1 = NanoAODQuantity(
    "HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy = NanoAODQuantity(
    "HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng = NanoAODQuantity(
    "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60 = NanoAODQuantity(
    "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75 = NanoAODQuantity(
    "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 = NanoAODQuantity(
    "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu2_Jpsi_LowPt = NanoAODQuantity("HLT_DoubleMu2_Jpsi_LowPt")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon = NanoAODQuantity(
    "HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu3_TkMu_DsTau3Mu = NanoAODQuantity("HLT_DoubleMu3_TkMu_DsTau3Mu")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu3_DCA_PFMET50_PFMHT60 = NanoAODQuantity("HLT_DoubleMu3_DCA_PFMET50_PFMHT60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0 = NanoAODQuantity(
    "HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA = NanoAODQuantity(
    "HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu3_DZ_PFMET50_PFMHT60 = NanoAODQuantity("HLT_DoubleMu3_DZ_PFMET50_PFMHT60")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu3_DZ_PFMET70_PFMHT70 = NanoAODQuantity("HLT_DoubleMu3_DZ_PFMET70_PFMHT70")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu3_DZ_PFMET90_PFMHT90 = NanoAODQuantity("HLT_DoubleMu3_DZ_PFMET90_PFMHT90")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu3_Trk_Tau3mu = NanoAODQuantity("HLT_DoubleMu3_Trk_Tau3mu")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass = NanoAODQuantity("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu4_JpsiTrkTrk_Displaced = NanoAODQuantity(
    "HLT_DoubleMu4_JpsiTrkTrk_Displaced"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_JpsiTrk_Bc = NanoAODQuantity("HLT_DoubleMu4_JpsiTrk_Bc")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_LowMass_Displaced = NanoAODQuantity("HLT_DoubleMu4_LowMass_Displaced")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_Mass3p8_DZ_PFHT350 = NanoAODQuantity("HLT_DoubleMu4_Mass3p8_DZ_PFHT350")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_MuMuTrk_Displaced = NanoAODQuantity("HLT_DoubleMu4_MuMuTrk_Displaced")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu4_3_Bs = NanoAODQuantity("HLT_DoubleMu4_3_Bs")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG = NanoAODQuantity(
    "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_3_Jpsi = NanoAODQuantity("HLT_DoubleMu4_3_Jpsi")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_3_LowMass = NanoAODQuantity("HLT_DoubleMu4_3_LowMass")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_3_LowMass_SS = NanoAODQuantity("HLT_DoubleMu4_3_LowMass_SS")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_3_Photon4_BsToMMG = NanoAODQuantity("HLT_DoubleMu4_3_Photon4_BsToMMG")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoubleMu4_Jpsi_Displaced = NanoAODQuantity("HLT_DoubleMu4_Jpsi_Displaced")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoubleMu4_Jpsi_NoVertexing = NanoAODQuantity("HLT_DoubleMu4_Jpsi_NoVertexing")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60 = NanoAODQuantity(
    "HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75 = NanoAODQuantity(
    "HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele15_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele15_IsoVVVL_PFHT450_PFMET50 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT450_PFMET50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele15_IsoVVVL_PFHT600 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT600")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity(
    "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele23_CaloIdM_TrackIdM_PFJet30 = NanoAODQuantity(
    "HLT_Ele23_CaloIdM_TrackIdM_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity(
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity(
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 = NanoAODQuantity(
    "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele28_HighEta_SC20_Mass55 = NanoAODQuantity("HLT_Ele28_HighEta_SC20_Mass55")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = NanoAODQuantity(
    "HLT_Ele28_eta2p1_WPTight_Gsf_HT150"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele30_WPTight_Gsf = NanoAODQuantity("HLT_Ele30_WPTight_Gsf")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = NanoAODQuantity(
    "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele32_WPTight_Gsf = NanoAODQuantity("HLT_Ele32_WPTight_Gsf")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele32_WPTight_Gsf_L1DoubleEG = NanoAODQuantity("HLT_Ele32_WPTight_Gsf_L1DoubleEG")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele50_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Ele50_IsoVVVL_PFHT450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = NanoAODQuantity(
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40 = NanoAODQuantity(
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40 = NanoAODQuantity(
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10 = NanoAODQuantity(
    "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity(
    "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Ele8_CaloIdM_TrackIdM_PFJet30 = NanoAODQuantity("HLT_Ele8_CaloIdM_TrackIdM_PFJet30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7 = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8 = NanoAODQuantity(
    "HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive = NanoAODQuantity(
    "HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive = NanoAODQuantity(
    "HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive = NanoAODQuantity(
    "HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT350_DelayedJet40_SingleDelay3nsInclusive = NanoAODQuantity(
    "HLT_HT350_DelayedJet40_SingleDelay3nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive = NanoAODQuantity(
    "HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive = NanoAODQuantity(
    "HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay2nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay2nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive = NanoAODQuantity(
    "HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_HT430_DisplacedDijet40_DisplacedTrack = NanoAODQuantity(
    "HLT_HT430_DisplacedDijet40_DisplacedTrack"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5 = NanoAODQuantity(
    "HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24 = NanoAODQuantity("HLT_IsoMu24")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_OneProng32 = NanoAODQuantity("HLT_IsoMu24_OneProng32")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_TwoProngs35 = NanoAODQuantity("HLT_IsoMu24_TwoProngs35")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1 = NanoAODQuantity("HLT_IsoMu24_eta2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1 = (
    NanoAODQuantity(
        "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1 = (
    NanoAODQuantity(
        "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1 = (
    NanoAODQuantity(
        "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_PFHT250 = NanoAODQuantity("HLT_IsoMu24_eta2p1_PFHT250")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring = (
    NanoAODQuantity(
        "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring = (
    NanoAODQuantity(
        "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1 = NanoAODQuantity(
    "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring = (
    NanoAODQuantity(
        "HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu50_AK8PFJet220_SoftDropMass40 = NanoAODQuantity(
    "HLT_IsoMu50_AK8PFJet220_SoftDropMass40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_IsoMu50_AK8PFJet230_SoftDropMass40 = NanoAODQuantity(
    "HLT_IsoMu50_AK8PFJet230_SoftDropMass40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06 = NanoAODQuantity(
    "HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10 = NanoAODQuantity(
    "HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L1CSCShower_DTCluster50 = NanoAODQuantity("HLT_L1CSCShower_DTCluster50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1CSCShower_DTCluster75 = NanoAODQuantity("HLT_L1CSCShower_DTCluster75")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L1MET_DTCluster50 = NanoAODQuantity("HLT_L1MET_DTCluster50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1MET_DTClusterNoMB1S50 = NanoAODQuantity("HLT_L1MET_DTClusterNoMB1S50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive = NanoAODQuantity(
    "HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L2Mu10NoVtx_2Cha = NanoAODQuantity("HLT_L2Mu10NoVtx_2Cha")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu10NoVtx_2Cha_CosmicSeed = NanoAODQuantity("HLT_L2Mu10NoVtx_2Cha_CosmicSeed")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L2Mu10_NoVertex_NoBPTX = NanoAODQuantity("HLT_L2Mu10_NoVertex_NoBPTX")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu10_NoVertex_NoBPTX3BX = NanoAODQuantity("HLT_L2Mu10_NoVertex_NoBPTX3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L2Mu23NoVtx_2Cha = NanoAODQuantity("HLT_L2Mu23NoVtx_2Cha")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu23NoVtx_2Cha_CosmicSeed = NanoAODQuantity("HLT_L2Mu23NoVtx_2Cha_CosmicSeed")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm = NanoAODQuantity(
    "HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu0_L1DoubleMu = NanoAODQuantity("HLT_Mu0_L1DoubleMu")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu0_Barrel = NanoAODQuantity("HLT_Mu0_Barrel")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP10 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP10")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP11 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP11")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP6 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP6_IP6 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP6_IP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP7 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP7")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP8 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP8")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu0_Barrel_L1HP9 = NanoAODQuantity("HLT_Mu0_Barrel_L1HP9")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu12_DoublePFJets100_PNetBTag_0p11 = NanoAODQuantity(
    "HLT_Mu12_DoublePFJets100_PNetBTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_DoublePFJets200_PNetBTag_0p11 = NanoAODQuantity(
    "HLT_Mu12_DoublePFJets200_PNetBTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_DoublePFJets350_PNetBTag_0p11 = NanoAODQuantity(
    "HLT_Mu12_DoublePFJets350_PNetBTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11 = NanoAODQuantity(
    "HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_DoublePFJets40_PNetBTag_0p11 = NanoAODQuantity(
    "HLT_Mu12_DoublePFJets40_PNetBTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11 = NanoAODQuantity(
    "HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53 = NanoAODQuantity(
    "HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity(
    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity(
    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu15_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu15_IsoVVVL_PFHT450_PFMET50 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT450_PFMET50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu15_IsoVVVL_PFHT600 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT600")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu17 = NanoAODQuantity("HLT_Mu17")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_Photon30_IsoCaloId = NanoAODQuantity("HLT_Mu17_Photon30_IsoCaloId")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_TrkIsoVVL = NanoAODQuantity("HLT_Mu17_TrkIsoVVL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = NanoAODQuantity("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30 = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30 = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30 = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30 = NanoAODQuantity(
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu19 = NanoAODQuantity("HLT_Mu19")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu19_TrkIsoVVL = NanoAODQuantity("HLT_Mu19_TrkIsoVVL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL = NanoAODQuantity("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ = NanoAODQuantity(
    "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 = NanoAODQuantity(
    "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 = NanoAODQuantity(
    "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity(
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity(
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu3_L1SingleMu5orSingleMu7 = NanoAODQuantity("HLT_Mu3_L1SingleMu5orSingleMu7")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu3_PFJet40 = NanoAODQuantity("HLT_Mu3_PFJet40")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu30_TkMu0_Psi = NanoAODQuantity("HLT_Mu30_TkMu0_Psi")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu30_TkMu0_Upsilon = NanoAODQuantity("HLT_Mu30_TkMu0_Upsilon")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu37_Ele27_CaloIdL_MW = NanoAODQuantity("HLT_Mu37_Ele27_CaloIdL_MW")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu37_TkMu27 = NanoAODQuantity("HLT_Mu37_TkMu27")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight = NanoAODQuantity(
    "HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight = NanoAODQuantity(
    "HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight = NanoAODQuantity(
    "HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight = NanoAODQuantity(
    "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight = NanoAODQuantity(
    "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight = NanoAODQuantity(
    "HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu50 = NanoAODQuantity("HLT_Mu50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu50_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Mu50_IsoVVVL_PFHT450")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu50_L1SingleMuShower = NanoAODQuantity("HLT_Mu50_L1SingleMuShower")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5 = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5 = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5 = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5 = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose = NanoAODQuantity(
    "HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu7p5_L2Mu2_Jpsi = NanoAODQuantity("HLT_Mu7p5_L2Mu2_Jpsi")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu7p5_L2Mu2_Upsilon = NanoAODQuantity("HLT_Mu7p5_L2Mu2_Upsilon")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8 = NanoAODQuantity("HLT_Mu8")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_Barrel_L1HP9_IP6 = NanoAODQuantity("HLT_Mu8_Barrel_L1HP9_IP6")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL = NanoAODQuantity("HLT_Mu8_TrkIsoVVL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_DiEle12_CaloIdL_TrackIdL = NanoAODQuantity("HLT_Mu8_DiEle12_CaloIdL_TrackIdL")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ = NanoAODQuantity(
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350 = NanoAODQuantity(
    "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ = NanoAODQuantity(
    "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5 = (
    NanoAODQuantity(
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50 = (
    NanoAODQuantity(
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50"
    )
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55 = NanoAODQuantity(
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT250 = NanoAODQuantity("HLT_PFHT250")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT250_QuadPFJet25 = NanoAODQuantity("HLT_PFHT250_QuadPFJet25")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50 = NanoAODQuantity(
    "HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55 = NanoAODQuantity(
    "HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50 = NanoAODQuantity(
    "HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55 = NanoAODQuantity(
    "HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60 = NanoAODQuantity(
    "HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT280_QuadPFJet30 = NanoAODQuantity("HLT_PFHT280_QuadPFJet30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50 = NanoAODQuantity(
    "HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55 = NanoAODQuantity(
    "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60 = NanoAODQuantity(
    "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT330PT30_QuadPFJet_75_60_45_40 = NanoAODQuantity(
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5 = NanoAODQuantity(
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0 = NanoAODQuantity(
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3 = NanoAODQuantity(
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT400_SixPFJet32 = NanoAODQuantity("HLT_PFHT400_SixPFJet32")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50 = NanoAODQuantity(
    "HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT400_FivePFJet_120_120_60_30_30 = NanoAODQuantity(
    "HLT_PFHT400_FivePFJet_120_120_60_30_30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3 = NanoAODQuantity(
    "HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6 = NanoAODQuantity(
    "HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT450_SixPFJet36 = NanoAODQuantity("HLT_PFHT450_SixPFJet36")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT450_SixPFJet36_PNetBTag0p35 = NanoAODQuantity(
    "HLT_PFHT450_SixPFJet36_PNetBTag0p35"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFHT500_PFMET100_PFMHT100_IDTight = NanoAODQuantity(
    "HLT_PFHT500_PFMET100_PFMHT100_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFHT500_PFMET110_PFMHT110_IDTight = NanoAODQuantity(
    "HLT_PFHT500_PFMET110_PFMHT110_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFJet200 = NanoAODQuantity("HLT_PFJet200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet200_TimeGt2p5ns = NanoAODQuantity("HLT_PFJet200_TimeGt2p5ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFJet200_TimeLtNeg2p5ns = NanoAODQuantity("HLT_PFJet200_TimeLtNeg2p5ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFMET120_PFMHT120_IDTight = NanoAODQuantity("HLT_PFMET120_PFMHT120_IDTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET120_PFMHT120_IDTight_PFHT60 = NanoAODQuantity(
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFMET200_BeamHaloCleaned = NanoAODQuantity("HLT_PFMET200_BeamHaloCleaned")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMET200_NotCleaned = NanoAODQuantity("HLT_PFMET200_NotCleaned")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = NanoAODQuantity(
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF = NanoAODQuantity(
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = NanoAODQuantity(
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = NanoAODQuantity(
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF = NanoAODQuantity(
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = NanoAODQuantity(
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF = NanoAODQuantity(
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Photon110EB_TightID_TightIso = NanoAODQuantity("HLT_Photon110EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon110EB_TightID_TightIso_AK8CaloJet30 = NanoAODQuantity(
    "HLT_Photon110EB_TightID_TightIso_AK8CaloJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon110EB_TightID_TightIso_AK8PFJet30 = NanoAODQuantity(
    "HLT_Photon110EB_TightID_TightIso_AK8PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon110EB_TightID_TightIso_CaloJet30 = NanoAODQuantity(
    "HLT_Photon110EB_TightID_TightIso_CaloJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon110EB_TightID_TightIso_PFJet30 = NanoAODQuantity(
    "HLT_Photon110EB_TightID_TightIso_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Photon50 = NanoAODQuantity("HLT_Photon50")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon50_R9Id90_HE10_IsoM")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50_TimeGt2p5ns = NanoAODQuantity("HLT_Photon50_TimeGt2p5ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50_TimeLtNeg2p5ns = NanoAODQuantity("HLT_Photon50_TimeLtNeg2p5ns")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Photon50EB_TightID_TightIso = NanoAODQuantity("HLT_Photon50EB_TightID_TightIso")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50EB_TightID_TightIso_AK8CaloJet30 = NanoAODQuantity(
    "HLT_Photon50EB_TightID_TightIso_AK8CaloJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50EB_TightID_TightIso_AK8PFJet30 = NanoAODQuantity(
    "HLT_Photon50EB_TightID_TightIso_AK8PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50EB_TightID_TightIso_CaloJet30 = NanoAODQuantity(
    "HLT_Photon50EB_TightID_TightIso_CaloJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon50EB_TightID_TightIso_PFJet30 = NanoAODQuantity(
    "HLT_Photon50EB_TightID_TightIso_PFJet30"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3 = NanoAODQuantity(
    "HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350 = NanoAODQuantity(
    "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380 = NanoAODQuantity(
    "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400 = NanoAODQuantity(
    "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Photon75_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon75_R9Id90_HE10_IsoM")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3 = NanoAODQuantity(
    "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_QuadPFJet100_88_70_30 = NanoAODQuantity("HLT_QuadPFJet100_88_70_30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight = NanoAODQuantity(
    "HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_QuadPFJet103_88_75_15 = NanoAODQuantity("HLT_QuadPFJet103_88_75_15")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1 = NanoAODQuantity(
    "HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2 = NanoAODQuantity(
    "HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_QuadPFJet105_88_75_30 = NanoAODQuantity("HLT_QuadPFJet105_88_75_30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight = NanoAODQuantity(
    "HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_QuadPFJet105_88_76_15 = NanoAODQuantity("HLT_QuadPFJet105_88_76_15")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1 = NanoAODQuantity(
    "HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2 = NanoAODQuantity(
    "HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_QuadPFJet111_90_80_30 = NanoAODQuantity("HLT_QuadPFJet111_90_80_30")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight = NanoAODQuantity(
    "HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_QuadPFJet111_90_80_15 = NanoAODQuantity("HLT_QuadPFJet111_90_80_15")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1 = NanoAODQuantity(
    "HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2 = NanoAODQuantity(
    "HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = NanoAODQuantity("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = NanoAODQuantity(
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15 = NanoAODQuantity("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1 = NanoAODQuantity(
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_TripleMu_10_5_5_DZ = NanoAODQuantity("HLT_TripleMu_10_5_5_DZ")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_TripleMu_12_10_5 = NanoAODQuantity("HLT_TripleMu_12_10_5")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_TripleMu_5_3_3_Mass3p8_DCA = NanoAODQuantity("HLT_TripleMu_5_3_3_Mass3p8_DCA")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_TripleMu_5_3_3_Mass3p8_DZ = NanoAODQuantity("HLT_TripleMu_5_3_3_Mass3p8_DZ")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_UncorrectedJetE30_NoBPTX = NanoAODQuantity("HLT_UncorrectedJetE30_NoBPTX")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_UncorrectedJetE30_NoBPTX3BX = NanoAODQuantity("HLT_UncorrectedJetE30_NoBPTX3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1 = NanoAODQuantity(
    "HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DoublePNetTauhPFJet20_eta2p2 = NanoAODQuantity(
    "HLT_VBF_DoublePNetTauhPFJet20_eta2p2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet125_45_Mjj1050 = NanoAODQuantity("HLT_VBF_DiPFJet125_45_Mjj1050")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet125_45_Mjj1200 = NanoAODQuantity("HLT_VBF_DiPFJet125_45_Mjj1200")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1 = NanoAODQuantity(
    "HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1 = NanoAODQuantity(
    "HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3 = NanoAODQuantity(
    "HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf = NanoAODQuantity(
    "HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet50_Mjj750_Photon22 = NanoAODQuantity("HLT_VBF_DiPFJet50_Mjj750_Photon22")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf = NanoAODQuantity(
    "HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet50_Mjj650_Photon22 = NanoAODQuantity("HLT_VBF_DiPFJet50_Mjj650_Photon22")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60 = NanoAODQuantity(
    "HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60 = NanoAODQuantity(
    "HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85 = NanoAODQuantity(
    "HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85 = NanoAODQuantity(
    "HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL = NanoAODQuantity(
    "HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL = NanoAODQuantity(
    "HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HLT_ZeroBias = NanoAODQuantity("HLT_ZeroBias")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_Alignment = NanoAODQuantity("HLT_ZeroBias_Alignment")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_Beamspot = NanoAODQuantity("HLT_ZeroBias_Beamspot")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_FirstBXAfterTrain = NanoAODQuantity("HLT_ZeroBias_FirstBXAfterTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_FirstCollisionAfterAbortGap = NanoAODQuantity(
    "HLT_ZeroBias_FirstCollisionAfterAbortGap"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_FirstCollisionInTrain = NanoAODQuantity(
    "HLT_ZeroBias_FirstCollisionInTrain"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_IsolatedBunches = NanoAODQuantity("HLT_ZeroBias_IsolatedBunches")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLT_ZeroBias_LastCollisionInTrain = NanoAODQuantity("HLT_ZeroBias_LastCollisionInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

HTXS_Mjj = NanoAODQuantity("HTXS_Mjj")
"""dtype: Float_t; description: invariant mass of the dijet (pt>30) system as identified in HTXS """
HTXS_V_pt = NanoAODQuantity("HTXS_V_pt")
"""dtype: Float_t; description: pt of the vector boson as identified in HTXS """
HTXS_dPhijj = NanoAODQuantity("HTXS_dPhijj")
"""dtype: Float_t; description: DeltaPhi between jets (pt>30) in dijet system as identified in HTXS """
HTXS_njets25 = NanoAODQuantity("HTXS_njets25")
"""dtype: UChar_t; description: number of jets with pt>25 GeV as identified in HTXS """
HTXS_njets30 = NanoAODQuantity("HTXS_njets30")
"""dtype: UChar_t; description: number of jets with pt>30 GeV as identified in HTXS """
HTXS_ptHjj = NanoAODQuantity("HTXS_ptHjj")
"""dtype: Float_t; description: pt of the dijet(pt>30)-plus-higgs system as identified in HTXS """

HTXS_Higgs_pt = NanoAODQuantity("HTXS_Higgs_pt")
"""dtype: Float_t; description: pt of the Higgs boson as identified in HTXS """
HTXS_Higgs_y = NanoAODQuantity("HTXS_Higgs_y")
"""dtype: Float_t; description: rapidity of the Higgs boson as identified in HTXS """

HTXS_stage_0 = NanoAODQuantity("HTXS_stage_0")
"""dtype: Int_t; description: HTXS stage-0 category """

HTXS_stage1_1_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_1_cat_pTjet25GeV")
"""dtype: Int_t; description: HTXS stage-1.1 category(jet pt>25 GeV) """
HTXS_stage1_1_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_1_cat_pTjet30GeV")
"""dtype: Int_t; description: HTXS stage-1.1 category(jet pt>30 GeV) """

HTXS_stage1_1_fine_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_1_fine_cat_pTjet25GeV")
"""dtype: Int_t; description: HTXS stage-1.1-fine category(jet pt>25 GeV) """
HTXS_stage1_1_fine_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_1_fine_cat_pTjet30GeV")
"""dtype: Int_t; description: HTXS stage-1.1-fine category(jet pt>30 GeV) """

HTXS_stage1_2_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_2_cat_pTjet25GeV")
"""dtype: Int_t; description: HTXS stage-1.2 category(jet pt>25 GeV) """
HTXS_stage1_2_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_2_cat_pTjet30GeV")
"""dtype: Int_t; description: HTXS stage-1.2 category(jet pt>30 GeV) """

HTXS_stage1_2_fine_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_2_fine_cat_pTjet25GeV")
"""dtype: Int_t; description: HTXS stage-1.2-fine category(jet pt>25 GeV) """
HTXS_stage1_2_fine_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_2_fine_cat_pTjet30GeV")
"""dtype: Int_t; description: HTXS stage-1.2-fine category(jet pt>30 GeV) """

HTXS_stage_1_pTjet25 = NanoAODQuantity("HTXS_stage_1_pTjet25")
"""dtype: Int_t; description: HTXS stage-1 category (jet pt>25 GeV) """
HTXS_stage_1_pTjet30 = NanoAODQuantity("HTXS_stage_1_pTjet30")
"""dtype: Int_t; description: HTXS stage-1 category (jet pt>30 GeV) """

nIsoTrack = NanoAODQuantity("nIsoTrack")
"""dtype: Int_t; description: isolated tracks after basic selection (((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && ((abs(dxy) < 0.2 && abs(dz) < 0.1) || pt>15) && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) and lepton veto """
IsoTrack_charge = NanoAODQuantity("IsoTrack_charge")
"""dtype: Short_t; description: electric charge """
IsoTrack_dxy = NanoAODQuantity("IsoTrack_dxy")
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm """
IsoTrack_dz = NanoAODQuantity("IsoTrack_dz")
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm """
IsoTrack_eta = NanoAODQuantity("IsoTrack_eta")
"""dtype: Float_t; description: eta """
IsoTrack_fromPV = NanoAODQuantity("IsoTrack_fromPV")
"""dtype: Short_t; description: isolated track comes from PV """
IsoTrack_isFromLostTrack = NanoAODQuantity("IsoTrack_isFromLostTrack")
"""dtype: Bool_t; description: if isolated track comes from a lost track """
IsoTrack_isHighPurityTrack = NanoAODQuantity("IsoTrack_isHighPurityTrack")
"""dtype: Bool_t; description: track is high purity """
IsoTrack_isPFcand = NanoAODQuantity("IsoTrack_isPFcand")
"""dtype: Bool_t; description: if isolated track is a PF candidate """
IsoTrack_pdgId = NanoAODQuantity("IsoTrack_pdgId")
"""dtype: Int_t; description: PDG id of PF cand """
IsoTrack_phi = NanoAODQuantity("IsoTrack_phi")
"""dtype: Float_t; description: phi """
IsoTrack_pt = NanoAODQuantity("IsoTrack_pt")
"""dtype: Float_t; description: pt """

IsoTrack_miniPFRelIso_all = NanoAODQuantity("IsoTrack_miniPFRelIso_all")
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections) """
IsoTrack_miniPFRelIso_chg = NanoAODQuantity("IsoTrack_miniPFRelIso_chg")
"""dtype: Float_t; description: mini PF relative isolation, charged component """

IsoTrack_pfRelIso03_all = NanoAODQuantity("IsoTrack_pfRelIso03_all")
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (deltaBeta corrections) """
IsoTrack_pfRelIso03_chg = NanoAODQuantity("IsoTrack_pfRelIso03_chg")
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component """

nJet = NanoAODQuantity("nJet")
"""dtype: Int_t; description: slimmedJetsPuppi, i.e. ak4 PFJets Puppi with JECs applied, after basic selection (pt > 15) """
Jet_PNetRegPtRawCorr = NanoAODQuantity("Jet_PNetRegPtRawCorr")
"""dtype: Float_t; description: ParticleNet universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT """
Jet_PNetRegPtRawCorrNeutrino = NanoAODQuantity("Jet_PNetRegPtRawCorrNeutrino")
"""dtype: Float_t; description: ParticleNet universal flavor-aware pT regression neutrino correction, relative to visible. To apply full regression, multiply raw jet pT by both PNetRegPtRawCorr and PNetRegPtRawCorrNeutrino. """
Jet_PNetRegPtRawRes = NanoAODQuantity("Jet_PNetRegPtRawRes")
"""dtype: Float_t; description: ParticleNet universal flavor-aware jet pT resolution estimator, (q84 - q16)/2 """
Jet_UParTAK4RegPtRawCorr = NanoAODQuantity("Jet_UParTAK4RegPtRawCorr")
"""dtype: Float_t; description: UnifiedParT universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT """
Jet_UParTAK4RegPtRawCorrNeutrino = NanoAODQuantity("Jet_UParTAK4RegPtRawCorrNeutrino")
"""dtype: Float_t; description: UnifiedParT universal flavor-aware pT regression neutrino correction, relative to visible. Correction relative to raw jet pT """
Jet_UParTAK4RegPtRawRes = NanoAODQuantity("Jet_UParTAK4RegPtRawRes")
"""dtype: Float_t; description: UnifiedParT universal flavor-aware jet pT resolution estimator, (q84 - q16)/2 """
Jet_UParTAK4V1RegPtRawCorr = NanoAODQuantity("Jet_UParTAK4V1RegPtRawCorr")
"""dtype: Float_t; description: UnifiedParT V1 universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT """
Jet_UParTAK4V1RegPtRawCorrNeutrino = NanoAODQuantity(
    "Jet_UParTAK4V1RegPtRawCorrNeutrino"
)
"""dtype: Float_t; description: UnifiedParT V1 universal flavor-aware pT regression neutrino correction, relative to visible. Correction relative to raw jet pT """
Jet_UParTAK4V1RegPtRawRes = NanoAODQuantity("Jet_UParTAK4V1RegPtRawRes")
"""dtype: Float_t; description: UnifiedParT V1 universal flavor-aware jet pT resolution estimator, (q84 - q16)/2 """
Jet_area = NanoAODQuantity("Jet_area")
"""dtype: Float_t; description: jet catchment area, for JECs """
Jet_btagDeepFlavB = NanoAODQuantity("Jet_btagDeepFlavB")
"""dtype: Float_t; description: DeepJet b+bb+lepb tag discriminator """
Jet_btagDeepFlavCvB = NanoAODQuantity("Jet_btagDeepFlavCvB")
"""dtype: Float_t; description: DeepJet c vs b+bb+lepb discriminator """
Jet_btagDeepFlavCvL = NanoAODQuantity("Jet_btagDeepFlavCvL")
"""dtype: Float_t; description: DeepJet c vs uds+g discriminator """
Jet_btagDeepFlavQG = NanoAODQuantity("Jet_btagDeepFlavQG")
"""dtype: Float_t; description: DeepJet g vs uds discriminator """
Jet_btagPNetB = NanoAODQuantity("Jet_btagPNetB")
"""dtype: Float_t; description: ParticleNet b vs. udscg """
Jet_btagPNetCvB = NanoAODQuantity("Jet_btagPNetCvB")
"""dtype: Float_t; description: ParticleNet c vs. b """
Jet_btagPNetCvL = NanoAODQuantity("Jet_btagPNetCvL")
"""dtype: Float_t; description: ParticleNet c vs. udsg """
Jet_btagPNetCvNotB = NanoAODQuantity("Jet_btagPNetCvNotB")
"""dtype: Float_t; description: ParticleNet C vs notB """
Jet_btagPNetQvG = NanoAODQuantity("Jet_btagPNetQvG")
"""dtype: Float_t; description: ParticleNet q (udsbc) vs. g """
Jet_btagPNetTauVJet = NanoAODQuantity("Jet_btagPNetTauVJet")
"""dtype: Float_t; description: ParticleNet tau vs. jet """
Jet_btagUParTAK4B = NanoAODQuantity("Jet_btagUParTAK4B")
"""dtype: Float_t; description: UnifiedParT b vs. udscg """
Jet_btagUParTAK4CvB = NanoAODQuantity("Jet_btagUParTAK4CvB")
"""dtype: Float_t; description: UnifiedParT c vs. b """
Jet_btagUParTAK4CvL = NanoAODQuantity("Jet_btagUParTAK4CvL")
"""dtype: Float_t; description: UnifiedParT c vs. udsg """
Jet_btagUParTAK4CvNotB = NanoAODQuantity("Jet_btagUParTAK4CvNotB")
"""dtype: Float_t; description: UnifiedParT c vs. not b """
Jet_btagUParTAK4Ele = NanoAODQuantity("Jet_btagUParTAK4Ele")
"""dtype: Float_t; description: UnifiedParT electron raw score """
Jet_btagUParTAK4Mu = NanoAODQuantity("Jet_btagUParTAK4Mu")
"""dtype: Float_t; description: UnifiedParT muon raw score """
Jet_btagUParTAK4QvG = NanoAODQuantity("Jet_btagUParTAK4QvG")
"""dtype: Float_t; description: UnifiedParT q (uds) vs. g """
Jet_btagUParTAK4SvCB = NanoAODQuantity("Jet_btagUParTAK4SvCB")
"""dtype: Float_t; description: UnifiedParT s vs. bc """
Jet_btagUParTAK4SvUDG = NanoAODQuantity("Jet_btagUParTAK4SvUDG")
"""dtype: Float_t; description: UnifiedParT s vs. udg """
Jet_btagUParTAK4TauVJet = NanoAODQuantity("Jet_btagUParTAK4TauVJet")
"""dtype: Float_t; description: UnifiedParT tau vs. jet """
Jet_btagUParTAK4UDG = NanoAODQuantity("Jet_btagUParTAK4UDG")
"""dtype: Float_t; description: UnifiedParT u+d+g raw score """
Jet_btagUParTAK4probb = NanoAODQuantity("Jet_btagUParTAK4probb")
"""dtype: Float_t; description: UnifiedParT b raw score """
Jet_btagUParTAK4probbb = NanoAODQuantity("Jet_btagUParTAK4probbb")
"""dtype: Float_t; description: UnifiedParT bb raw score """
Jet_chEmEF = NanoAODQuantity("Jet_chEmEF")
"""dtype: Float_t; description: charged Electromagnetic Energy Fraction """
Jet_chHEF = NanoAODQuantity("Jet_chHEF")
"""dtype: Float_t; description: charged Hadron Energy Fraction """
Jet_chMultiplicity = NanoAODQuantity("Jet_chMultiplicity")
"""dtype: UChar_t; description: (Puppi-weighted) Number of charged particles in the jet """
Jet_electronIdx1 = NanoAODQuantity("Jet_electronIdx1")
"""dtype: Short_t; description: index of first matching electron """
Jet_electronIdx2 = NanoAODQuantity("Jet_electronIdx2")
"""dtype: Short_t; description: index of second matching electron """
Jet_eta = NanoAODQuantity("Jet_eta")
"""dtype: Float_t; description: eta """
Jet_genJetIdx = NanoAODQuantity("Jet_genJetIdx")
"""dtype: Short_t; description: index of matched gen jet """
Jet_hadronFlavour = NanoAODQuantity("Jet_hadronFlavour")
"""dtype: UChar_t; description: flavour from hadron ghost clustering """
Jet_hfEmEF = NanoAODQuantity("Jet_hfEmEF")
"""dtype: Float_t; description: electromagnetic Energy Fraction in HF """
Jet_hfHEF = NanoAODQuantity("Jet_hfHEF")
"""dtype: Float_t; description: hadronic Energy Fraction in HF """
Jet_hfadjacentEtaStripsSize = NanoAODQuantity("Jet_hfadjacentEtaStripsSize")
"""dtype: Int_t; description: eta size of the strips next to the central tower strip in HF (noise discriminating variable) """
Jet_hfcentralEtaStripSize = NanoAODQuantity("Jet_hfcentralEtaStripSize")
"""dtype: Int_t; description: eta size of the central tower strip in HF (noise discriminating variable) """
Jet_hfsigmaEtaEta = NanoAODQuantity("Jet_hfsigmaEtaEta")
"""dtype: Float_t; description: sigmaEtaEta for HF jets (noise discriminating variable) """
Jet_hfsigmaPhiPhi = NanoAODQuantity("Jet_hfsigmaPhiPhi")
"""dtype: Float_t; description: sigmaPhiPhi for HF jets (noise discriminating variable) """
Jet_mass = NanoAODQuantity("Jet_mass")
"""dtype: Float_t; description: mass """
Jet_muEF = NanoAODQuantity("Jet_muEF")
"""dtype: Float_t; description: muon Energy Fraction """
Jet_muonIdx1 = NanoAODQuantity("Jet_muonIdx1")
"""dtype: Short_t; description: index of first matching muon """
Jet_muonIdx2 = NanoAODQuantity("Jet_muonIdx2")
"""dtype: Short_t; description: index of second matching muon """
Jet_muonSubtrDeltaEta = NanoAODQuantity("Jet_muonSubtrDeltaEta")
"""dtype: Float_t; description: muon-subtracted raw eta - eta """
Jet_muonSubtrDeltaPhi = NanoAODQuantity("Jet_muonSubtrDeltaPhi")
"""dtype: Float_t; description: muon-subtracted raw phi - phi """
Jet_muonSubtrFactor = NanoAODQuantity("Jet_muonSubtrFactor")
"""dtype: Float_t; description: 1-(muon-subtracted raw pt)/(raw pt) """
Jet_nConstituents = NanoAODQuantity("Jet_nConstituents")
"""dtype: UChar_t; description: Number of particles in the jet """
Jet_nElectrons = NanoAODQuantity("Jet_nElectrons")
"""dtype: UChar_t; description: number of electrons in the jet """
Jet_nMuons = NanoAODQuantity("Jet_nMuons")
"""dtype: UChar_t; description: number of muons in the jet """
Jet_nSVs = NanoAODQuantity("Jet_nSVs")
"""dtype: UChar_t; description: number of secondary vertices in the jet """
Jet_neEmEF = NanoAODQuantity("Jet_neEmEF")
"""dtype: Float_t; description: neutral Electromagnetic Energy Fraction """
Jet_neHEF = NanoAODQuantity("Jet_neHEF")
"""dtype: Float_t; description: neutral Hadron Energy Fraction """
Jet_neMultiplicity = NanoAODQuantity("Jet_neMultiplicity")
"""dtype: UChar_t; description: (Puppi-weighted) Number of neutral particles in the jet """
Jet_partonFlavour = NanoAODQuantity("Jet_partonFlavour")
"""dtype: Short_t; description: flavour from parton matching """
Jet_phi = NanoAODQuantity("Jet_phi")
"""dtype: Float_t; description: phi """
Jet_pt = NanoAODQuantity("Jet_pt")
"""dtype: Float_t; description: pt """
Jet_puIdDisc = NanoAODQuantity("Jet_puIdDisc")
"""dtype: Float_t; description: Pileup ID BDT discriminant with 133X Winter24 PuppiV18 training """
Jet_rawFactor = NanoAODQuantity("Jet_rawFactor")
"""dtype: Float_t; description: 1 - Factor to get back to raw pT """
Jet_svIdx1 = NanoAODQuantity("Jet_svIdx1")
"""dtype: Short_t; description: index of first matching secondary vertex """
Jet_svIdx2 = NanoAODQuantity("Jet_svIdx2")
"""dtype: Short_t; description: index of second matching secondary vertex """

L1_AlwaysTrue = NanoAODQuantity("L1_AlwaysTrue")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BptxMinus = NanoAODQuantity("L1_BptxMinus")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BptxOR = NanoAODQuantity("L1_BptxOR")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BptxPlus = NanoAODQuantity("L1_BptxPlus")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BptxXOR = NanoAODQuantity("L1_BptxXOR")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 = NanoAODQuantity(
    "L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG11_er1p2_dR_Max0p6 = NanoAODQuantity("L1_DoubleEG11_er1p2_dR_Max0p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau28er2p1 = NanoAODQuantity("L1_DoubleIsoTau28er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau30er2p1 = NanoAODQuantity("L1_DoubleIsoTau30er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau32er2p1 = NanoAODQuantity("L1_DoubleIsoTau32er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau32er2p1_Mass_Max80 = NanoAODQuantity("L1_DoubleIsoTau32er2p1_Mass_Max80")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau34er2p1 = NanoAODQuantity("L1_DoubleIsoTau34er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau35er2p1 = NanoAODQuantity("L1_DoubleIsoTau35er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau36er2p1 = NanoAODQuantity("L1_DoubleIsoTau36er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet100er2p3_dEta_Max1p6 = NanoAODQuantity("L1_DoubleJet100er2p3_dEta_Max1p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet100er2p5 = NanoAODQuantity("L1_DoubleJet100er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet112er2p3_dEta_Max1p6 = NanoAODQuantity("L1_DoubleJet112er2p3_dEta_Max1p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet120er2p5 = NanoAODQuantity("L1_DoubleJet120er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet120er2p5_Mu3_dR_Max0p8 = NanoAODQuantity(
    "L1_DoubleJet120er2p5_Mu3_dR_Max0p8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet150er2p5 = NanoAODQuantity("L1_DoubleJet150er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet16er2p5_Mu3_dR_Max0p4 = NanoAODQuantity("L1_DoubleJet16er2p5_Mu3_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet35er2p5_Mu3_dR_Max0p4 = NanoAODQuantity("L1_DoubleJet35er2p5_Mu3_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet40er2p5 = NanoAODQuantity("L1_DoubleJet40er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet60er2p5_Mu3_dR_Max0p4 = NanoAODQuantity("L1_DoubleJet60er2p5_Mu3_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet80er2p5_Mu3_dR_Max0p4 = NanoAODQuantity("L1_DoubleJet80er2p5_Mu3_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleLLPJet40 = NanoAODQuantity("L1_DoubleLLPJet40")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleLooseIsoEG22er2p1 = NanoAODQuantity("L1_DoubleLooseIsoEG22er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleLooseIsoEG24er2p1 = NanoAODQuantity("L1_DoubleLooseIsoEG24er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu18er2p1_SQ = NanoAODQuantity("L1_DoubleMu18er2p1_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6 = NanoAODQuantity("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6 = NanoAODQuantity("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu6_Upt6_SQ_er2p0 = NanoAODQuantity("L1_DoubleMu6_Upt6_SQ_er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu7_Upt7_SQ_er2p0 = NanoAODQuantity("L1_DoubleMu7_Upt7_SQ_er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu9_SQ = NanoAODQuantity("L1_DoubleMu9_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleTau70er2p1 = NanoAODQuantity("L1_DoubleTau70er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETM120 = NanoAODQuantity("L1_ETM120")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETM150 = NanoAODQuantity("L1_ETM150")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF100 = NanoAODQuantity("L1_ETMHF100")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF100_HTT60er = NanoAODQuantity("L1_ETMHF100_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF110 = NanoAODQuantity("L1_ETMHF110")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF110_HTT60er = NanoAODQuantity("L1_ETMHF110_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF120 = NanoAODQuantity("L1_ETMHF120")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF120_HTT60er = NanoAODQuantity("L1_ETMHF120_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF130 = NanoAODQuantity("L1_ETMHF130")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF130_HTT60er = NanoAODQuantity("L1_ETMHF130_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF140 = NanoAODQuantity("L1_ETMHF140")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF150 = NanoAODQuantity("L1_ETMHF150")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF70 = NanoAODQuantity("L1_ETMHF70")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF70_HTT60er = NanoAODQuantity("L1_ETMHF70_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETT2000 = NanoAODQuantity("L1_ETT2000")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_FirstBunchAfterTrain = NanoAODQuantity("L1_FirstBunchAfterTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_FirstBunchBeforeTrain = NanoAODQuantity("L1_FirstBunchBeforeTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_FirstBunchInTrain = NanoAODQuantity("L1_FirstBunchInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_FirstCollisionInOrbit = NanoAODQuantity("L1_FirstCollisionInOrbit")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_FirstCollisionInTrain = NanoAODQuantity("L1_FirstCollisionInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTMHF100 = NanoAODQuantity("L1_HTMHF100")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTMHF120 = NanoAODQuantity("L1_HTMHF120")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTMHF125 = NanoAODQuantity("L1_HTMHF125")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTMHF130 = NanoAODQuantity("L1_HTMHF130")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTMHF150 = NanoAODQuantity("L1_HTMHF150")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT120_SingleLLPJet40 = NanoAODQuantity("L1_HTT120_SingleLLPJet40")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT120er = NanoAODQuantity("L1_HTT120er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT160_SingleLLPJet50 = NanoAODQuantity("L1_HTT160_SingleLLPJet50")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT160er = NanoAODQuantity("L1_HTT160er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT200_SingleLLPJet60 = NanoAODQuantity("L1_HTT200_SingleLLPJet60")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT200er = NanoAODQuantity("L1_HTT200er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT240_SingleLLPJet70 = NanoAODQuantity("L1_HTT240_SingleLLPJet70")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT255er = NanoAODQuantity("L1_HTT255er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT280er = NanoAODQuantity("L1_HTT280er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT280er_QuadJet_70_55_40_35_er2p5 = NanoAODQuantity(
    "L1_HTT280er_QuadJet_70_55_40_35_er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT320er = NanoAODQuantity("L1_HTT320er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT360er = NanoAODQuantity("L1_HTT360er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT400er = NanoAODQuantity("L1_HTT400er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT450er = NanoAODQuantity("L1_HTT450er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_IsoEG32er2p5_Mt40 = NanoAODQuantity("L1_IsoEG32er2p5_Mt40")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_IsoTau52er2p1_QuadJet36er2p5 = NanoAODQuantity("L1_IsoTau52er2p1_QuadJet36er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_IsolatedBunch = NanoAODQuantity("L1_IsolatedBunch")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LastBunchInTrain = NanoAODQuantity("L1_LastBunchInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LastCollisionInTrain = NanoAODQuantity("L1_LastCollisionInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG14er2p5_HTT200er = NanoAODQuantity("L1_LooseIsoEG14er2p5_HTT200er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG16er2p5_HTT200er = NanoAODQuantity("L1_LooseIsoEG16er2p5_HTT200er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = NanoAODQuantity(
    "L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_MinimumBiasHF0 = NanoAODQuantity("L1_MinimumBiasHF0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_MinimumBiasHF0_AND_BptxAND = NanoAODQuantity("L1_MinimumBiasHF0_AND_BptxAND")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6 = NanoAODQuantity(
    "L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu12_HTT150er = NanoAODQuantity("L1_Mu12_HTT150er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6 = NanoAODQuantity(
    "L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu14_HTT150er = NanoAODQuantity("L1_Mu14_HTT150er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu20_EG10er2p5 = NanoAODQuantity("L1_Mu20_EG10er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_NotBptxOR = NanoAODQuantity("L1_NotBptxOR")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_QuadJet60er2p5 = NanoAODQuantity("L1_QuadJet60er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0 = NanoAODQuantity(
    "L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SecondBunchInTrain = NanoAODQuantity("L1_SecondBunchInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SecondLastBunchInTrain = NanoAODQuantity("L1_SecondLastBunchInTrain")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG10er2p5 = NanoAODQuantity("L1_SingleEG10er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG15er2p5 = NanoAODQuantity("L1_SingleEG15er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG26er2p5 = NanoAODQuantity("L1_SingleEG26er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG28_FWD2p5 = NanoAODQuantity("L1_SingleEG28_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG28er1p5 = NanoAODQuantity("L1_SingleEG28er1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG28er2p1 = NanoAODQuantity("L1_SingleEG28er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG28er2p5 = NanoAODQuantity("L1_SingleEG28er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG34er2p5 = NanoAODQuantity("L1_SingleEG34er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG36er2p5 = NanoAODQuantity("L1_SingleEG36er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG38er2p5 = NanoAODQuantity("L1_SingleEG38er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG40er2p5 = NanoAODQuantity("L1_SingleEG40er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG42er2p5 = NanoAODQuantity("L1_SingleEG42er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG45er2p5 = NanoAODQuantity("L1_SingleEG45er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG50 = NanoAODQuantity("L1_SingleEG50")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG60 = NanoAODQuantity("L1_SingleEG60")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleEG8er2p5 = NanoAODQuantity("L1_SingleEG8er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG24er2p1 = NanoAODQuantity("L1_SingleIsoEG24er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG26er2p1 = NanoAODQuantity("L1_SingleIsoEG26er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG26er2p5 = NanoAODQuantity("L1_SingleIsoEG26er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG28_FWD2p5 = NanoAODQuantity("L1_SingleIsoEG28_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG28er1p5 = NanoAODQuantity("L1_SingleIsoEG28er1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG28er2p1 = NanoAODQuantity("L1_SingleIsoEG28er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG28er2p5 = NanoAODQuantity("L1_SingleIsoEG28er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG30er2p1 = NanoAODQuantity("L1_SingleIsoEG30er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG30er2p5 = NanoAODQuantity("L1_SingleIsoEG30er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG32er2p1 = NanoAODQuantity("L1_SingleIsoEG32er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG32er2p5 = NanoAODQuantity("L1_SingleIsoEG32er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleIsoEG34er2p5 = NanoAODQuantity("L1_SingleIsoEG34er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet10erHE = NanoAODQuantity("L1_SingleJet10erHE")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet120er1p3 = NanoAODQuantity("L1_SingleJet120er1p3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet120er2p5 = NanoAODQuantity("L1_SingleJet120er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet12erHE = NanoAODQuantity("L1_SingleJet12erHE")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet140er2p5 = NanoAODQuantity("L1_SingleJet140er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet160er2p5 = NanoAODQuantity("L1_SingleJet160er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet180 = NanoAODQuantity("L1_SingleJet180")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet180er2p5 = NanoAODQuantity("L1_SingleJet180er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet200 = NanoAODQuantity("L1_SingleJet200")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet35er1p3 = NanoAODQuantity("L1_SingleJet35er1p3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet35er2p5 = NanoAODQuantity("L1_SingleJet35er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet43er2p5_NotBptxOR_3BX = NanoAODQuantity("L1_SingleJet43er2p5_NotBptxOR_3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet46er2p5_NotBptxOR_3BX = NanoAODQuantity("L1_SingleJet46er2p5_NotBptxOR_3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet60 = NanoAODQuantity("L1_SingleJet60")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet60_FWD2p5 = NanoAODQuantity("L1_SingleJet60_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet8erHE = NanoAODQuantity("L1_SingleJet8erHE")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet90 = NanoAODQuantity("L1_SingleJet90")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet90_FWD2p5 = NanoAODQuantity("L1_SingleJet90_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG26er1p5 = NanoAODQuantity("L1_SingleLooseIsoEG26er1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG26er2p5 = NanoAODQuantity("L1_SingleLooseIsoEG26er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG28_FWD2p5 = NanoAODQuantity("L1_SingleLooseIsoEG28_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG28er1p5 = NanoAODQuantity("L1_SingleLooseIsoEG28er1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG28er2p1 = NanoAODQuantity("L1_SingleLooseIsoEG28er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG28er2p5 = NanoAODQuantity("L1_SingleLooseIsoEG28er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG30er1p5 = NanoAODQuantity("L1_SingleLooseIsoEG30er1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleLooseIsoEG30er2p5 = NanoAODQuantity("L1_SingleLooseIsoEG30er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu10_SQ14_BMTF = NanoAODQuantity("L1_SingleMu10_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu11_SQ14_BMTF = NanoAODQuantity("L1_SingleMu11_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu15_DQ = NanoAODQuantity("L1_SingleMu15_DQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu18 = NanoAODQuantity("L1_SingleMu18")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu20 = NanoAODQuantity("L1_SingleMu20")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu25 = NanoAODQuantity("L1_SingleMu25")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu3 = NanoAODQuantity("L1_SingleMu3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu5 = NanoAODQuantity("L1_SingleMu5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu5_SQ14_BMTF = NanoAODQuantity("L1_SingleMu5_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu6_SQ14_BMTF = NanoAODQuantity("L1_SingleMu6_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu8_SQ14_BMTF = NanoAODQuantity("L1_SingleMu8_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu9_SQ14_BMTF = NanoAODQuantity("L1_SingleMu9_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleTau120er2p1 = NanoAODQuantity("L1_SingleTau120er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleTau130er2p1 = NanoAODQuantity("L1_SingleTau130er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu3 = NanoAODQuantity("L1_TripleMu3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu3_SQ = NanoAODQuantity("L1_TripleMu3_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TwoMuShower_Loose = NanoAODQuantity("L1_TwoMuShower_Loose")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_UnpairedBunchBptxMinus = NanoAODQuantity("L1_UnpairedBunchBptxMinus")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_UnpairedBunchBptxPlus = NanoAODQuantity("L1_UnpairedBunchBptxPlus")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ZeroBias = NanoAODQuantity("L1_ZeroBias")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ZeroBias_copy = NanoAODQuantity("L1_ZeroBias_copy")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_AXO_Loose = NanoAODQuantity("L1_AXO_Loose")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_AXO_Nominal = NanoAODQuantity("L1_AXO_Nominal")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_AXO_Tight = NanoAODQuantity("L1_AXO_Tight")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_AXO_VLoose = NanoAODQuantity("L1_AXO_VLoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_AXO_VTight = NanoAODQuantity("L1_AXO_VTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_BPTX_NotOR_VME = NanoAODQuantity("L1_BPTX_NotOR_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_RefAND_VME = NanoAODQuantity("L1_BPTX_RefAND_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_BPTX_AND_Ref1_VME = NanoAODQuantity("L1_BPTX_AND_Ref1_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_AND_Ref3_VME = NanoAODQuantity("L1_BPTX_AND_Ref3_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_AND_Ref4_VME = NanoAODQuantity("L1_BPTX_AND_Ref4_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_BPTX_BeamGas_B1_VME = NanoAODQuantity("L1_BPTX_BeamGas_B1_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_BeamGas_B2_VME = NanoAODQuantity("L1_BPTX_BeamGas_B2_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_BeamGas_Ref1_VME = NanoAODQuantity("L1_BPTX_BeamGas_Ref1_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_BeamGas_Ref2_VME = NanoAODQuantity("L1_BPTX_BeamGas_Ref2_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_BPTX_OR_Ref3_VME = NanoAODQuantity("L1_BPTX_OR_Ref3_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_BPTX_OR_Ref4_VME = NanoAODQuantity("L1_BPTX_OR_Ref4_VME")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_CICADA_Loose = NanoAODQuantity("L1_CICADA_Loose")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_CICADA_Medium = NanoAODQuantity("L1_CICADA_Medium")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_CICADA_Tight = NanoAODQuantity("L1_CICADA_Tight")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_CICADA_VLoose = NanoAODQuantity("L1_CICADA_VLoose")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_CICADA_VTight = NanoAODQuantity("L1_CICADA_VTight")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleEG_15_10_er2p5 = NanoAODQuantity("L1_DoubleEG_15_10_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_20_10_er2p5 = NanoAODQuantity("L1_DoubleEG_20_10_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_22_10_er2p5 = NanoAODQuantity("L1_DoubleEG_22_10_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_27_14_er2p5 = NanoAODQuantity("L1_DoubleEG_27_14_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_LooseIso16_LooseIso12_er1p5 = NanoAODQuantity(
    "L1_DoubleEG_LooseIso16_LooseIso12_er1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_LooseIso18_LooseIso12_er1p5 = NanoAODQuantity(
    "L1_DoubleEG_LooseIso18_LooseIso12_er1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_LooseIso20_LooseIso12_er1p5 = NanoAODQuantity(
    "L1_DoubleEG_LooseIso20_LooseIso12_er1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleEG8er2p5_HTT280er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT280er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG8er2p5_HTT300er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT300er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG8er2p5_HTT320er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT320er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleEG_25_12_er2p5 = NanoAODQuantity("L1_DoubleEG_25_12_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_25_14_er2p5 = NanoAODQuantity("L1_DoubleEG_25_14_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleEG_LooseIso22_12_er2p5 = NanoAODQuantity("L1_DoubleEG_LooseIso22_12_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_LooseIso22_LooseIso12_er1p5 = NanoAODQuantity(
    "L1_DoubleEG_LooseIso22_LooseIso12_er1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleEG_LooseIso25_12_er2p5 = NanoAODQuantity("L1_DoubleEG_LooseIso25_12_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleEG_LooseIso25_LooseIso12_er1p5 = NanoAODQuantity(
    "L1_DoubleEG_LooseIso25_LooseIso12_er1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5 = NanoAODQuantity(
    "L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5 = NanoAODQuantity(
    "L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5 = NanoAODQuantity(
    "L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 = NanoAODQuantity(
    "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 = NanoAODQuantity(
    "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5 = NanoAODQuantity(
    "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5 = NanoAODQuantity(
    "L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2 = NanoAODQuantity(
    "L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5 = NanoAODQuantity(
    "L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2 = NanoAODQuantity(
    "L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet_110_35_DoubleJet35_Mass_Min800 = NanoAODQuantity(
    "L1_DoubleJet_110_35_DoubleJet35_Mass_Min800"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet_110_35_DoubleJet35_Mass_Min850 = NanoAODQuantity(
    "L1_DoubleJet_110_35_DoubleJet35_Mass_Min850"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50 = NanoAODQuantity(
    "L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50 = NanoAODQuantity(
    "L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65 = NanoAODQuantity(
    "L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65 = NanoAODQuantity(
    "L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ = NanoAODQuantity(
    "L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ = NanoAODQuantity(
    "L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu_12_5 = NanoAODQuantity("L1_DoubleMu_12_5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0 = NanoAODQuantity("L1_DoubleMu0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_Mass_Min1 = NanoAODQuantity("L1_DoubleMu0_Mass_Min1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_OQ = NanoAODQuantity("L1_DoubleMu0_OQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_SQ = NanoAODQuantity("L1_DoubleMu0_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_SQ_OS = NanoAODQuantity("L1_DoubleMu0_SQ_OS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_Upt15_Upt7 = NanoAODQuantity("L1_DoubleMu0_Upt15_Upt7")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_Upt5_Upt5 = NanoAODQuantity("L1_DoubleMu0_Upt5_Upt5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_Upt7_SQ_er2p0 = NanoAODQuantity("L1_DoubleMu0_Upt7_SQ_er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_Upt8_SQ_er2p0 = NanoAODQuantity("L1_DoubleMu0_Upt8_SQ_er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8 = NanoAODQuantity(
    "L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0_Upt6_IP_Min1_Upt4 = NanoAODQuantity("L1_DoubleMu0_Upt6_IP_Min1_Upt4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0_Upt6_SQ_er2p0 = NanoAODQuantity("L1_DoubleMu0_Upt6_SQ_er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6 = NanoAODQuantity(
    "L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2 = NanoAODQuantity(
    "L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0er1p5_SQ = NanoAODQuantity("L1_DoubleMu0er1p5_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er1p5_SQ_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er1p5_SQ_dR_Max1p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0er1p5_SQ_OS = NanoAODQuantity("L1_DoubleMu0er1p5_SQ_OS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2 = NanoAODQuantity(
    "L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2 = NanoAODQuantity(
    "L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5 = NanoAODQuantity(
    "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6 = NanoAODQuantity(
    "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu0er2p0_SQ_dEta_Max1p5 = NanoAODQuantity("L1_DoubleMu0er2p0_SQ_dEta_Max1p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu0er2p0_SQ_dEta_Max1p6 = NanoAODQuantity("L1_DoubleMu0er2p0_SQ_dEta_Max1p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20 = NanoAODQuantity(
    "L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5 = NanoAODQuantity(
    "L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu3_SQ_HTT220er = NanoAODQuantity("L1_DoubleMu3_SQ_HTT220er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu3_SQ_ETMHF30_HTT60er = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF30_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5 = NanoAODQuantity(
    "L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu3_SQ_ETMHF40_HTT60er = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF40_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5 = NanoAODQuantity(
    "L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu3_SQ_ETMHF50_HTT60er = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF50_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5 = NanoAODQuantity(
    "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5 = NanoAODQuantity(
    "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu4_SQ_EG9er2p5 = NanoAODQuantity("L1_DoubleMu4_SQ_EG9er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu4_SQ_OS = NanoAODQuantity("L1_DoubleMu4_SQ_OS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu4_SQ_OS_dR_Max1p2 = NanoAODQuantity("L1_DoubleMu4_SQ_OS_dR_Max1p2")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu4p5_SQ_OS = NanoAODQuantity("L1_DoubleMu4p5_SQ_OS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu4p5_SQ_OS_dR_Max1p2 = NanoAODQuantity("L1_DoubleMu4p5_SQ_OS_dR_Max1p2")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu4p5er2p0_SQ_OS = NanoAODQuantity("L1_DoubleMu4p5er2p0_SQ_OS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18 = NanoAODQuantity(
    "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 = NanoAODQuantity(
    "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20 = NanoAODQuantity(
    "L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu5_SQ_EG9er2p5 = NanoAODQuantity("L1_DoubleMu5_SQ_EG9er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu5_SQ_OS_dR_Max1p6 = NanoAODQuantity("L1_DoubleMu5_SQ_OS_dR_Max1p6")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu8_SQ = NanoAODQuantity("L1_DoubleMu8_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu8_Upt8_SQ_er2p0 = NanoAODQuantity("L1_DoubleMu8_Upt8_SQ_er2p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu_15_5_SQ = NanoAODQuantity("L1_DoubleMu_15_5_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_DoubleMu_15_7 = NanoAODQuantity("L1_DoubleMu_15_7")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu_15_7_Mass_Min1 = NanoAODQuantity("L1_DoubleMu_15_7_Mass_Min1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_DoubleMu_15_7_SQ = NanoAODQuantity("L1_DoubleMu_15_7_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_ETMHF80 = NanoAODQuantity("L1_ETMHF80")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF80_HTT60er = NanoAODQuantity("L1_ETMHF80_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1 = NanoAODQuantity(
    "L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6 = NanoAODQuantity(
    "L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_ETMHF90 = NanoAODQuantity("L1_ETMHF90")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF90_HTT60er = NanoAODQuantity("L1_ETMHF90_HTT60er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1 = NanoAODQuantity(
    "L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6 = NanoAODQuantity(
    "L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_FinalOR_BXmin1 = NanoAODQuantity("L1_FinalOR_BXmin1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_FinalOR_BXmin2 = NanoAODQuantity("L1_FinalOR_BXmin2")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_HCAL_LaserMon_Trig = NanoAODQuantity("L1_HCAL_LaserMon_Trig")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HCAL_LaserMon_Veto = NanoAODQuantity("L1_HCAL_LaserMon_Veto")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_HTT320er_QuadJet_70_55_40_40_er2p5 = NanoAODQuantity(
    "L1_HTT320er_QuadJet_70_55_40_40_er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 = NanoAODQuantity(
    "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = NanoAODQuantity(
    "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3 = NanoAODQuantity(
    "L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3 = NanoAODQuantity(
    "L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_LooseIsoEG24er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG24er2p1_HTT100er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3 = NanoAODQuantity(
    "L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_LooseIsoEG26er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG26er2p1_HTT100er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = NanoAODQuantity(
    "L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_LooseIsoEG28er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG28er2p1_HTT100er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = NanoAODQuantity(
    "L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu18er2p1_Tau24er2p1 = NanoAODQuantity("L1_Mu18er2p1_Tau24er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu18er2p1_Tau26er2p1 = NanoAODQuantity("L1_Mu18er2p1_Tau26er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu18er2p1_Tau26er2p1_Jet55 = NanoAODQuantity("L1_Mu18er2p1_Tau26er2p1_Jet55")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu18er2p1_Tau26er2p1_Jet70 = NanoAODQuantity("L1_Mu18er2p1_Tau26er2p1_Jet70")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu22er2p1_IsoTau30er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau30er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu22er2p1_IsoTau32er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau32er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu22er2p1_IsoTau34er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau34er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu22er2p1_IsoTau40er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau40er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu22er2p1_Tau70er2p1 = NanoAODQuantity("L1_Mu22er2p1_Tau70er2p1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu3_Jet120er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet120er2p5_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu3_Jet16er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet16er2p5_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu3_Jet30er2p5 = NanoAODQuantity("L1_Mu3_Jet30er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu3_Jet60er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet60er2p5_dR_Max0p4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu3er1p5_Jet100er2p5_ETMHF30 = NanoAODQuantity("L1_Mu3er1p5_Jet100er2p5_ETMHF30")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu3er1p5_Jet100er2p5_ETMHF40 = NanoAODQuantity("L1_Mu3er1p5_Jet100er2p5_ETMHF40")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu3er1p5_Jet100er2p5_ETMHF50 = NanoAODQuantity("L1_Mu3er1p5_Jet100er2p5_ETMHF50")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu5_EG23er2p5 = NanoAODQuantity("L1_Mu5_EG23er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu5_LooseIsoEG20er2p5 = NanoAODQuantity("L1_Mu5_LooseIsoEG20er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu6_DoubleEG10er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG10er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu6_DoubleEG12er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG12er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu6_DoubleEG15er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG15er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu6_DoubleEG17er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG17er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu6_HTT240er = NanoAODQuantity("L1_Mu6_HTT240er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu6_HTT250er = NanoAODQuantity("L1_Mu6_HTT250er")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_Mu7_EG20er2p5 = NanoAODQuantity("L1_Mu7_EG20er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu7_EG23er2p5 = NanoAODQuantity("L1_Mu7_EG23er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu7_LooseIsoEG20er2p5 = NanoAODQuantity("L1_Mu7_LooseIsoEG20er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_Mu7_LooseIsoEG23er2p5 = NanoAODQuantity("L1_Mu7_LooseIsoEG23er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_QuadMu0 = NanoAODQuantity("L1_QuadMu0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_QuadMu0_OQ = NanoAODQuantity("L1_QuadMu0_OQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_QuadMu0_SQ = NanoAODQuantity("L1_QuadMu0_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleJet120 = NanoAODQuantity("L1_SingleJet120")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet120_FWD2p5 = NanoAODQuantity("L1_SingleJet120_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet120_FWD3p0 = NanoAODQuantity("L1_SingleJet120_FWD3p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleJet20er2p5_NotBptxOR = NanoAODQuantity("L1_SingleJet20er2p5_NotBptxOR")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet20er2p5_NotBptxOR_3BX = NanoAODQuantity("L1_SingleJet20er2p5_NotBptxOR_3BX")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleJet35 = NanoAODQuantity("L1_SingleJet35")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet35_FWD2p5 = NanoAODQuantity("L1_SingleJet35_FWD2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleJet35_FWD3p0 = NanoAODQuantity("L1_SingleJet35_FWD3p0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu0_BMTF = NanoAODQuantity("L1_SingleMu0_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_DQ = NanoAODQuantity("L1_SingleMu0_DQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_EMTF = NanoAODQuantity("L1_SingleMu0_EMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_OMTF = NanoAODQuantity("L1_SingleMu0_OMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_SQ13_BMTF = NanoAODQuantity("L1_SingleMu0_SQ13_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_SQ14_BMTF = NanoAODQuantity("L1_SingleMu0_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_SQ15_BMTF = NanoAODQuantity("L1_SingleMu0_SQ15_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt15_SQ14_BMTF = NanoAODQuantity("L1_SingleMu0_Upt15_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt20_SQ14_BMTF = NanoAODQuantity("L1_SingleMu0_Upt20_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt25_SQ14_BMTF = NanoAODQuantity("L1_SingleMu0_Upt25_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu0_Upt10 = NanoAODQuantity("L1_SingleMu0_Upt10")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt10_BMTF = NanoAODQuantity("L1_SingleMu0_Upt10_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt10_EMTF = NanoAODQuantity("L1_SingleMu0_Upt10_EMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt10_OMTF = NanoAODQuantity("L1_SingleMu0_Upt10_OMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu0_Upt10_SQ14_BMTF = NanoAODQuantity("L1_SingleMu0_Upt10_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu12_DQ_BMTF = NanoAODQuantity("L1_SingleMu12_DQ_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu12_DQ_EMTF = NanoAODQuantity("L1_SingleMu12_DQ_EMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu12_DQ_OMTF = NanoAODQuantity("L1_SingleMu12_DQ_OMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu22 = NanoAODQuantity("L1_SingleMu22")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_DQ = NanoAODQuantity("L1_SingleMu22_DQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_OQ = NanoAODQuantity("L1_SingleMu22_OQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu22_BMTF = NanoAODQuantity("L1_SingleMu22_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_BMTF_NEG = NanoAODQuantity("L1_SingleMu22_BMTF_NEG")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_BMTF_POS = NanoAODQuantity("L1_SingleMu22_BMTF_POS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu22_EMTF = NanoAODQuantity("L1_SingleMu22_EMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_EMTF_NEG = NanoAODQuantity("L1_SingleMu22_EMTF_NEG")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_EMTF_POS = NanoAODQuantity("L1_SingleMu22_EMTF_POS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu22_OMTF = NanoAODQuantity("L1_SingleMu22_OMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_OMTF_NEG = NanoAODQuantity("L1_SingleMu22_OMTF_NEG")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu22_OMTF_POS = NanoAODQuantity("L1_SingleMu22_OMTF_POS")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMu7 = NanoAODQuantity("L1_SingleMu7")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu7_DQ = NanoAODQuantity("L1_SingleMu7_DQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMu7_SQ14_BMTF = NanoAODQuantity("L1_SingleMu7_SQ14_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMuCosmics = NanoAODQuantity("L1_SingleMuCosmics")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuCosmics_BMTF = NanoAODQuantity("L1_SingleMuCosmics_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuCosmics_EMTF = NanoAODQuantity("L1_SingleMuCosmics_EMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuCosmics_OMTF = NanoAODQuantity("L1_SingleMuCosmics_OMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMuOpen = NanoAODQuantity("L1_SingleMuOpen")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuOpen_BMTF = NanoAODQuantity("L1_SingleMuOpen_BMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuOpen_EMTF = NanoAODQuantity("L1_SingleMuOpen_EMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuOpen_NotBptxOR = NanoAODQuantity("L1_SingleMuOpen_NotBptxOR")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuOpen_OMTF = NanoAODQuantity("L1_SingleMuOpen_OMTF")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuOpen_er1p1_NotBptxOR_3BX = NanoAODQuantity(
    "L1_SingleMuOpen_er1p1_NotBptxOR_3BX"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuOpen_er1p4_NotBptxOR_3BX = NanoAODQuantity(
    "L1_SingleMuOpen_er1p4_NotBptxOR_3BX"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_SingleMuShower_Nominal = NanoAODQuantity("L1_SingleMuShower_Nominal")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_SingleMuShower_Tight = NanoAODQuantity("L1_SingleMuShower_Tight")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TOTEM_1 = NanoAODQuantity("L1_TOTEM_1")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TOTEM_2 = NanoAODQuantity("L1_TOTEM_2")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TOTEM_3 = NanoAODQuantity("L1_TOTEM_3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TOTEM_4 = NanoAODQuantity("L1_TOTEM_4")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleEG_18_17_8_er2p5 = NanoAODQuantity("L1_TripleEG_18_17_8_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleEG_18_18_12_er2p5 = NanoAODQuantity("L1_TripleEG_18_18_12_er2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5 = NanoAODQuantity(
    "L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5 = NanoAODQuantity(
    "L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 = NanoAODQuantity(
    "L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12 = NanoAODQuantity(
    "L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu0 = NanoAODQuantity("L1_TripleMu0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu0_OQ = NanoAODQuantity("L1_TripleMu0_OQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu0_SQ = NanoAODQuantity("L1_TripleMu0_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu_3SQ_2p5SQ_0 = NanoAODQuantity("L1_TripleMu_3SQ_2p5SQ_0")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12 = NanoAODQuantity(
    "L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12 = NanoAODQuantity(
    "L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17 = NanoAODQuantity(
    "L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_5_5_3 = NanoAODQuantity("L1_TripleMu_5_5_3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu_5SQ_3SQ_0OQ = NanoAODQuantity("L1_TripleMu_5SQ_3SQ_0OQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9 = NanoAODQuantity(
    "L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9 = NanoAODQuantity(
    "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu_5_3_3 = NanoAODQuantity("L1_TripleMu_5_3_3")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_5_3_3_SQ = NanoAODQuantity("L1_TripleMu_5_3_3_SQ")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_TripleMu_5_3p5_2p5 = NanoAODQuantity("L1_TripleMu_5_3p5_2p5")
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17 = NanoAODQuantity(
    "L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

L1_UnprefireableEvent_FirstBxInTrain = NanoAODQuantity(
    "L1_UnprefireableEvent_FirstBxInTrain"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """
L1_UnprefireableEvent_TriggerRules = NanoAODQuantity(
    "L1_UnprefireableEvent_TriggerRules"
)
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO) """

LHE_AlphaS = NanoAODQuantity("LHE_AlphaS")
"""dtype: Float_t; description: Per-event alphaS """
LHE_HT = NanoAODQuantity("LHE_HT")
"""dtype: Float_t; description: HT, scalar sum of parton pTs at LHE step """
LHE_HTIncoming = NanoAODQuantity("LHE_HTIncoming")
"""dtype: Float_t; description: HT, scalar sum of parton pTs at LHE step, restricted to partons """
LHE_Nb = NanoAODQuantity("LHE_Nb")
"""dtype: UChar_t; description: Number of b partons at LHE step """
LHE_Nc = NanoAODQuantity("LHE_Nc")
"""dtype: UChar_t; description: Number of c partons at LHE step """
LHE_Nglu = NanoAODQuantity("LHE_Nglu")
"""dtype: UChar_t; description: Number of gluon partons at LHE step """
LHE_Njets = NanoAODQuantity("LHE_Njets")
"""dtype: UChar_t; description: Number of jets (partons) at LHE step """
LHE_NpLO = NanoAODQuantity("LHE_NpLO")
"""dtype: UChar_t; description: number of partons at LO """
LHE_NpNLO = NanoAODQuantity("LHE_NpNLO")
"""dtype: UChar_t; description: number of partons at NLO """
LHE_Nuds = NanoAODQuantity("LHE_Nuds")
"""dtype: UChar_t; description: Number of u,d,s partons at LHE step """
LHE_Vpt = NanoAODQuantity("LHE_Vpt")
"""dtype: Float_t; description: pT of the W or Z boson at LHE step """

nLHEPart = NanoAODQuantity("nLHEPart")
"""dtype: Int_t; description:  """
LHEPart_eta = NanoAODQuantity("LHEPart_eta")
"""dtype: Float_t; description: Pseodorapidity of LHE particles """
LHEPart_firstMotherIdx = NanoAODQuantity("LHEPart_firstMotherIdx")
"""dtype: Short_t; description: Index of this particle's first mother in the LHEPart collection """
LHEPart_incomingpz = NanoAODQuantity("LHEPart_incomingpz")
"""dtype: Float_t; description: Pz of incoming LHE particles """
LHEPart_lastMotherIdx = NanoAODQuantity("LHEPart_lastMotherIdx")
"""dtype: Short_t; description: Index of this particle's last mother in the LHEPart collection """
LHEPart_mass = NanoAODQuantity("LHEPart_mass")
"""dtype: Float_t; description: Mass of LHE particles """
LHEPart_pdgId = NanoAODQuantity("LHEPart_pdgId")
"""dtype: Int_t; description: PDG ID of LHE particles """
LHEPart_phi = NanoAODQuantity("LHEPart_phi")
"""dtype: Float_t; description: Phi of LHE particles """
LHEPart_pt = NanoAODQuantity("LHEPart_pt")
"""dtype: Float_t; description: Pt of LHE particles """
LHEPart_spin = NanoAODQuantity("LHEPart_spin")
"""dtype: Int_t; description: Spin of LHE particles """
LHEPart_status = NanoAODQuantity("LHEPart_status")
"""dtype: Int_t; description: LHE particle status; -1:incoming, 1:outgoing """

nLHEPdfWeight = NanoAODQuantity("nLHEPdfWeight")
"""dtype: Int_t; description:  """
LHEPdfWeight = NanoAODQuantity("LHEPdfWeight")
"""dtype: Float_t; description: LHE pdf variation weights (w_var / w_nominal) for LHA IDs 325300 - 325402 """

nLHEReweightingWeight = NanoAODQuantity("nLHEReweightingWeight")
"""dtype: Int_t; description:  """
LHEReweightingWeight = NanoAODQuantity("LHEReweightingWeight")
"""dtype: Float_t; description:  """

nLHEScaleWeight = NanoAODQuantity("nLHEScaleWeight")
"""dtype: Int_t; description:  """
LHEScaleWeight = NanoAODQuantity("LHEScaleWeight")
"""dtype: Float_t; description: LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0  """

nLowPtElectron = NanoAODQuantity("nLowPtElectron")
"""dtype: Int_t; description: slimmedLowPtElectrons after basic selection (pt > 1. && electronID('ID') > -0.25) """
LowPtElectron_ID = NanoAODQuantity("LowPtElectron_ID")
"""dtype: Float_t; description: ID, BDT (raw) score """
LowPtElectron_charge = NanoAODQuantity("LowPtElectron_charge")
"""dtype: Int_t; description: electric charge """
LowPtElectron_convVeto = NanoAODQuantity("LowPtElectron_convVeto")
"""dtype: Bool_t; description: pass conversion veto """
LowPtElectron_convVtxRadius = NanoAODQuantity("LowPtElectron_convVtxRadius")
"""dtype: Float_t; description: conversion vertex radius (cm) """
LowPtElectron_convWP = NanoAODQuantity("LowPtElectron_convWP")
"""dtype: UChar_t; description: conversion flag bit map: 1=Veto, 2=Loose, 3=Tight """
LowPtElectron_deltaEtaSC = NanoAODQuantity("LowPtElectron_deltaEtaSC")
"""dtype: Float_t; description: delta eta (SC,ele) with sign """
LowPtElectron_dxy = NanoAODQuantity("LowPtElectron_dxy")
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm """
LowPtElectron_dxyErr = NanoAODQuantity("LowPtElectron_dxyErr")
"""dtype: Float_t; description: dxy uncertainty, in cm """
LowPtElectron_dz = NanoAODQuantity("LowPtElectron_dz")
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm """
LowPtElectron_dzErr = NanoAODQuantity("LowPtElectron_dzErr")
"""dtype: Float_t; description: dz uncertainty, in cm """
LowPtElectron_eInvMinusPInv = NanoAODQuantity("LowPtElectron_eInvMinusPInv")
"""dtype: Float_t; description: 1/E_SC - 1/p_trk """
LowPtElectron_electronIdx = NanoAODQuantity("LowPtElectron_electronIdx")
"""dtype: Short_t; description: index of the overlapping PF electron (-1 if none) """
LowPtElectron_energyErr = NanoAODQuantity("LowPtElectron_energyErr")
"""dtype: Float_t; description: energy error of the cluster-track combination """
LowPtElectron_eta = NanoAODQuantity("LowPtElectron_eta")
"""dtype: Float_t; description: eta """
LowPtElectron_genPartFlav = NanoAODQuantity("LowPtElectron_genPartFlav")
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 electrons or photons: 1 = prompt electron (including gamma*->mu mu), 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched """
LowPtElectron_genPartIdx = NanoAODQuantity("LowPtElectron_genPartIdx")
"""dtype: Short_t; description: Index into genParticle list for MC matching to status==1 electrons or photons """
LowPtElectron_hoe = NanoAODQuantity("LowPtElectron_hoe")
"""dtype: Float_t; description: H over E """
LowPtElectron_lostHits = NanoAODQuantity("LowPtElectron_lostHits")
"""dtype: UChar_t; description: number of missing inner hits """
LowPtElectron_mass = NanoAODQuantity("LowPtElectron_mass")
"""dtype: Float_t; description: mass """
LowPtElectron_pdgId = NanoAODQuantity("LowPtElectron_pdgId")
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth) """
LowPtElectron_phi = NanoAODQuantity("LowPtElectron_phi")
"""dtype: Float_t; description: phi """
LowPtElectron_photonIdx = NanoAODQuantity("LowPtElectron_photonIdx")
"""dtype: Short_t; description: index of the first associated photon (-1 if none) """
LowPtElectron_pt = NanoAODQuantity("LowPtElectron_pt")
"""dtype: Float_t; description: pt """
LowPtElectron_ptbiased = NanoAODQuantity("LowPtElectron_ptbiased")
"""dtype: Float_t; description: ElectronSeed, pT- and dxy- dependent BDT (raw) score """
LowPtElectron_r9 = NanoAODQuantity("LowPtElectron_r9")
"""dtype: Float_t; description: R9 of the SC, calculated with full 5x5 region """
LowPtElectron_scEtOverPt = NanoAODQuantity("LowPtElectron_scEtOverPt")
"""dtype: Float_t; description: (SC energy)/pt-1 """
LowPtElectron_sieie = NanoAODQuantity("LowPtElectron_sieie")
"""dtype: Float_t; description: sigma_IetaIeta of the SC, calculated with full 5x5 region """
LowPtElectron_unbiased = NanoAODQuantity("LowPtElectron_unbiased")
"""dtype: Float_t; description: ElectronSeed, pT- and dxy- agnostic BDT (raw) score """

LowPtElectron_miniPFRelIso_all = NanoAODQuantity("LowPtElectron_miniPFRelIso_all")
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections) """
LowPtElectron_miniPFRelIso_chg = NanoAODQuantity("LowPtElectron_miniPFRelIso_chg")
"""dtype: Float_t; description: mini PF relative isolation, charged component """

nMuon = NanoAODQuantity("nMuon")
"""dtype: Int_t; description: slimmedMuons after basic selection (pt > 15 || (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))) """
Muon_IPx = NanoAODQuantity("Muon_IPx")
"""dtype: Float_t; description: x coordinate of impact parameter vector """
Muon_IPy = NanoAODQuantity("Muon_IPy")
"""dtype: Float_t; description: y coordinate of impact parameter vector """
Muon_IPz = NanoAODQuantity("Muon_IPz")
"""dtype: Float_t; description: z coordinate of impact parameter vector """
Muon_bestTrackType = NanoAODQuantity("Muon_bestTrackType")
"""dtype: UChar_t; description: Type of track used (1=inner, 2=STA, 3=global, 4=TPFMS, 5=Picky, 6=DYT) """
Muon_bsConstrainedChi2 = NanoAODQuantity("Muon_bsConstrainedChi2")
"""dtype: Float_t; description: chi2 of beamspot constraint """
Muon_bsConstrainedPt = NanoAODQuantity("Muon_bsConstrainedPt")
"""dtype: Float_t; description: pT with beamspot constraint """
Muon_bsConstrainedPtErr = NanoAODQuantity("Muon_bsConstrainedPtErr")
"""dtype: Float_t; description: pT error with beamspot constraint  """
Muon_charge = NanoAODQuantity("Muon_charge")
"""dtype: Int_t; description: electric charge """
Muon_dxy = NanoAODQuantity("Muon_dxy")
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm """
Muon_dxyErr = NanoAODQuantity("Muon_dxyErr")
"""dtype: Float_t; description: dxy uncertainty, in cm """
Muon_dxybs = NanoAODQuantity("Muon_dxybs")
"""dtype: Float_t; description: dxy (with sign) wrt the beam spot, in cm """
Muon_dxybsErr = NanoAODQuantity("Muon_dxybsErr")
"""dtype: Float_t; description: dxy uncertainty wrt the beam spot, in cm """
Muon_dz = NanoAODQuantity("Muon_dz")
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm """
Muon_dzErr = NanoAODQuantity("Muon_dzErr")
"""dtype: Float_t; description: dz uncertainty, in cm """
Muon_eta = NanoAODQuantity("Muon_eta")
"""dtype: Float_t; description: eta """
Muon_fsrPhotonIdx = NanoAODQuantity("Muon_fsrPhotonIdx")
"""dtype: Short_t; description: Index of the lowest-dR/ET2 among associated FSR photons """
Muon_genPartFlav = NanoAODQuantity("Muon_genPartFlav")
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 muons: 1 = prompt muon (including gamma*->mu mu), 15 = muon from prompt tau, 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched """
Muon_genPartIdx = NanoAODQuantity("Muon_genPartIdx")
"""dtype: Short_t; description: Index into genParticle list for MC matching to status==1 muons """
Muon_highPtId = NanoAODQuantity("Muon_highPtId")
"""dtype: UChar_t; description: high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT) """
Muon_highPurity = NanoAODQuantity("Muon_highPurity")
"""dtype: Bool_t; description: inner track is high purity """
Muon_inTimeMuon = NanoAODQuantity("Muon_inTimeMuon")
"""dtype: Bool_t; description: inTimeMuon ID """
Muon_ip3d = NanoAODQuantity("Muon_ip3d")
"""dtype: Float_t; description: 3D impact parameter wrt first PV, in cm """
Muon_ipLengthSig = NanoAODQuantity("Muon_ipLengthSig")
"""dtype: Float_t; description: significance of impact parameter """
Muon_isGlobal = NanoAODQuantity("Muon_isGlobal")
"""dtype: Bool_t; description: muon is global muon """
Muon_isPFcand = NanoAODQuantity("Muon_isPFcand")
"""dtype: Bool_t; description: muon is PF candidate """
Muon_isStandalone = NanoAODQuantity("Muon_isStandalone")
"""dtype: Bool_t; description: muon is a standalone muon """
Muon_isTracker = NanoAODQuantity("Muon_isTracker")
"""dtype: Bool_t; description: muon is tracker muon """
Muon_jetDF = NanoAODQuantity("Muon_jetDF")
"""dtype: Float_t; description: value of the DEEPJET b tagging algorithm discriminator of the associated jet (0 if none) """
Muon_jetIdx = NanoAODQuantity("Muon_jetIdx")
"""dtype: Short_t; description: index of the associated jet (-1 if none) """
Muon_jetNDauCharged = NanoAODQuantity("Muon_jetNDauCharged")
"""dtype: UChar_t; description: number of charged daughters of the closest jet """
Muon_jetPtRelv2 = NanoAODQuantity("Muon_jetPtRelv2")
"""dtype: Float_t; description: Relative momentum of the lepton with respect to the closest jet after subtracting the lepton """
Muon_jetRelIso = NanoAODQuantity("Muon_jetRelIso")
"""dtype: Float_t; description: Relative isolation in matched jet (1/ptRatio-1), -1 if none """
Muon_looseId = NanoAODQuantity("Muon_looseId")
"""dtype: Bool_t; description: muon is loose muon """
Muon_mass = NanoAODQuantity("Muon_mass")
"""dtype: Float_t; description: mass """
Muon_mediumId = NanoAODQuantity("Muon_mediumId")
"""dtype: Bool_t; description: cut-based ID, medium WP """
Muon_mediumPromptId = NanoAODQuantity("Muon_mediumPromptId")
"""dtype: Bool_t; description: cut-based ID, medium prompt WP """
Muon_miniIsoId = NanoAODQuantity("Muon_miniIsoId")
"""dtype: UChar_t; description: MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight) """
Muon_multiIsoId = NanoAODQuantity("Muon_multiIsoId")
"""dtype: UChar_t; description: MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium) """
Muon_mvaLowPt = NanoAODQuantity("Muon_mvaLowPt")
"""dtype: Float_t; description: Low pt muon ID score """
Muon_mvaMuID = NanoAODQuantity("Muon_mvaMuID")
"""dtype: Float_t; description: MVA-based ID score """
Muon_mvaMuID_WP = NanoAODQuantity("Muon_mvaMuID_WP")
"""dtype: UChar_t; description: MVA-based ID selector WPs (1=MVAIDwpMedium,2=MVAIDwpTight) """
Muon_nStations = NanoAODQuantity("Muon_nStations")
"""dtype: UChar_t; description: number of matched stations with default arbitration (segment & track) """
Muon_nTrackerLayers = NanoAODQuantity("Muon_nTrackerLayers")
"""dtype: UChar_t; description: number of layers in the tracker """
Muon_pdgId = NanoAODQuantity("Muon_pdgId")
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth) """
Muon_pfIsoId = NanoAODQuantity("Muon_pfIsoId")
"""dtype: UChar_t; description: PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight) """
Muon_pfRelIso04_all = NanoAODQuantity("Muon_pfRelIso04_all")
"""dtype: Float_t; description: PF relative isolation dR=0.4, total (deltaBeta corrections) """
Muon_phi = NanoAODQuantity("Muon_phi")
"""dtype: Float_t; description: phi """
Muon_promptMVA = NanoAODQuantity("Muon_promptMVA")
"""dtype: Float_t; description: Prompt MVA lepton ID score. Corresponds to the previous mvaTTH """
Muon_pt = NanoAODQuantity("Muon_pt")
"""dtype: Float_t; description: pt """
Muon_ptErr = NanoAODQuantity("Muon_ptErr")
"""dtype: Float_t; description: ptError of the muon track """
Muon_puppiIsoId = NanoAODQuantity("Muon_puppiIsoId")
"""dtype: UChar_t; description: PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight) """
Muon_segmentComp = NanoAODQuantity("Muon_segmentComp")
"""dtype: Float_t; description: muon segment compatibility """
Muon_sip3d = NanoAODQuantity("Muon_sip3d")
"""dtype: Float_t; description: 3D impact parameter significance wrt first PV """
Muon_softId = NanoAODQuantity("Muon_softId")
"""dtype: Bool_t; description: soft cut-based ID """
Muon_softMva = NanoAODQuantity("Muon_softMva")
"""dtype: Float_t; description: soft MVA ID score """
Muon_softMvaId = NanoAODQuantity("Muon_softMvaId")
"""dtype: Bool_t; description: soft MVA ID """
Muon_softMvaRun3 = NanoAODQuantity("Muon_softMvaRun3")
"""dtype: Float_t; description: soft MVA Run3 ID score """
Muon_svIdx = NanoAODQuantity("Muon_svIdx")
"""dtype: Short_t; description: index of matching secondary vertex """
Muon_tightCharge = NanoAODQuantity("Muon_tightCharge")
"""dtype: UChar_t; description: Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass) """
Muon_tightId = NanoAODQuantity("Muon_tightId")
"""dtype: Bool_t; description: cut-based ID, tight WP """
Muon_tkIsoId = NanoAODQuantity("Muon_tkIsoId")
"""dtype: UChar_t; description: TkIso ID (1=TkIsoLoose, 2=TkIsoTight) """
Muon_tkRelIso = NanoAODQuantity("Muon_tkRelIso")
"""dtype: Float_t; description: Tracker-based relative isolation dR=0.3 for highPt, trkIso/pt """
Muon_triggerIdLoose = NanoAODQuantity("Muon_triggerIdLoose")
"""dtype: Bool_t; description: TriggerIdLoose ID """
Muon_tunepRelPt = NanoAODQuantity("Muon_tunepRelPt")
"""dtype: Float_t; description: TuneP relative pt, tunePpt/pt """

Muon_VXBS_Cov00 = NanoAODQuantity("Muon_VXBS_Cov00")
"""dtype: Float_t; description: 0, 0 element of the VXBS Covariance matrix """
Muon_VXBS_Cov03 = NanoAODQuantity("Muon_VXBS_Cov03")
"""dtype: Float_t; description: 0, 3 element of the VXBS Covariance matrix """
Muon_VXBS_Cov33 = NanoAODQuantity("Muon_VXBS_Cov33")
"""dtype: Float_t; description: 3, 3 element of the VXBS Covariance matrix """

Muon_miniPFRelIso_all = NanoAODQuantity("Muon_miniPFRelIso_all")
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections) """
Muon_miniPFRelIso_chg = NanoAODQuantity("Muon_miniPFRelIso_chg")
"""dtype: Float_t; description: mini PF relative isolation, charged component """

Muon_pfRelIso03_all = NanoAODQuantity("Muon_pfRelIso03_all")
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (deltaBeta corrections) """
Muon_pfRelIso03_chg = NanoAODQuantity("Muon_pfRelIso03_chg")
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component """

Muon_pnScore_heavy = NanoAODQuantity("Muon_pnScore_heavy")
"""dtype: Float_t; description: PNet muon ID score for lepton from B or D hadrons """
Muon_pnScore_light = NanoAODQuantity("Muon_pnScore_light")
"""dtype: Float_t; description: PNet muon ID score for lepton from hadrons w/o b or c quarks OR w/o generator matching """
Muon_pnScore_prompt = NanoAODQuantity("Muon_pnScore_prompt")
"""dtype: Float_t; description: PNet muon ID score for lepton from W/Z/H bosons """
Muon_pnScore_tau = NanoAODQuantity("Muon_pnScore_tau")
"""dtype: Float_t; description: PNet muon ID score for decay of tau to light leptons (mu) """

Muon_tuneP_charge = NanoAODQuantity("Muon_tuneP_charge")
"""dtype: Float_t; description: tunePMuonBestTrack() charge """
Muon_tuneP_pterr = NanoAODQuantity("Muon_tuneP_pterr")
"""dtype: Float_t; description: pTerr from tunePMuonBestTrack """

nOtherPV = NanoAODQuantity("nOtherPV")
"""dtype: Int_t; description:  """
OtherPV_score = NanoAODQuantity("OtherPV_score")
"""dtype: Float_t; description: scores of other primary vertices, excluding the main PV """
OtherPV_z = NanoAODQuantity("OtherPV_z")
"""dtype: Float_t; description: Z position of other primary vertices, excluding the main PV """

nPFCand = NanoAODQuantity("nPFCand")
"""dtype: Int_t; description: PF candidate constituents of AK8 puppi jets (FatJet) with |eta| <= 2.5 """
PFCand_eta = NanoAODQuantity("PFCand_eta")
"""dtype: Float_t; description: eta """
PFCand_mass = NanoAODQuantity("PFCand_mass")
"""dtype: Float_t; description: Puppi-weighted mass """
PFCand_pdgId = NanoAODQuantity("PFCand_pdgId")
"""dtype: Int_t; description: PF candidate type (+/-211 = ChgHad, 130 = NeuHad, 22 = Photon, +/-11 = Electron, +/-13 = Muon, 1 = HFHad, 2 = HFEM) """
PFCand_phi = NanoAODQuantity("PFCand_phi")
"""dtype: Float_t; description: phi """
PFCand_pt = NanoAODQuantity("PFCand_pt")
"""dtype: Float_t; description: Puppi-weighted pt """

PFMET_covXX = NanoAODQuantity("PFMET_covXX")
"""dtype: Float_t; description: xx element of met covariance matrix """
PFMET_covXY = NanoAODQuantity("PFMET_covXY")
"""dtype: Float_t; description: xy element of met covariance matrix """
PFMET_covYY = NanoAODQuantity("PFMET_covYY")
"""dtype: Float_t; description: yy element of met covariance matrix """
PFMET_phi = NanoAODQuantity("PFMET_phi")
"""dtype: Float_t; description: phi """
PFMET_phiUnclusteredDown = NanoAODQuantity("PFMET_phiUnclusteredDown")
"""dtype: Float_t; description: Unclustered down phi """
PFMET_phiUnclusteredUp = NanoAODQuantity("PFMET_phiUnclusteredUp")
"""dtype: Float_t; description: Unclustered up phi """
PFMET_pt = NanoAODQuantity("PFMET_pt")
"""dtype: Float_t; description: pt """
PFMET_ptUnclusteredDown = NanoAODQuantity("PFMET_ptUnclusteredDown")
"""dtype: Float_t; description: Unclustered down pt """
PFMET_ptUnclusteredUp = NanoAODQuantity("PFMET_ptUnclusteredUp")
"""dtype: Float_t; description: Unclustered up pt """
PFMET_significance = NanoAODQuantity("PFMET_significance")
"""dtype: Float_t; description: MET significance """
PFMET_sumEt = NanoAODQuantity("PFMET_sumEt")
"""dtype: Float_t; description: scalar sum of Et """
PFMET_sumPtUnclustered = NanoAODQuantity("PFMET_sumPtUnclustered")
"""dtype: Float_t; description: sumPt used for MET significance """

nPSWeight = NanoAODQuantity("nPSWeight")
"""dtype: Int_t; description:  """
PSWeight = NanoAODQuantity("PSWeight")
"""dtype: Float_t; description: PS weights (w_var / w_nominal) [0] isr.murfac=2.0; [1] fsr.murfac=2.0; [2] isr.murfac=0.5; [3] fsr.murfac=0.5;  """

PV_chi2 = NanoAODQuantity("PV_chi2")
"""dtype: Float_t; description: main primary vertex reduced chi2 """
PV_ndof = NanoAODQuantity("PV_ndof")
"""dtype: Float_t; description: main primary vertex number of degree of freedom """
PV_npvs = NanoAODQuantity("PV_npvs")
"""dtype: UChar_t; description: total number of reconstructed primary vertices """
PV_npvsGood = NanoAODQuantity("PV_npvsGood")
"""dtype: UChar_t; description: number of good reconstructed primary vertices. selection:!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2 """
PV_score = NanoAODQuantity("PV_score")
"""dtype: Float_t; description: main primary vertex score, i.e. sum pt2 of clustered objects """
PV_sumpt2 = NanoAODQuantity("PV_sumpt2")
"""dtype: Float_t; description: sum pt2 of pf charged candidates for the main primary vertex """
PV_sumpx = NanoAODQuantity("PV_sumpx")
"""dtype: Float_t; description: sum px of pf charged candidates for the main primary vertex """
PV_sumpy = NanoAODQuantity("PV_sumpy")
"""dtype: Float_t; description: sum py of pf charged candidates for the main primary vertex """
PV_x = NanoAODQuantity("PV_x")
"""dtype: Float_t; description: main primary vertex position x coordinate """
PV_y = NanoAODQuantity("PV_y")
"""dtype: Float_t; description: main primary vertex position y coordinate """
PV_z = NanoAODQuantity("PV_z")
"""dtype: Float_t; description: main primary vertex position z coordinate """

nPVBS = NanoAODQuantity("nPVBS")
"""dtype: Int_t; description: main primary vertex with beam-spot """
PVBS_chi2 = NanoAODQuantity("PVBS_chi2")
"""dtype: Float_t; description: reduced chi2, i.e. chi2/ndof """
PVBS_cov00 = NanoAODQuantity("PVBS_cov00")
"""dtype: Float_t; description: vertex covariance (0,0) """
PVBS_cov10 = NanoAODQuantity("PVBS_cov10")
"""dtype: Float_t; description: vertex covariance (1,0) """
PVBS_cov11 = NanoAODQuantity("PVBS_cov11")
"""dtype: Float_t; description: vertex covariance (1,1) """
PVBS_cov20 = NanoAODQuantity("PVBS_cov20")
"""dtype: Float_t; description: vertex covariance (2,0) """
PVBS_cov21 = NanoAODQuantity("PVBS_cov21")
"""dtype: Float_t; description: vertex covariance (2,1) """
PVBS_cov22 = NanoAODQuantity("PVBS_cov22")
"""dtype: Float_t; description: vertex covariance (2,2) """
PVBS_x = NanoAODQuantity("PVBS_x")
"""dtype: Float_t; description: position x coordinate, in cm """
PVBS_y = NanoAODQuantity("PVBS_y")
"""dtype: Float_t; description: position y coordinate, in cm """
PVBS_z = NanoAODQuantity("PVBS_z")
"""dtype: Float_t; description: position z coordinate, in cm """

nPhoton = NanoAODQuantity("nPhoton")
"""dtype: Int_t; description: slimmedPhotons after basic selection (pt > 5 ) """
Photon_cutBased = NanoAODQuantity("Photon_cutBased")
"""dtype: UChar_t; description: cut-based ID bitmap, RunIIIWinter22V1: fail ==0, loose >=1 , medium >=2, tight >=3 """
Photon_ecalPFClusterIso = NanoAODQuantity("Photon_ecalPFClusterIso")
"""dtype: Float_t; description: sum pt of ecal clusters, vetoing clusters part of photon """
Photon_electronIdx = NanoAODQuantity("Photon_electronIdx")
"""dtype: Short_t; description: index of the associated electron (-1 if none) """
Photon_electronVeto = NanoAODQuantity("Photon_electronVeto")
"""dtype: Bool_t; description: pass electron veto """
Photon_energyErr = NanoAODQuantity("Photon_energyErr")
"""dtype: Float_t; description: energy error of the cluster from regression """
Photon_energyRaw = NanoAODQuantity("Photon_energyRaw")
"""dtype: Float_t; description: raw energy of photon supercluster """
Photon_esEffSigmaRR = NanoAODQuantity("Photon_esEffSigmaRR")
"""dtype: Float_t; description: preshower sigmaRR """
Photon_esEnergyOverRawE = NanoAODQuantity("Photon_esEnergyOverRawE")
"""dtype: Float_t; description: ratio of preshower energy to raw supercluster energy """
Photon_eta = NanoAODQuantity("Photon_eta")
"""dtype: Float_t; description: eta """
Photon_etaWidth = NanoAODQuantity("Photon_etaWidth")
"""dtype: Float_t; description: Width of the photon supercluster in eta """
Photon_genPartFlav = NanoAODQuantity("Photon_genPartFlav")
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 photons or electrons: 1 = prompt photon, 11 = prompt electron, 0 = unknown or unmatched """
Photon_genPartIdx = NanoAODQuantity("Photon_genPartIdx")
"""dtype: Short_t; description: Index into genParticle list for MC matching to status==1 photons or electrons """
Photon_haloTaggerMVAVal = NanoAODQuantity("Photon_haloTaggerMVAVal")
"""dtype: Float_t; description: Value of MVA based BDT based  beam halo tagger in the Ecal endcap (valid for pT > 200 GeV) """
Photon_hasConversionTracks = NanoAODQuantity("Photon_hasConversionTracks")
"""dtype: Bool_t; description: Variable specifying if photon has associated conversion tracks (one-legged or two-legged) """
Photon_hcalPFClusterIso = NanoAODQuantity("Photon_hcalPFClusterIso")
"""dtype: Float_t; description: sum pt of hcal clusters, vetoing clusters part of photon """
Photon_isScEtaEB = NanoAODQuantity("Photon_isScEtaEB")
"""dtype: Bool_t; description: is supercluster eta within barrel acceptance """
Photon_isScEtaEE = NanoAODQuantity("Photon_isScEtaEE")
"""dtype: Bool_t; description: is supercluster eta within endcap acceptance """
Photon_jetIdx = NanoAODQuantity("Photon_jetIdx")
"""dtype: Short_t; description: index of the associated jet (-1 if none) """
Photon_pfChargedIso = NanoAODQuantity("Photon_pfChargedIso")
"""dtype: Float_t; description: PF absolute isolation dR=0.3, charged component with dxy,dz match to PV """
Photon_pfChargedIsoPFPV = NanoAODQuantity("Photon_pfChargedIsoPFPV")
"""dtype: Float_t; description: PF absolute isolation dR=0.3, charged component (PF PV only) """
Photon_pfChargedIsoWorstVtx = NanoAODQuantity("Photon_pfChargedIsoWorstVtx")
"""dtype: Float_t; description: PF absolute isolation dR=0.3, charged component (Vertex with largest isolation) """
Photon_pfPhoIso03 = NanoAODQuantity("Photon_pfPhoIso03")
"""dtype: Float_t; description: PF absolute isolation dR=0.3, photon component (uncorrected) """
Photon_phi = NanoAODQuantity("Photon_phi")
"""dtype: Float_t; description: phi """
Photon_phiWidth = NanoAODQuantity("Photon_phiWidth")
"""dtype: Float_t; description: Width of the photon supercluster in phi """
Photon_pixelSeed = NanoAODQuantity("Photon_pixelSeed")
"""dtype: Bool_t; description: has pixel seed """
Photon_pt = NanoAODQuantity("Photon_pt")
"""dtype: Float_t; description: pt """
Photon_r9 = NanoAODQuantity("Photon_r9")
"""dtype: Float_t; description: R9 of the supercluster, calculated with full 5x5 region """
Photon_s4 = NanoAODQuantity("Photon_s4")
"""dtype: Float_t; description: e2x2/e5x5 of the supercluster, calculated with full 5x5 region """
Photon_seedGain = NanoAODQuantity("Photon_seedGain")
"""dtype: UChar_t; description: Gain of the seed crystal """
Photon_seediEtaOriX = NanoAODQuantity("Photon_seediEtaOriX")
"""dtype: Short_t; description: iEta or iX of seed crystal. iEta is barrel-only, iX is endcap-only. iEta runs from -85 to +85, with no crystal at iEta=0. iX runs from 1 to 100. """
Photon_seediPhiOriY = NanoAODQuantity("Photon_seediPhiOriY")
"""dtype: Short_t; description: iPhi or iY of seed crystal. iPhi is barrel-only, iY is endcap-only. iPhi runs from 1 to 360. iY runs from 1 to 100. """
Photon_sieie = NanoAODQuantity("Photon_sieie")
"""dtype: Float_t; description: sigma_IetaIeta of the supercluster, calculated with full 5x5 region """
Photon_sieip = NanoAODQuantity("Photon_sieip")
"""dtype: Float_t; description: sigma_IetaIphi of the supercluster, calculated with full 5x5 region """
Photon_sipip = NanoAODQuantity("Photon_sipip")
"""dtype: Float_t; description: sigmaIphiIphi of the supercluster """
Photon_superclusterEta = NanoAODQuantity("Photon_superclusterEta")
"""dtype: Float_t; description: supercluster eta """
Photon_trkSumPtHollowConeDR03 = NanoAODQuantity("Photon_trkSumPtHollowConeDR03")
"""dtype: Float_t; description: Sum of track pT in a hollow cone of outer radius, inner radius """
Photon_trkSumPtSolidConeDR04 = NanoAODQuantity("Photon_trkSumPtSolidConeDR04")
"""dtype: Float_t; description: Sum of track pT in a cone of dR=0.4 """
Photon_vidNestedWPBitmap = NanoAODQuantity("Photon_vidNestedWPBitmap")
"""dtype: Int_t; description: RunIIIWinter22V1 VID compressed bitmap (MinPtCut,PhoSCEtaMultiRangeCut,PhoFull5x5SigmaIEtaIEtaCut,PhoGenericQuadraticRhoPtScaledCut,PhoGenericQuadraticRhoPtScaledCut,PhoGenericQuadraticRhoPtScaledCut,PhoGenericQuadraticRhoPtScaledCut), 2 bits per cut """
Photon_x_calo = NanoAODQuantity("Photon_x_calo")
"""dtype: Float_t; description: photon supercluster position on calorimeter, x coordinate (cm) """
Photon_y_calo = NanoAODQuantity("Photon_y_calo")
"""dtype: Float_t; description: photon supercluster position on calorimeter, y coordinate (cm) """
Photon_z_calo = NanoAODQuantity("Photon_z_calo")
"""dtype: Float_t; description: photon supercluster position on calorimeter, z coordinate (cm) """

Photon_hoe = NanoAODQuantity("Photon_hoe")
"""dtype: Float_t; description: H over E """
Photon_hoe_PUcorr = NanoAODQuantity("Photon_hoe_PUcorr")
"""dtype: Float_t; description: PU corrected H/E (cone-based with quadraticEA*rho*rho + linearEA*rho Winter22V1 corrections) """
Photon_hoe_Tower = NanoAODQuantity("Photon_hoe_Tower")
"""dtype: Float_t; description: H over E Tower based calculation """

Photon_mvaID = NanoAODQuantity("Photon_mvaID")
"""dtype: Float_t; description: MVA ID score, Winter22V1 """
Photon_mvaID_WP80 = NanoAODQuantity("Photon_mvaID_WP80")
"""dtype: Bool_t; description: MVA ID WP80, Winter22V1 """
Photon_mvaID_WP90 = NanoAODQuantity("Photon_mvaID_WP90")
"""dtype: Bool_t; description: MVA ID WP90, Winter22V1 """

Photon_pfRelIso03_all_quadratic = NanoAODQuantity("Photon_pfRelIso03_all_quadratic")
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (with quadraticEA*rho*rho + linearEA*rho Winter22V1 corrections) """
Photon_pfRelIso03_chg_quadratic = NanoAODQuantity("Photon_pfRelIso03_chg_quadratic")
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged hadron component (with quadraticEA*rho*rho + linearEA*rho Winter22V1 corrections) """

Pileup_gpudensity = NanoAODQuantity("Pileup_gpudensity")
"""dtype: Float_t; description: Generator-level PU vertices / mm """
Pileup_nPU = NanoAODQuantity("Pileup_nPU")
"""dtype: Int_t; description: the number of pileup interactions that have been added to the event in the current bunch crossing """
Pileup_nTrueInt = NanoAODQuantity("Pileup_nTrueInt")
"""dtype: Float_t; description: the true mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled """
Pileup_pthatmax = NanoAODQuantity("Pileup_pthatmax")
"""dtype: Float_t; description: Maximum pt-hat """
Pileup_pudensity = NanoAODQuantity("Pileup_pudensity")
"""dtype: Float_t; description: PU vertices / mm """
Pileup_sumEOOT = NanoAODQuantity("Pileup_sumEOOT")
"""dtype: Int_t; description: number of early out of time pileup """
Pileup_sumLOOT = NanoAODQuantity("Pileup_sumLOOT")
"""dtype: Int_t; description: number of late out of time pileup """

PuppiMET_covXX = NanoAODQuantity("PuppiMET_covXX")
"""dtype: Float_t; description: xx element of met covariance matrix """
PuppiMET_covXY = NanoAODQuantity("PuppiMET_covXY")
"""dtype: Float_t; description: xy element of met covariance matrix """
PuppiMET_covYY = NanoAODQuantity("PuppiMET_covYY")
"""dtype: Float_t; description: yy element of met covariance matrix """
PuppiMET_phi = NanoAODQuantity("PuppiMET_phi")
"""dtype: Float_t; description: phi """
PuppiMET_phiUnclusteredDown = NanoAODQuantity("PuppiMET_phiUnclusteredDown")
"""dtype: Float_t; description: Unclustered down phi """
PuppiMET_phiUnclusteredUp = NanoAODQuantity("PuppiMET_phiUnclusteredUp")
"""dtype: Float_t; description: Unclustered up phi """
PuppiMET_pt = NanoAODQuantity("PuppiMET_pt")
"""dtype: Float_t; description: pt """
PuppiMET_ptUnclusteredDown = NanoAODQuantity("PuppiMET_ptUnclusteredDown")
"""dtype: Float_t; description: Unclustered down pt """
PuppiMET_ptUnclusteredUp = NanoAODQuantity("PuppiMET_ptUnclusteredUp")
"""dtype: Float_t; description: Unclustered up pt """
PuppiMET_significance = NanoAODQuantity("PuppiMET_significance")
"""dtype: Float_t; description: MET significance """
PuppiMET_sumEt = NanoAODQuantity("PuppiMET_sumEt")
"""dtype: Float_t; description: scalar sum of Et """
PuppiMET_sumPtUnclustered = NanoAODQuantity("PuppiMET_sumPtUnclustered")
"""dtype: Float_t; description: sumPt used for MET significance """

RawPFMET_phi = NanoAODQuantity("RawPFMET_phi")
"""dtype: Float_t; description: phi """
RawPFMET_pt = NanoAODQuantity("RawPFMET_pt")
"""dtype: Float_t; description: pt """
RawPFMET_sumEt = NanoAODQuantity("RawPFMET_sumEt")
"""dtype: Float_t; description: scalar sum of Et """

RawPuppiMET_phi = NanoAODQuantity("RawPuppiMET_phi")
"""dtype: Float_t; description: phi """
RawPuppiMET_pt = NanoAODQuantity("RawPuppiMET_pt")
"""dtype: Float_t; description: pt """
RawPuppiMET_sumEt = NanoAODQuantity("RawPuppiMET_sumEt")
"""dtype: Float_t; description: scalar sum of Et """

Rho_fixedGridRhoAll = NanoAODQuantity("Rho_fixedGridRhoAll")
"""dtype: Float_t; description: rho from all PF Candidates, no foreground removal (for isolation of prompt photons) """
Rho_fixedGridRhoFastjetAll = NanoAODQuantity("Rho_fixedGridRhoFastjetAll")
"""dtype: Float_t; description: rho from all PF Candidates, used e.g. for JECs """
Rho_fixedGridRhoFastjetCentral = NanoAODQuantity("Rho_fixedGridRhoFastjetCentral")
"""dtype: Float_t; description: rho from all PF Candidates for central region, used e.g. for JECs """
Rho_fixedGridRhoFastjetCentralCalo = NanoAODQuantity(
    "Rho_fixedGridRhoFastjetCentralCalo"
)
"""dtype: Float_t; description: rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation """
Rho_fixedGridRhoFastjetCentralChargedPileUp = NanoAODQuantity(
    "Rho_fixedGridRhoFastjetCentralChargedPileUp"
)
"""dtype: Float_t; description: rho from charged PF Candidates for central region, used e.g. for JECs """
Rho_fixedGridRhoFastjetCentralNeutral = NanoAODQuantity(
    "Rho_fixedGridRhoFastjetCentralNeutral"
)
"""dtype: Float_t; description: rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations """

nSV = NanoAODQuantity("nSV")
"""dtype: Int_t; description: secondary vertices from IVF algorithm """
SV_charge = NanoAODQuantity("SV_charge")
"""dtype: Short_t; description: sum of the charge of the SV tracks """
SV_chi2 = NanoAODQuantity("SV_chi2")
"""dtype: Float_t; description: reduced chi2, i.e. chi/ndof """
SV_dlen = NanoAODQuantity("SV_dlen")
"""dtype: Float_t; description: decay length in cm """
SV_dlenSig = NanoAODQuantity("SV_dlenSig")
"""dtype: Float_t; description: decay length significance """
SV_dxy = NanoAODQuantity("SV_dxy")
"""dtype: Float_t; description: 2D decay length in cm """
SV_dxySig = NanoAODQuantity("SV_dxySig")
"""dtype: Float_t; description: 2D decay length significance """
SV_eta = NanoAODQuantity("SV_eta")
"""dtype: Float_t; description: eta """
SV_mass = NanoAODQuantity("SV_mass")
"""dtype: Float_t; description: mass """
SV_ndof = NanoAODQuantity("SV_ndof")
"""dtype: Float_t; description: number of degrees of freedom """
SV_ntracks = NanoAODQuantity("SV_ntracks")
"""dtype: UChar_t; description: number of tracks """
SV_pAngle = NanoAODQuantity("SV_pAngle")
"""dtype: Float_t; description: pointing angle, i.e. acos(p_SV * (SV - PV))  """
SV_phi = NanoAODQuantity("SV_phi")
"""dtype: Float_t; description: phi """
SV_pt = NanoAODQuantity("SV_pt")
"""dtype: Float_t; description: pt """
SV_x = NanoAODQuantity("SV_x")
"""dtype: Float_t; description: secondary vertex X position, in cm """
SV_y = NanoAODQuantity("SV_y")
"""dtype: Float_t; description: secondary vertex Y position, in cm """
SV_z = NanoAODQuantity("SV_z")
"""dtype: Float_t; description: secondary vertex Z position, in cm """

nSoftActivityJet = NanoAODQuantity("nSoftActivityJet")
"""dtype: Int_t; description: jets clustered from charged candidates compatible with primary vertex (charge()!=0 && pvAssociationQuality()>=5 && vertexRef().key()==0) """
SoftActivityJet_eta = NanoAODQuantity("SoftActivityJet_eta")
"""dtype: Float_t; description: eta """
SoftActivityJet_phi = NanoAODQuantity("SoftActivityJet_phi")
"""dtype: Float_t; description: phi """
SoftActivityJet_pt = NanoAODQuantity("SoftActivityJet_pt")
"""dtype: Float_t; description: pt """

nSubGenJetAK8 = NanoAODQuantity("nSubGenJetAK8")
"""dtype: Int_t; description: slimmedGenJetsAK8SoftDropSubJets, i.e. subjets of ak8 Jets made with visible genparticles """
SubGenJetAK8_eta = NanoAODQuantity("SubGenJetAK8_eta")
"""dtype: Float_t; description: eta """
SubGenJetAK8_mass = NanoAODQuantity("SubGenJetAK8_mass")
"""dtype: Float_t; description: mass """
SubGenJetAK8_phi = NanoAODQuantity("SubGenJetAK8_phi")
"""dtype: Float_t; description: phi """
SubGenJetAK8_pt = NanoAODQuantity("SubGenJetAK8_pt")
"""dtype: Float_t; description: pt """

nSubJet = NanoAODQuantity("nSubJet")
"""dtype: Int_t; description: slimmedJetsAK8PFPuppiSoftDropPacked::SubJets, i.e. soft-drop subjets for ak8 fat jets for boosted analysis """
SubJet_UParTAK4RegPtRawCorr = NanoAODQuantity("SubJet_UParTAK4RegPtRawCorr")
"""dtype: Float_t; description: UnifiedParT universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT """
SubJet_UParTAK4RegPtRawCorrNeutrino = NanoAODQuantity(
    "SubJet_UParTAK4RegPtRawCorrNeutrino"
)
"""dtype: Float_t; description: UnifiedParT universal flavor-aware pT regression neutrino correction, relative to visible. Correction relative to raw jet pT """
SubJet_UParTAK4RegPtRawRes = NanoAODQuantity("SubJet_UParTAK4RegPtRawRes")
"""dtype: Float_t; description: UnifiedParT universal flavor-aware jet pT resolution estimator, (q84 - q16)/2 """
SubJet_UParTAK4V1RegPtRawCorr = NanoAODQuantity("SubJet_UParTAK4V1RegPtRawCorr")
"""dtype: Float_t; description: UnifiedParT V1 universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT """
SubJet_UParTAK4V1RegPtRawCorrNeutrino = NanoAODQuantity(
    "SubJet_UParTAK4V1RegPtRawCorrNeutrino"
)
"""dtype: Float_t; description: UnifiedParT V1 universal flavor-aware pT regression neutrino correction, relative to visible. Correction relative to raw jet pT """
SubJet_UParTAK4V1RegPtRawRes = NanoAODQuantity("SubJet_UParTAK4V1RegPtRawRes")
"""dtype: Float_t; description: UnifiedParT V1 universal flavor-aware jet pT resolution estimator, (q84 - q16)/2 """
SubJet_area = NanoAODQuantity("SubJet_area")
"""dtype: Float_t; description: jet catchment area, for JECs """
SubJet_btagDeepFlavB = NanoAODQuantity("SubJet_btagDeepFlavB")
"""dtype: Float_t; description: DeepJet b+bb+lepb tag discriminator """
SubJet_btagUParTAK4B = NanoAODQuantity("SubJet_btagUParTAK4B")
"""dtype: Float_t; description: UnifiedParT b vs. udscg """
SubJet_eta = NanoAODQuantity("SubJet_eta")
"""dtype: Float_t; description: eta """
SubJet_hadronFlavour = NanoAODQuantity("SubJet_hadronFlavour")
"""dtype: UChar_t; description: flavour from hadron ghost clustering """
SubJet_mass = NanoAODQuantity("SubJet_mass")
"""dtype: Float_t; description: mass """
SubJet_n2b1 = NanoAODQuantity("SubJet_n2b1")
"""dtype: Float_t; description: N2 with beta=1 """
SubJet_n3b1 = NanoAODQuantity("SubJet_n3b1")
"""dtype: Float_t; description: N3 with beta=1 """
SubJet_nBHadrons = NanoAODQuantity("SubJet_nBHadrons")
"""dtype: UChar_t; description: number of b-hadrons """
SubJet_nCHadrons = NanoAODQuantity("SubJet_nCHadrons")
"""dtype: UChar_t; description: number of c-hadrons """
SubJet_phi = NanoAODQuantity("SubJet_phi")
"""dtype: Float_t; description: phi """
SubJet_pt = NanoAODQuantity("SubJet_pt")
"""dtype: Float_t; description: pt """
SubJet_rawFactor = NanoAODQuantity("SubJet_rawFactor")
"""dtype: Float_t; description: 1 - Factor to get back to raw pT """
SubJet_subGenJetAK8Idx = NanoAODQuantity("SubJet_subGenJetAK8Idx")
"""dtype: Short_t; description: index of matched gen subjet in SubGenJetAK8 """
SubJet_tau1 = NanoAODQuantity("SubJet_tau1")
"""dtype: Float_t; description: Nsubjettiness (1 axis) """
SubJet_tau2 = NanoAODQuantity("SubJet_tau2")
"""dtype: Float_t; description: Nsubjettiness (2 axis) """
SubJet_tau3 = NanoAODQuantity("SubJet_tau3")
"""dtype: Float_t; description: Nsubjettiness (3 axis) """
SubJet_tau4 = NanoAODQuantity("SubJet_tau4")
"""dtype: Float_t; description: Nsubjettiness (4 axis) """

nTau = NanoAODQuantity("nTau")
"""dtype: Int_t; description: slimmedTaus after basic selection (pt > 18 && ((tauID('decayModeFindingNewDMs') > 0.5 && (tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') || (tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum'))<2.5) || tauID('byVVVLooseDeepTau2018v2p5VSjet'))) || (?isTauIDAvailable('byUTagCHSVSjetraw')?tauID('byUTagCHSVSjetraw'):-1) > 0.05 || (?isTauIDAvailable('byUTagPUPPIVSjetraw')?tauID('byUTagPUPPIVSjetraw'):-1) > 0.05)) """
Tau_IPx = NanoAODQuantity("Tau_IPx")
"""dtype: Float_t; description: x coordinate of impact parameter vector """
Tau_IPy = NanoAODQuantity("Tau_IPy")
"""dtype: Float_t; description: y coordinate of impact parameter vector """
Tau_IPz = NanoAODQuantity("Tau_IPz")
"""dtype: Float_t; description: z coordinate of impact parameter vector """
Tau_charge = NanoAODQuantity("Tau_charge")
"""dtype: Short_t; description: electric charge """
Tau_chargedIso = NanoAODQuantity("Tau_chargedIso")
"""dtype: Float_t; description: charged isolation """
Tau_decayMode = NanoAODQuantity("Tau_decayMode")
"""dtype: UChar_t; description: decayMode() """
Tau_decayModePNet = NanoAODQuantity("Tau_decayModePNet")
"""dtype: Short_t; description: decay mode of the highest tau score of ParticleNet (CHS Jets) """
Tau_decayModeUParT = NanoAODQuantity("Tau_decayModeUParT")
"""dtype: Short_t; description: decay mode of the highest tau score of Unified ParT 2024 (PUPPI Jets) """
Tau_dxy = NanoAODQuantity("Tau_dxy")
"""dtype: Float_t; description: d_{xy} of lead track with respect to PV, in cm (with sign) """
Tau_dz = NanoAODQuantity("Tau_dz")
"""dtype: Float_t; description: d_{z} of lead track with respect to PV, in cm (with sign) """
Tau_eleIdx = NanoAODQuantity("Tau_eleIdx")
"""dtype: Short_t; description: index of first matching electron """
Tau_eta = NanoAODQuantity("Tau_eta")
"""dtype: Float_t; description: eta """
Tau_genPartFlav = NanoAODQuantity("Tau_genPartFlav")
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==2 taus: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched """
Tau_genPartIdx = NanoAODQuantity("Tau_genPartIdx")
"""dtype: Short_t; description: Index into genParticle list for MC matching to status==2 taus """
Tau_hasRefitSV = NanoAODQuantity("Tau_hasRefitSV")
"""dtype: Bool_t; description: has SV refit using miniAOD quantities """
Tau_idAntiEleDeadECal = NanoAODQuantity("Tau_idAntiEleDeadECal")
"""dtype: Bool_t; description: Anti-electron dead-ECal discriminator """
Tau_idAntiMu = NanoAODQuantity("Tau_idAntiMu")
"""dtype: UChar_t; description: Anti-muon discriminator V3: : 1 = Loose, 2 = Tight """
Tau_idDecayModeNewDMs = NanoAODQuantity("Tau_idDecayModeNewDMs")
"""dtype: Bool_t; description: (?isTauIDAvailable('decayModeFindingNewDMs')?tauID('decayModeFindingNewDMs'):-1) > 0 """
Tau_idDecayModeOldDMs = NanoAODQuantity("Tau_idDecayModeOldDMs")
"""dtype: Bool_t; description: (?isTauIDAvailable('decayModeFinding')?tauID('decayModeFinding'):-1) > 0 """
Tau_idDeepTau2018v2p5VSe = NanoAODQuantity("Tau_idDeepTau2018v2p5VSe")
"""dtype: UChar_t; description: byDeepTau2018v2p5VSe ID working points (deepTau2018v2p5): 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight """
Tau_idDeepTau2018v2p5VSjet = NanoAODQuantity("Tau_idDeepTau2018v2p5VSjet")
"""dtype: UChar_t; description: byDeepTau2018v2p5VSjet ID working points (deepTau2018v2p5): 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight """
Tau_idDeepTau2018v2p5VSmu = NanoAODQuantity("Tau_idDeepTau2018v2p5VSmu")
"""dtype: UChar_t; description: byDeepTau2018v2p5VSmu ID working points (deepTau2018v2p5): 1 = VLoose, 2 = Loose, 3 = Medium, 4 = Tight """
Tau_ipLengthSig = NanoAODQuantity("Tau_ipLengthSig")
"""dtype: Float_t; description: significance of impact parameter """
Tau_jetIdx = NanoAODQuantity("Tau_jetIdx")
"""dtype: Short_t; description: index of the associated jet (-1 if none) """
Tau_leadTkDeltaEta = NanoAODQuantity("Tau_leadTkDeltaEta")
"""dtype: Float_t; description: eta of the leading track, minus tau eta """
Tau_leadTkDeltaPhi = NanoAODQuantity("Tau_leadTkDeltaPhi")
"""dtype: Float_t; description: phi of the leading track, minus tau phi """
Tau_leadTkPtOverTauPt = NanoAODQuantity("Tau_leadTkPtOverTauPt")
"""dtype: Float_t; description: pt of the leading track divided by tau pt """
Tau_mass = NanoAODQuantity("Tau_mass")
"""dtype: Float_t; description: mass """
Tau_muIdx = NanoAODQuantity("Tau_muIdx")
"""dtype: Short_t; description: index of first matching muon """
Tau_nSVs = NanoAODQuantity("Tau_nSVs")
"""dtype: UChar_t; description: number of secondary vertices in the tau """
Tau_neutralIso = NanoAODQuantity("Tau_neutralIso")
"""dtype: Float_t; description: neutral (photon) isolation """
Tau_phi = NanoAODQuantity("Tau_phi")
"""dtype: Float_t; description: phi """
Tau_photonsOutsideSignalCone = NanoAODQuantity("Tau_photonsOutsideSignalCone")
"""dtype: Float_t; description: sum of photons outside signal cone """
Tau_probDM0PNet = NanoAODQuantity("Tau_probDM0PNet")
"""dtype: Float_t; description: normalised probablity of decayMode 0, 1h+0pi0 (PNet 2023 - CHS Jets) """
Tau_probDM0UParT = NanoAODQuantity("Tau_probDM0UParT")
"""dtype: Float_t; description: normalised probablity of decayMode 0, 1h+0pi0 (Unified ParT 2024 - PUPPI Jets) """
Tau_probDM10PNet = NanoAODQuantity("Tau_probDM10PNet")
"""dtype: Float_t; description: normalised probablity of decayMode 10, 3h+0pi0 (PNet 2023 - CHS Jets) """
Tau_probDM10UParT = NanoAODQuantity("Tau_probDM10UParT")
"""dtype: Float_t; description: normalised probablity of decayMode 10, 3h+0pi0 (Unified ParT 2024 - PUPPI Jets) """
Tau_probDM11PNet = NanoAODQuantity("Tau_probDM11PNet")
"""dtype: Float_t; description: normalised probablity of decayMode 11, 3h+1pi0 (PNet 2023 - CHS Jets) """
Tau_probDM11UParT = NanoAODQuantity("Tau_probDM11UParT")
"""dtype: Float_t; description: normalised probablity of decayMode 11, 3h+1pi0 (Unified ParT 2024 - PUPPI Jets) """
Tau_probDM1PNet = NanoAODQuantity("Tau_probDM1PNet")
"""dtype: Float_t; description: normalised probablity of decayMode 1, 1h+1pi0 (PNet 2023 - CHS Jets) """
Tau_probDM1UParT = NanoAODQuantity("Tau_probDM1UParT")
"""dtype: Float_t; description: normalised probablity of decayMode 1, 1h+1pi0 (Unified ParT 2024 - PUPPI Jets) """
Tau_probDM2PNet = NanoAODQuantity("Tau_probDM2PNet")
"""dtype: Float_t; description: normalised probablity of decayMode 2, 1h+2pi0 (PNet 2023 - CHS Jets) """
Tau_probDM2UParT = NanoAODQuantity("Tau_probDM2UParT")
"""dtype: Float_t; description: normalised probablity of decayMode 2, 1h+2pi0 (Unified ParT 2024 - PUPPI Jets) """
Tau_pt = NanoAODQuantity("Tau_pt")
"""dtype: Float_t; description: pt """
Tau_ptCorrPNet = NanoAODQuantity("Tau_ptCorrPNet")
"""dtype: Float_t; description: pt correction (PNet 2023 - CHS Jets) """
Tau_ptCorrUParT = NanoAODQuantity("Tau_ptCorrUParT")
"""dtype: Float_t; description: pt correction (Unified ParT 2024 - PUPPI Jets) """
Tau_puCorr = NanoAODQuantity("Tau_puCorr")
"""dtype: Float_t; description: pileup correction """
Tau_qConfPNet = NanoAODQuantity("Tau_qConfPNet")
"""dtype: Float_t; description: signed charge confidence (PNet 2023 - CHS Jets) """
Tau_qConfUParT = NanoAODQuantity("Tau_qConfUParT")
"""dtype: Float_t; description: signed charge confidence (Unified ParT 2024 - PUPPI Jets) """
Tau_rawDeepTau2018v2p5VSe = NanoAODQuantity("Tau_rawDeepTau2018v2p5VSe")
"""dtype: Float_t; description: byDeepTau2018v2p5VSe raw output discriminator (deepTau2018v2p5) """
Tau_rawDeepTau2018v2p5VSjet = NanoAODQuantity("Tau_rawDeepTau2018v2p5VSjet")
"""dtype: Float_t; description: byDeepTau2018v2p5VSjet raw output discriminator (deepTau2018v2p5) """
Tau_rawDeepTau2018v2p5VSmu = NanoAODQuantity("Tau_rawDeepTau2018v2p5VSmu")
"""dtype: Float_t; description: byDeepTau2018v2p5VSmu raw output discriminator (deepTau2018v2p5) """
Tau_rawIso = NanoAODQuantity("Tau_rawIso")
"""dtype: Float_t; description: combined isolation (deltaBeta corrections) """
Tau_rawIsodR03 = NanoAODQuantity("Tau_rawIsodR03")
"""dtype: Float_t; description: combined isolation (deltaBeta corrections, dR=0.3) """
Tau_rawPNetVSe = NanoAODQuantity("Tau_rawPNetVSe")
"""dtype: Float_t; description: raw output of ParticleNetVsE discriminator (PNet 2023 - CHS Jets) """
Tau_rawPNetVSjet = NanoAODQuantity("Tau_rawPNetVSjet")
"""dtype: Float_t; description: raw output of ParticleNetVsJet discriminator (PNet 2023 - CHS Jets) """
Tau_rawPNetVSmu = NanoAODQuantity("Tau_rawPNetVSmu")
"""dtype: Float_t; description: raw output of ParticleNetVsMu discriminator (PNet 2023 - CHS Jets) """
Tau_rawUParTVSe = NanoAODQuantity("Tau_rawUParTVSe")
"""dtype: Float_t; description: raw output of UParTVsE discriminator (Unified ParT 2024 - PUPPI Jets) """
Tau_rawUParTVSjet = NanoAODQuantity("Tau_rawUParTVSjet")
"""dtype: Float_t; description: raw output of UParTVsJet discriminator (Unified ParT 2024 - PUPPI Jets) """
Tau_rawUParTVSmu = NanoAODQuantity("Tau_rawUParTVSmu")
"""dtype: Float_t; description: raw output of UParTVsMu discriminator (Unified ParT 2024 - PUPPI Jets) """
Tau_refitSVchi2 = NanoAODQuantity("Tau_refitSVchi2")
"""dtype: Float_t; description: reduced chi2, i.e. chi2/ndof, of SV fit """
Tau_refitSVcov00 = NanoAODQuantity("Tau_refitSVcov00")
"""dtype: Float_t; description: Covariance of SV (0,0) """
Tau_refitSVcov10 = NanoAODQuantity("Tau_refitSVcov10")
"""dtype: Float_t; description: Covariance of SV (1,0) """
Tau_refitSVcov11 = NanoAODQuantity("Tau_refitSVcov11")
"""dtype: Float_t; description: Covariance of SV (1,1) """
Tau_refitSVcov20 = NanoAODQuantity("Tau_refitSVcov20")
"""dtype: Float_t; description: Covariance of SV (2,0) """
Tau_refitSVcov21 = NanoAODQuantity("Tau_refitSVcov21")
"""dtype: Float_t; description: Covariance of SV (2,1) """
Tau_refitSVcov22 = NanoAODQuantity("Tau_refitSVcov22")
"""dtype: Float_t; description: Covariance of SV (2,2) """
Tau_refitSVx = NanoAODQuantity("Tau_refitSVx")
"""dtype: Float_t; description: x coordinate of SV """
Tau_refitSVy = NanoAODQuantity("Tau_refitSVy")
"""dtype: Float_t; description: y coordinate of SV """
Tau_refitSVz = NanoAODQuantity("Tau_refitSVz")
"""dtype: Float_t; description: z coordinate of SV """
Tau_svIdx1 = NanoAODQuantity("Tau_svIdx1")
"""dtype: Short_t; description: index of first matching secondary vertex """
Tau_svIdx2 = NanoAODQuantity("Tau_svIdx2")
"""dtype: Short_t; description: index of second matching secondary vertex """

nTauProd = NanoAODQuantity("nTauProd")
"""dtype: Int_t; description: tau signal candidates """
TauProd_eta = NanoAODQuantity("TauProd_eta")
"""dtype: Float_t; description: eta """
TauProd_pdgId = NanoAODQuantity("TauProd_pdgId")
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth) """
TauProd_phi = NanoAODQuantity("TauProd_phi")
"""dtype: Float_t; description: phi """
TauProd_pt = NanoAODQuantity("TauProd_pt")
"""dtype: Float_t; description: pt """
TauProd_tauIdx = NanoAODQuantity("TauProd_tauIdx")
"""dtype: Short_t; description: index of the mother tau """

TauSpinner_weight_cp_0 = NanoAODQuantity("TauSpinner_weight_cp_0")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0 """
TauSpinner_weight_cp_0_alt = NanoAODQuantity("TauSpinner_weight_cp_0_alt")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0 (alternative hadronic currents) """
TauSpinner_weight_cp_0p25 = NanoAODQuantity("TauSpinner_weight_cp_0p25")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0p25 """
TauSpinner_weight_cp_0p25_alt = NanoAODQuantity("TauSpinner_weight_cp_0p25_alt")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0p25 (alternative hadronic currents) """
TauSpinner_weight_cp_0p375 = NanoAODQuantity("TauSpinner_weight_cp_0p375")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0p375 """
TauSpinner_weight_cp_0p375_alt = NanoAODQuantity("TauSpinner_weight_cp_0p375_alt")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0p375 (alternative hadronic currents) """
TauSpinner_weight_cp_0p5 = NanoAODQuantity("TauSpinner_weight_cp_0p5")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0p5 """
TauSpinner_weight_cp_0p5_alt = NanoAODQuantity("TauSpinner_weight_cp_0p5_alt")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = 0p5 (alternative hadronic currents) """
TauSpinner_weight_cp_minus0p25 = NanoAODQuantity("TauSpinner_weight_cp_minus0p25")
"""dtype: Double_t; description: TauSpinner weight for theta_CP = minus0p25 """
TauSpinner_weight_cp_minus0p25_alt = NanoAODQuantity(
    "TauSpinner_weight_cp_minus0p25_alt"
)
"""dtype: Double_t; description: TauSpinner weight for theta_CP = minus0p25 (alternative hadronic currents) """

nTrackGenJetAK4 = NanoAODQuantity("nTrackGenJetAK4")
"""dtype: Int_t; description: AK4 GenJets made with charged particles only """
TrackGenJetAK4_eta = NanoAODQuantity("TrackGenJetAK4_eta")
"""dtype: Float_t; description: eta """
TrackGenJetAK4_phi = NanoAODQuantity("TrackGenJetAK4_phi")
"""dtype: Float_t; description: phi """
TrackGenJetAK4_pt = NanoAODQuantity("TrackGenJetAK4_pt")
"""dtype: Float_t; description: pt """

nTrigObj = NanoAODQuantity("nTrigObj")
"""dtype: Int_t; description:  """
TrigObj_eta = NanoAODQuantity("TrigObj_eta")
"""dtype: Float_t; description: eta """
TrigObj_filterBits = NanoAODQuantity("TrigObj_filterBits")
"""dtype: ULong64_t; description: extra bits of associated information: 0 => HLT_AK8PFJetX_SoftDropMassY_PFAK8ParticleNetTauTau0p30, 1 => HLT_AK8PFJetX_SoftDropMassY_PNetTauTau0p03, 2 => HLT_AK8PFJetX_SoftDropMassY_PNetTauTau0p05 for BoostedTau; 0 => CaloIdL_TrackIdL_IsoVL, 1 => 1e (WPTight with possibile contribution from Xtriggers besides singleElectron), 2 => 1e (WPLoose), 3 => OverlapFilter PFTau, 4 => 2e (Leg 1), 5 => 2e (Leg 2), 6 => 1e-1mu, 7 => 1e-1tau, 8 => 3e, 9 => 2e-1mu, 10 => 1e-2mu, 11 => 1e (32_L1DoubleEG_AND_L1SingleEGOr), 12 => 1e (CaloIdVT_GsfTrkIdT), 13 => 1e (PFJet), 14 => 1e (Photon175_OR_Photon200), 15 => 2e (CaloIdL_MW seeded), 16 => 2e (CaloIdL_MW unseeded), 17 => 1e-1tau PNet, 18 => 1e (HLT30WPTightGSfTrackIso), 19 => WPTightGsfTrackIso for VBF for Electron; 0 => , 1 => hltAK8SingleCaloJet200, 2 => , 3 => hltAK8SinglePFJets230SoftDropMass40BTagParticleNetBB0p35 OR hltAK8SinglePFJets250SoftDropMass40BTagParticleNetBB0p35 OR hltAK8SinglePFJets275SoftDropMass40BTagParticleNetBB0p35, 4 => hltAK8DoublePFJetSDModMass30, 5 => hltAK8DoublePFJetSDModMass50, 6 => hltAK8SinglePFJets*SoftDropMass*PNetBBTag0p06 for FatJet; 0 => hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet, 1 => hltL1sQuadJetC50IorQuadJetC60IorHTT280IorHTT300IorHTT320IorTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBF OR hltL1sQuadJetCIorTripleJetVBFIorHTT, 2 => hltL1sQuadJetC60IorHTT380IorHTT280QuadJetIorHTT300QuadJet OR hltL1sQuadJetC50to60IorHTT280to500IorHTT250to340QuadJet, 3 => hltCaloQuadJet30HT300 OR hltCaloQuadJet30HT320, 4 => hltPFCentralJetsLooseIDQuad30HT300 OR hltPFCentralJetsLooseIDQuad30HT330, 5 => hltPFHT280Jet30 for HT; 0 => hlt4PixelOnlyPFCentralJetTightIDPt20, 1 => hlt3PixelOnlyPFCentralJetTightIDPt30, 2 => hltPFJetFilterTwoC30, 3 => hlt4PFCentralJetTightIDPt30, 4 => hlt4PFCentralJetTightIDPt35, 5 => hltQuadCentralJet30, 6 => hlt2PixelOnlyPFCentralJetTightIDPt40, 7 => hltL1sTripleJet1008572VBFIorHTTIorDoubleJetCIorSingleJet OR hltL1sTripleJet1058576VBFIorHTTIorDoubleJetCIorSingleJet OR hltL1sTripleJetVBFIorHTTIorSingleJet, 8 => hlt3PFCentralJetTightIDPt40, 9 => hlt3PFCentralJetTightIDPt45, 10 => hltL1sQuadJetC60IorHTT380IorHTT280QuadJetIorHTT300QuadJet OR hltL1sQuadJetC50to60IorHTT280to500IorHTT250to340QuadJet, 11 => hltBTagCaloDeepCSVp17Double, 12 => hltPFCentralJetLooseIDQuad30, 13 => hlt1PFCentralJetLooseID75, 14 => hlt2PFCentralJetLooseID60, 15 => hlt3PFCentralJetLooseID45, 16 => hlt4PFCentralJetLooseID40, 17 => (DiTau+Jet (Jet) Signal), 18 => (VBF DiTau Jets), 19 => (VBF SingleTau Jets), 20 => Muon+Tau+Jet (Jet) Monitoring, 21 => hlt2PFCentralJetTightIDPt50, 22 => hlt1PixelOnlyPFCentralJetTightIDPt60, 23 => hlt1PFCentralJetTightIDPt70, 24 => hltBTagPFDeepJet1p5Single, 25 => hltBTagPFDeepJet4p5Triple, 26 => hltBTagCentralJetPt35PFParticleNet2BTagSum0p65 OR hltBTagCentralJetPt30PFParticleNet2BTagSum0p65 OR hltPFJetTwoC30PFBTagParticleNet2BTagSum0p65 OR hltPFCentralJetPt30PNet2BTagMean0p55, 27 => hlt2PixelOnlyPFCentralJetTightIDPt20 OR hlt1PixelOnlyPFCentralJetTightIDPt50, 28 => hlt2PFCentralJetTightIDPt30 OR hltPF2CentralJetTightIDPt30, 29 => hlt1PFCentralJetTightIDPt60, 30 => hltPF2CentralJetPt30PNet2BTagMean0p50, 31 => hlt4PFCentralJetPt25, 32 => hltPFCentralJetNoIDPt25PNet1BTag0p20, 33 => hltPFCentralJetNoIDPt25PNet1TauHTag0p50, 34 => hlt4PFCentralJetTightIDPt25, 35 => hltPFCentralJetPt25PNet2BTagMean0p55, 36 => hltL1PFJetCategoriesVBFinclLoose OR hltL1PFJetCategoriesVBFinclLooseTripleJet OR hltL1PFJetCategoriesVBFinclTight1050, 37 => hltL1PFJetCategoriesVBFdijetQuadjet OR hltL1PFJetCategoriesVBFdijetFivejets OR hltL1PFJetCategoriesVBFdijetSixjets OR hltL1PFJetCategoriesVBFdijetTightQuadjet800, 38 => hltL1PFJetCategoriesVBFMET OR hltL1PFJetCategoriesVBFMETTripleJet OR hltL1PFJetCategoriesVBFMETTight650, 39 => hltL1PFJetCategoriesVBFMu OR hltL1PFJetCategoriesVBFMuTripleJet OR hltL1PFJetCategoriesVBFMuTight750, 40 => hltOverlapFilterDoublePFJet45Photon12 OR hltOverlapFilterDoublePFJet45Photon17 OR hltDiPFJet50Photon22OverlapFilter, 41 => hltOverlapFilterDoublePFJet45Ele12 OR hltOverlapFilterDoublePFJet45Ele17 OR hltDiPFJet50Ele22OverlapFilter, 42 => SinglePFJetX, 43 => SinglePFJetFwdX, 44 => DiPFJetAveX, 45 => DiPFJetAveX_HFJEC for Jet; 0 => hltCaloQuadJet30HT300 OR hltCaloQuadJet30HT320, 1 => hltPFCentralJetsLooseIDQuad30HT300 OR hltPFCentralJetsLooseIDQuad30HT330 for MHT; 0 => TrkIsoVVL, 1 => Iso, 2 => OverlapFilter PFTau, 3 => 1mu, 4 => 2mu, 5 => 1mu-1e, 6 => 1mu-1tau, 7 => 3mu, 8 => 2mu-1e, 9 => 1mu-2e, 10 => 1mu (Mu50), 11 => 1mu (Mu100), 12 => 1mu-1photon, 13 => 1mu-1tau PNet for Muon; 0 => hltEG33L1EG26HEFilter, 1 => hltEG50HEFilter, 2 => hltEG75HEFilter, 3 => hltEG90HEFilter, 4 => hltEG120HEFilter, 5 => hltEG150HEFilter, 6 => hltEG175HEFilter, 7 => hltEG200HEFilter, 8 => hltHtEcal800, 9 => hltEG45EBTightIDTightIsoTrackIsoFilter, 10 => hltEG50EBTightIDTightIsoTrackIsoFilter, 11 => hltEG110EBTightIDTightIsoTrackIsoFilter, 12 => hltEG120EBTightIDTightIsoTrackIsoFilter, 13 => 1mu-1photon, 14 => hltEG30LR9Id85b90eHE12R9Id50b80eR9IdLastFilter, 15 => hltEG30LIso60CaloId15b35eHE12R9Id50b80eEcalIsoLastFilter, 16 => hltEG22Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter, 17 => hltEG22R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter, 18 => hltEG30Iso60CaloId15b35eR9Id50b90eHE12b10eR9Id50b80eEcalIsoFilter, 19 => hltEG18TrackIso60Iso60CaloId15b35eR9Id50b90eHE12b10eR9Id50b80eTrackIsoUnseededFilter, 20 => hltEG*L1VBFLooseIsoEGHEFilter for Photon; 0 => Loose, 1 => Medium, 2 => Tight, 3 => DeepTau no spec WP, 4 => PNet no specified WP, 5 => ChargedIso, 6 => Dxy, 7 => e-tau inside filter, 8 => mu-tau inside filter, 9 => Single Tau, 10 => VBF DiTau, 11 => di-tau, 12 => e-tau, 13 => mu-tau, 14 => di-tau + PFJet, 15 => e-tau displaced, 16 => mu-tau displaced, 17 => di-tau displaced, 18 => Monitoring, 19 => VBF SingleTau Monitoring, 20 => DiTau+Jet Monitoring, 21 => Monitoring muTau displaced, 22 => OneProng, 23 => DiTau Monitoring, 24 => OverlapFilter, 25 => VBF DiTau monitoring, 26 => SingleTau Monitoring, 27 => MatchL1HLT, 28 => HPS, 29 => single PF-tau inside filter, 30 => VBF SingleTau for Tau;  """
TrigObj_id = NanoAODQuantity("TrigObj_id")
"""dtype: UShort_t; description: ID of the object: 1515 = BoostedTau, 11 = Electron(PixelMatched e/gamma), 6 = FatJet, 3 = HT, 1 = Jet, 2 = MET, 4 = MHT, 13 = Muon, 22 = Photon, 15 = Tau """
TrigObj_l1charge = NanoAODQuantity("TrigObj_l1charge")
"""dtype: Short_t; description: charge of associated L1 seed """
TrigObj_l1iso = NanoAODQuantity("TrigObj_l1iso")
"""dtype: Int_t; description: iso of associated L1 seed """
TrigObj_l1pt = NanoAODQuantity("TrigObj_l1pt")
"""dtype: Float_t; description: pt of associated L1 seed """
TrigObj_l1pt_2 = NanoAODQuantity("TrigObj_l1pt_2")
"""dtype: Float_t; description: pt of associated secondary L1 seed """
TrigObj_l2pt = NanoAODQuantity("TrigObj_l2pt")
"""dtype: Float_t; description: pt of associated 'L2' seed (i.e. HLT before tracking/PF) """
TrigObj_phi = NanoAODQuantity("TrigObj_phi")
"""dtype: Float_t; description: phi """
TrigObj_pt = NanoAODQuantity("TrigObj_pt")
"""dtype: Float_t; description: pt """

TrkMET_phi = NanoAODQuantity("TrkMET_phi")
"""dtype: Float_t; description: raw track MET phi """
TrkMET_pt = NanoAODQuantity("TrkMET_pt")
"""dtype: Float_t; description: raw track MET pt """
TrkMET_sumEt = NanoAODQuantity("TrkMET_sumEt")
"""dtype: Float_t; description: raw track scalar sum of Et """

nboostedTau = NanoAODQuantity("nboostedTau")
"""dtype: Int_t; description: slimmedBoostedTaus after basic selection (pt > 25 && tauID('decayModeFindingNewDMs') && (tauID('byVVLooseIsolationMVArun2DBoldDMwLT') || tauID('byVVLooseIsolationMVArun2DBnewDMwLT') || tauID('byBoostedDeepTau20161718v2p0VSjetraw') > 0.82)) """
boostedTau_charge = NanoAODQuantity("boostedTau_charge")
"""dtype: Int_t; description: electric charge """
boostedTau_chargedIso = NanoAODQuantity("boostedTau_chargedIso")
"""dtype: Float_t; description: charged isolation """
boostedTau_decayMode = NanoAODQuantity("boostedTau_decayMode")
"""dtype: Int_t; description: decayMode() """
boostedTau_eta = NanoAODQuantity("boostedTau_eta")
"""dtype: Float_t; description: eta """
boostedTau_genPartFlav = NanoAODQuantity("boostedTau_genPartFlav")
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==2 taus: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched """
boostedTau_genPartIdx = NanoAODQuantity("boostedTau_genPartIdx")
"""dtype: Short_t; description: Index into genParticle list for MC matching to status==2 taus """
boostedTau_idAntiEle2018 = NanoAODQuantity("boostedTau_idAntiEle2018")
"""dtype: UChar_t; description: Anti-electron MVA discriminator V6 (2018): 1 = VLoose, 2 = Loose, 3 = Medium, 4 = Tight, 5 = VTight """
boostedTau_idAntiMu = NanoAODQuantity("boostedTau_idAntiMu")
"""dtype: UChar_t; description: Anti-muon discriminator V3: : 1 = Loose, 2 = Tight """
boostedTau_idMVAnewDM2017v2 = NanoAODQuantity("boostedTau_idMVAnewDM2017v2")
"""dtype: UChar_t; description: IsolationMVArun2DBnewDMwLT ID working point (2017v2): 1 = VVLoose, 2 = VLoose, 3 = Loose, 4 = Medium, 5 = Tight, 6 = VTight, 7 = VVTight """
boostedTau_idMVAoldDM2017v2 = NanoAODQuantity("boostedTau_idMVAoldDM2017v2")
"""dtype: UChar_t; description: IsolationMVArun2DBoldDMwLT ID working point (2017v2): 1 = VVLoose, 2 = VLoose, 3 = Loose, 4 = Medium, 5 = Tight, 6 = VTight, 7 = VVTight """
boostedTau_jetIdx = NanoAODQuantity("boostedTau_jetIdx")
"""dtype: Short_t; description: index of the associated jet (-1 if none) """
boostedTau_leadTkDeltaEta = NanoAODQuantity("boostedTau_leadTkDeltaEta")
"""dtype: Float_t; description: eta of the leading track, minus tau eta """
boostedTau_leadTkDeltaPhi = NanoAODQuantity("boostedTau_leadTkDeltaPhi")
"""dtype: Float_t; description: phi of the leading track, minus tau phi """
boostedTau_leadTkPtOverTauPt = NanoAODQuantity("boostedTau_leadTkPtOverTauPt")
"""dtype: Float_t; description: pt of the leading track divided by tau pt """
boostedTau_mass = NanoAODQuantity("boostedTau_mass")
"""dtype: Float_t; description: mass """
boostedTau_neutralIso = NanoAODQuantity("boostedTau_neutralIso")
"""dtype: Float_t; description: neutral (photon) isolation """
boostedTau_phi = NanoAODQuantity("boostedTau_phi")
"""dtype: Float_t; description: phi """
boostedTau_photonsOutsideSignalCone = NanoAODQuantity(
    "boostedTau_photonsOutsideSignalCone"
)
"""dtype: Float_t; description: sum of photons outside signal cone """
boostedTau_pt = NanoAODQuantity("boostedTau_pt")
"""dtype: Float_t; description: pt """
boostedTau_puCorr = NanoAODQuantity("boostedTau_puCorr")
"""dtype: Float_t; description: pileup correction """
boostedTau_rawAntiEle2018 = NanoAODQuantity("boostedTau_rawAntiEle2018")
"""dtype: Float_t; description: Anti-electron MVA discriminator V6 raw output discriminator (2018) """
boostedTau_rawAntiEleCat2018 = NanoAODQuantity("boostedTau_rawAntiEleCat2018")
"""dtype: Short_t; description: Anti-electron MVA discriminator V6 category (2018) """
boostedTau_rawBoostedDeepTauRunIIv2p0VSe = NanoAODQuantity(
    "boostedTau_rawBoostedDeepTauRunIIv2p0VSe"
)
"""dtype: Float_t; description: BoostedDeepTau(v2p0) tagger for boostedTaus raw scores Vs e """
boostedTau_rawBoostedDeepTauRunIIv2p0VSjet = NanoAODQuantity(
    "boostedTau_rawBoostedDeepTauRunIIv2p0VSjet"
)
"""dtype: Float_t; description: BoostedDeepTau(v2p0) tagger for boostedTaus raw scores Vs jet """
boostedTau_rawBoostedDeepTauRunIIv2p0VSmu = NanoAODQuantity(
    "boostedTau_rawBoostedDeepTauRunIIv2p0VSmu"
)
"""dtype: Float_t; description: BoostedDeepTau(v2p0) tagger for boostedTaus raw scores Vs mu """
boostedTau_rawIso = NanoAODQuantity("boostedTau_rawIso")
"""dtype: Float_t; description: combined isolation (deltaBeta corrections) """
boostedTau_rawIsodR03 = NanoAODQuantity("boostedTau_rawIsodR03")
"""dtype: Float_t; description: combined isolation (deltaBeta corrections, dR=0.3) """
boostedTau_rawMVAnewDM2017v2 = NanoAODQuantity("boostedTau_rawMVAnewDM2017v2")
"""dtype: Float_t; description: byIsolationMVArun2DBnewDMwLT raw output discriminator (2017v2) """
boostedTau_rawMVAoldDM2017v2 = NanoAODQuantity("boostedTau_rawMVAoldDM2017v2")
"""dtype: Float_t; description: byIsolationMVArun2DBoldDMwLT raw output discriminator (2017v2) """

run = NanoAODQuantity("run")
"""dtype: UInt_t; description: run/i """
luminosityBlock = NanoAODQuantity("luminosityBlock")
"""dtype: UInt_t; description: luminosityBlock/i """
event = NanoAODQuantity("event")
"""dtype: ULong64_t; description: event/l """
bunchCrossing = NanoAODQuantity("bunchCrossing")
"""dtype: UInt_t; description: bunchCrossing/i """
orbitNumber = NanoAODQuantity("orbitNumber")
"""dtype: UInt_t; description: orbitNumber/i """
genWeight = NanoAODQuantity("genWeight")
"""dtype: Float_t; description: generator weight """
LHEWeight_originalXWGTUP = NanoAODQuantity("LHEWeight_originalXWGTUP")
"""dtype: Float_t; description: Nominal event weight in the LHE file """
SoftActivityJetNjets10 = NanoAODQuantity("SoftActivityJetNjets10")
"""dtype: Int_t; description: number of soft activity jet pt, pt >2 """
SoftActivityJetNjets2 = NanoAODQuantity("SoftActivityJetNjets2")
"""dtype: Int_t; description: number of soft activity jet pt, pt >10 """
SoftActivityJetNjets5 = NanoAODQuantity("SoftActivityJetNjets5")
"""dtype: Int_t; description: number of soft activity jet pt, pt >5 """
SoftActivityJetHT = NanoAODQuantity("SoftActivityJetHT")
"""dtype: Float_t; description: scalar sum of soft activity jet pt, pt>1 """
SoftActivityJetHT10 = NanoAODQuantity("SoftActivityJetHT10")
"""dtype: Float_t; description: scalar sum of soft activity jet pt , pt >10 """
SoftActivityJetHT2 = NanoAODQuantity("SoftActivityJetHT2")
"""dtype: Float_t; description: scalar sum of soft activity jet pt, pt >2 """
SoftActivityJetHT5 = NanoAODQuantity("SoftActivityJetHT5")
"""dtype: Float_t; description: scalar sum of soft activity jet pt, pt>5 """
genTtbarId = NanoAODQuantity("genTtbarId")
"""dtype: Int_t; description: ttbar categorization """
L1Reco_step = NanoAODQuantity("L1Reco_step")
"""dtype: Bool_t; description: Trigger/flag bit (process: RECO) """
L1simulation_step = NanoAODQuantity("L1simulation_step")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLTriggerFirstPath = NanoAODQuantity("HLTriggerFirstPath")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
MC_PFScouting = NanoAODQuantity("MC_PFScouting")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """
HLTriggerFinalPath = NanoAODQuantity("HLTriggerFinalPath")
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT) """

TauEmbedding_SelectionNewMass = NanoAODQuantity("TauEmbedding_SelectionNewMass")
"""dtype: Float_t; description: Mass of the Dimuon pair using the new selection algorithm (for internal studies only) """
TauEmbedding_SelectionOldMass = NanoAODQuantity("TauEmbedding_SelectionOldMass")
"""dtype: Float_t; description: Mass of the Dimuon pair using the old selection algorithm (for internal studies only) """
TauEmbedding_initialMETEt = NanoAODQuantity("TauEmbedding_initialMETEt")
"""dtype: Float_t; description: MET Et of selected event """
TauEmbedding_initialMETphi = NanoAODQuantity("TauEmbedding_initialMETphi")
"""dtype: Float_t; description: MET phi of selected event """
TauEmbedding_initialPuppiMETEt = NanoAODQuantity("TauEmbedding_initialPuppiMETEt")
"""dtype: Float_t; description: PuppiMET Et of selected event """
TauEmbedding_initialPuppiMETphi = NanoAODQuantity("TauEmbedding_initialPuppiMETphi")
"""dtype: Float_t; description: PuppiMET phi of selected event """
TauEmbedding_isMediumLeadingMuon = NanoAODQuantity("TauEmbedding_isMediumLeadingMuon")
"""dtype: Bool_t; description: leading muon ID (medium) """
TauEmbedding_isMediumTrailingMuon = NanoAODQuantity("TauEmbedding_isMediumTrailingMuon")
"""dtype: Bool_t; description: trailing muon ID (medium) """
TauEmbedding_isTightLeadingMuon = NanoAODQuantity("TauEmbedding_isTightLeadingMuon")
"""dtype: Bool_t; description: leading muon ID (tight) """
TauEmbedding_isTightTrailingMuon = NanoAODQuantity("TauEmbedding_isTightTrailingMuon")
"""dtype: Bool_t; description: trailing muon ID (tight) """
TauEmbedding_nInitialPairCandidates = NanoAODQuantity(
    "TauEmbedding_nInitialPairCandidates"
)
"""dtype: Float_t; description: number of muons pairs suitable for selection (for internal studies only) """
TauEmbedding_InitialPairCandidates = NanoAODQuantity(
    "TauEmbedding_InitialPairCandidates"
)
"""dtype: Float_t; description: muons pairs suitable for selection (for internal studies only) """