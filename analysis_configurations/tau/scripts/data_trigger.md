# LuminosityBlocks Content

| Collection | Description |
| - | - |
| [**betaStar**](#betastar) | LHC beta star |
| [**crossingAngle**](#crossingangle) | LHC crossing angle |
| [**energy**](#energy) | LHC beam energy |
| [**luminosityBlock**](#luminosityblock) | luminosityBlock/i |
| [**run**](#run) | run/i |

# LuminosityBlocks detail

## <a id='betastar'></a>betaStar [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **betaStar** | Double_t| LHC beta star |

## <a id='crossingangle'></a>crossingAngle [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **crossingAngle** | Double_t| LHC crossing angle |

## <a id='energy'></a>energy [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **energy** | Double_t| LHC beam energy |

## <a id='luminosityblock'></a>luminosityBlock [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **luminosityBlock** | UInt_t| luminosityBlock/i |

## <a id='run'></a>run [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

# Events Content

| Collection | Description |
| - | - |
| [**CaloMET**](#calomet) | phi |
| [**ChsMET**](#chsmet) | raw chs PF MET phi |
| [**CorrT1METJet**](#corrt1metjet) | Additional low-pt jets for Type-1 MET re-correction |
| [**DeepMETResolutionTune**](#deepmetresolutiontune) | DeepmET ResolutionTune phi |
| [**DeepMETResponseTune**](#deepmetresponsetune) | DeepMET ResponseTune phi |
| [**Electron**](#electron) | slimmedElectrons after basic selection (pt > 5 ) |
| [**FatJet**](#fatjet) | slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis |
| [**Flag**](#flag) | Trigger/flag bit (process: PAT) |
| [**FsrPhoton**](#fsrphoton) | Final state radiation photons emitted by muons |
| [**HLT**](#hlt) | Trigger/flag bit (process: HLT) |
| [**HLTriggerFinalPath**](#hltriggerfinalpath) | Trigger/flag bit (process: HLT) |
| [**HLTriggerFirstPath**](#hltriggerfirstpath) | Trigger/flag bit (process: HLT) |
| [**IsoTrack**](#isotrack) | isolated tracks after basic selection (((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && ((abs(dxy) < 0.2 && abs(dz) < 0.1) || pt>15) && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) and lepton veto |
| [**Jet**](#jet) | slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15) |
| [**L1**](#l1) | Trigger/flag bit (process: NANO) |
| [**L1PreFiringWeight**](#l1prefiringweight) | L1 pre-firing event correction weight (1-probability), down var. |
| [**L1Reco**](#l1reco) | Trigger/flag bit (process: RECO) |
| [**LowPtElectron**](#lowptelectron) | slimmedLowPtElectrons after basic selection (pt > 1. && userFloat('ID') > -0.25) |
| [**MET**](#met) | Delta (METx_mod-METx) Unclustered Energy Up |
| [**Muon**](#muon) | slimmedMuons after basic selection (pt > 15 || (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))) |
| [**OtherPV**](#otherpv) | Z position of other primary vertices, excluding the main PV |
| [**PPSLocalTrack**](#ppslocaltrack) | ppsLocalTrack variables |
| [**PV**](#pv) | main primary vertex number of degree of freedom |
| [**Photon**](#photon) | slimmedPhotons after basic selection (pt > 5 ) |
| [**Proton**](#proton) | Mandelstam variable t |
| [**PuppiMET**](#puppimet) | phi |
| [**RawMET**](#rawmet) | phi |
| [**RawPuppiMET**](#rawpuppimet) | phi |
| [**SV**](#sv) | decay length in cm |
| [**SoftActivityJet**](#softactivityjet) | jets clustered from charged candidates compatible with primary vertex (charge()!=0 && pvAssociationQuality()>=5 && vertexRef().key()==0) |
| [**SoftActivityJetHT**](#softactivityjetht) | scalar sum of soft activity jet pt, pt>1 |
| [**SoftActivityJetHT10**](#softactivityjetht10) | scalar sum of soft activity jet pt , pt >10 |
| [**SoftActivityJetHT2**](#softactivityjetht2) | scalar sum of soft activity jet pt, pt >2 |
| [**SoftActivityJetHT5**](#softactivityjetht5) | scalar sum of soft activity jet pt, pt>5 |
| [**SoftActivityJetNjets10**](#softactivityjetnjets10) | number of soft activity jet pt, pt >2 |
| [**SoftActivityJetNjets2**](#softactivityjetnjets2) | number of soft activity jet pt, pt >10 |
| [**SoftActivityJetNjets5**](#softactivityjetnjets5) | number of soft activity jet pt, pt >5 |
| [**SubJet**](#subjet) | slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis |
| [**Tau**](#tau) | slimmedTaus after basic selection (pt > 18 && tauID('decayModeFindingNewDMs') && (tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') || (tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum'))<2.5) || tauID('byVVVLooseDeepTau2017v2p1VSjet'))) |
| [**TkMET**](#tkmet) | raw track MET phi |
| [**TrigObj**](#trigobj) | pt |
| [**boostedTau**](#boostedtau) | slimmedBoostedTaus after basic selection (pt > 40 && tauID('decayModeFindingNewDMs') && (tauID('byVVLooseIsolationMVArun2017v2DBoldDMwLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBnewDMwLT2017'))) |
| [**event**](#event) | event/l |
| [**fixedGridRhoFastjetAll**](#fixedgridrhofastjetall) | rho from all PF Candidates, used e.g. for JECs |
| [**fixedGridRhoFastjetCentral**](#fixedgridrhofastjetcentral) | rho from all PF Candidates for central region, used e.g. for JECs |
| [**fixedGridRhoFastjetCentralCalo**](#fixedgridrhofastjetcentralcalo) | rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation |
| [**fixedGridRhoFastjetCentralChargedPileUp**](#fixedgridrhofastjetcentralchargedpileup) | rho from charged PF Candidates for central region, used e.g. for JECs |
| [**fixedGridRhoFastjetCentralNeutral**](#fixedgridrhofastjetcentralneutral) | rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations |
| [**luminosityBlock**](#luminosityblock) | luminosityBlock/i |
| [**run**](#run) | run/i |

# Events detail

## <a id='calomet'></a>CaloMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **CaloMET_phi** | Float_t| phi |
| **CaloMET_pt** | Float_t| pt |
| **CaloMET_sumEt** | Float_t| scalar sum of Et |

## <a id='chsmet'></a>ChsMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **ChsMET_phi** | Float_t| raw chs PF MET phi |
| **ChsMET_pt** | Float_t| raw chs PF MET pt |
| **ChsMET_sumEt** | Float_t| raw chs PF scalar sum of Et |

## <a id='corrt1metjet'></a>CorrT1METJet [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **CorrT1METJet_area** | Float_t| jet catchment area, for JECs |
| **CorrT1METJet_eta** | Float_t| eta |
| **CorrT1METJet_muonSubtrFactor** | Float_t| 1-(muon-subtracted raw pt)/(raw pt) |
| **CorrT1METJet_phi** | Float_t| phi |
| **CorrT1METJet_rawPt** | Float_t| pt()*jecFactor('Uncorrected') |
| **nCorrT1METJet** | UInt_t| Additional low-pt jets for Type-1 MET re-correction |

## <a id='deepmetresolutiontune'></a>DeepMETResolutionTune [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **DeepMETResolutionTune_phi** | Float_t| DeepmET ResolutionTune phi |
| **DeepMETResolutionTune_pt** | Float_t| DeepMET ResolutionTune pt |

## <a id='deepmetresponsetune'></a>DeepMETResponseTune [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **DeepMETResponseTune_phi** | Float_t| DeepMET ResponseTune phi |
| **DeepMETResponseTune_pt** | Float_t| DeepMET ResponseTune pt |

## <a id='electron'></a>Electron [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Electron_charge** | Int_t| electric charge |
| **Electron_cleanmask** | UChar_t| simple cleaning mask with priority to leptons |
| **Electron_convVeto** | Bool_t| pass conversion veto |
| **Electron_cutBased** | Int_t| cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight) |
| **Electron_cutBased_HEEP** | Bool_t| cut-based HEEP ID |
| **Electron_dEscaleDown** | Float_t| ecal energy scale shifted 1 sigma down (adding gain/stat/syst in quadrature) |
| **Electron_dEscaleUp** | Float_t| ecal energy scale shifted 1 sigma up(adding gain/stat/syst in quadrature) |
| **Electron_dEsigmaDown** | Float_t| ecal energy smearing value shifted 1 sigma up |
| **Electron_dEsigmaUp** | Float_t| ecal energy smearing value shifted 1 sigma up |
| **Electron_deltaEtaSC** | Float_t| delta eta (SC,ele) with sign |
| **Electron_dr03EcalRecHitSumEt** | Float_t| Non-PF Ecal isolation within a delta R cone of 0.3 with electron pt > 35 GeV |
| **Electron_dr03HcalDepth1TowerSumEt** | Float_t| Non-PF Hcal isolation within a delta R cone of 0.3 with electron pt > 35 GeV |
| **Electron_dr03TkSumPt** | Float_t| Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV |
| **Electron_dr03TkSumPtHEEP** | Float_t| Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV used in HEEP ID |
| **Electron_dxy** | Float_t| dxy (with sign) wrt first PV, in cm |
| **Electron_dxyErr** | Float_t| dxy uncertainty, in cm |
| **Electron_dz** | Float_t| dz (with sign) wrt first PV, in cm |
| **Electron_dzErr** | Float_t| dz uncertainty, in cm |
| **Electron_eCorr** | Float_t| ratio of the calibrated energy/miniaod energy |
| **Electron_eInvMinusPInv** | Float_t| 1/E_SC - 1/p_trk |
| **Electron_energyErr** | Float_t| energy error of the cluster-track combination |
| **Electron_eta** | Float_t| eta |
| **Electron_hoe** | Float_t| H over E |
| **Electron_ip3d** | Float_t| 3D impact parameter wrt first PV, in cm |
| **Electron_isPFcand** | Bool_t| electron is PF candidate |
| **Electron_jetIdx** | Int_t(index to Jet)| index of the associated jet (-1 if none) |
| **Electron_jetNDauCharged** | UChar_t| number of charged daughters of the closest jet |
| **Electron_jetPtRelv2** | Float_t| Relative momentum of the lepton with respect to the closest jet after subtracting the lepton |
| **Electron_jetRelIso** | Float_t| Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet) |
| **Electron_lostHits** | UChar_t| number of missing inner hits |
| **Electron_mass** | Float_t| mass |
| **Electron_miniPFRelIso_all** | Float_t| mini PF relative isolation, total (with scaled rho*EA PU corrections) |
| **Electron_miniPFRelIso_chg** | Float_t| mini PF relative isolation, charged component |
| **Electron_mvaFall17V2Iso** | Float_t| MVA Iso ID V2 score |
| **Electron_mvaFall17V2Iso_WP80** | Bool_t| MVA Iso ID V2 WP80 |
| **Electron_mvaFall17V2Iso_WP90** | Bool_t| MVA Iso ID V2 WP90 |
| **Electron_mvaFall17V2Iso_WPL** | Bool_t| MVA Iso ID V2 loose WP |
| **Electron_mvaFall17V2noIso** | Float_t| MVA noIso ID V2 score |
| **Electron_mvaFall17V2noIso_WP80** | Bool_t| MVA noIso ID V2 WP80 |
| **Electron_mvaFall17V2noIso_WP90** | Bool_t| MVA noIso ID V2 WP90 |
| **Electron_mvaFall17V2noIso_WPL** | Bool_t| MVA noIso ID V2 loose WP |
| **Electron_mvaTTH** | Float_t| TTH MVA lepton ID score |
| **Electron_pdgId** | Int_t| PDG code assigned by the event reconstruction (not by MC truth) |
| **Electron_pfRelIso03_all** | Float_t| PF relative isolation dR=0.3, total (with rho*EA PU corrections) |
| **Electron_pfRelIso03_chg** | Float_t| PF relative isolation dR=0.3, charged component |
| **Electron_phi** | Float_t| phi |
| **Electron_photonIdx** | Int_t(index to Photon)| index of the associated photon (-1 if none) |
| **Electron_pt** | Float_t| p_{T} |
| **Electron_r9** | Float_t| R9 of the supercluster, calculated with full 5x5 region |
| **Electron_scEtOverPt** | Float_t| (supercluster transverse energy)/pt-1 |
| **Electron_seedGain** | UChar_t| Gain of the seed crystal |
| **Electron_sieie** | Float_t| sigma_IetaIeta of the supercluster, calculated with full 5x5 region |
| **Electron_sip3d** | Float_t| 3D impact parameter significance wrt first PV, in cm |
| **Electron_tightCharge** | Int_t| Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent) |
| **Electron_vidNestedWPBitmap** | Int_t| VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaCut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEInverseMinusPInverseCut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 3 bits per cut |
| **Electron_vidNestedWPBitmapHEEP** | Int_t| VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaWithSatCut,GsfEleFull5x5E2x5OverE5x5WithSatCut,GsfEleHadronicOverEMLinearCut,GsfEleTrkPtIsoCut,GsfEleEmHadD1IsoRhoCut,GsfEleDxyCut,GsfEleMissingHitsCut,GsfEleEcalDrivenCut), 1 bits per cut |
| **nElectron** | UInt_t| slimmedElectrons after basic selection (pt > 5 ) |

## <a id='fatjet'></a>FatJet [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **FatJet_area** | Float_t| jet catchment area, for JECs |
| **FatJet_btagCSVV2** | Float_t|  pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2) |
| **FatJet_btagDDBvLV2** | Float_t| DeepDoubleX V2(mass-decorrelated) discriminator for H(Z)->bb vs QCD |
| **FatJet_btagDDCvBV2** | Float_t| DeepDoubleX V2 (mass-decorrelated) discriminator for H(Z)->cc vs H(Z)->bb |
| **FatJet_btagDDCvLV2** | Float_t| DeepDoubleX V2 (mass-decorrelated) discriminator for H(Z)->cc vs QCD |
| **FatJet_btagDeepB** | Float_t| DeepCSV b+bb tag discriminator |
| **FatJet_btagHbb** | Float_t| Higgs to BB tagger discriminator |
| **FatJet_deepTagMD_H4qvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger H->4q vs QCD discriminator |
| **FatJet_deepTagMD_HbbvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger H->bb vs QCD discriminator |
| **FatJet_deepTagMD_TvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger top vs QCD discriminator |
| **FatJet_deepTagMD_WvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger W vs QCD discriminator |
| **FatJet_deepTagMD_ZHbbvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger Z/H->bb vs QCD discriminator |
| **FatJet_deepTagMD_ZHccvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger Z/H->cc vs QCD discriminator |
| **FatJet_deepTagMD_ZbbvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger Z->bb vs QCD discriminator |
| **FatJet_deepTagMD_ZvsQCD** | Float_t| Mass-decorrelated DeepBoostedJet tagger Z vs QCD discriminator |
| **FatJet_deepTagMD_bbvsLight** | Float_t| Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->bb vs light flavour discriminator |
| **FatJet_deepTagMD_ccvsLight** | Float_t| Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->cc vs light flavour discriminator |
| **FatJet_deepTag_H** | Float_t| DeepBoostedJet tagger H(bb,cc,4q) sum |
| **FatJet_deepTag_QCD** | Float_t| DeepBoostedJet tagger QCD(bb,cc,b,c,others) sum |
| **FatJet_deepTag_QCDothers** | Float_t| DeepBoostedJet tagger QCDothers value |
| **FatJet_deepTag_TvsQCD** | Float_t| DeepBoostedJet tagger top vs QCD discriminator |
| **FatJet_deepTag_WvsQCD** | Float_t| DeepBoostedJet tagger W vs QCD discriminator |
| **FatJet_deepTag_ZvsQCD** | Float_t| DeepBoostedJet tagger Z vs QCD discriminator |
| **FatJet_electronIdx3SJ** | Int_t(index to Electron)| index of electron matched to jet |
| **FatJet_eta** | Float_t| eta |
| **FatJet_jetId** | Int_t| Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto |
| **FatJet_lsf3** | Float_t| Lepton Subjet Fraction (3 subjets) |
| **FatJet_mass** | Float_t| mass |
| **FatJet_msoftdrop** | Float_t| Corrected soft drop mass with PUPPI |
| **FatJet_muonIdx3SJ** | Int_t(index to Muon)| index of muon matched to jet |
| **FatJet_n2b1** | Float_t| N2 with beta=1 |
| **FatJet_n3b1** | Float_t| N3 with beta=1 |
| **FatJet_nConstituents** | UChar_t| Number of particles in the jet |
| **FatJet_particleNetMD_QCD** | Float_t| Mass-decorrelated ParticleNet tagger raw QCD score |
| **FatJet_particleNetMD_Xbb** | Float_t| Mass-decorrelated ParticleNet tagger raw X->bb score. For X->bb vs QCD tagging, use Xbb/(Xbb+QCD) |
| **FatJet_particleNetMD_Xcc** | Float_t| Mass-decorrelated ParticleNet tagger raw X->cc score. For X->cc vs QCD tagging, use Xcc/(Xcc+QCD) |
| **FatJet_particleNetMD_Xqq** | Float_t| Mass-decorrelated ParticleNet tagger raw X->qq (uds) score. For X->qq vs QCD tagging, use Xqq/(Xqq+QCD). For W vs QCD tagging, use (Xcc+Xqq)/(Xcc+Xqq+QCD) |
| **FatJet_particleNet_H4qvsQCD** | Float_t| ParticleNet tagger H(->VV->qqqq) vs QCD discriminator |
| **FatJet_particleNet_HbbvsQCD** | Float_t| ParticleNet tagger H(->bb) vs QCD discriminator |
| **FatJet_particleNet_HccvsQCD** | Float_t| ParticleNet tagger H(->cc) vs QCD discriminator |
| **FatJet_particleNet_QCD** | Float_t| ParticleNet tagger QCD(bb,cc,b,c,others) sum |
| **FatJet_particleNet_TvsQCD** | Float_t| ParticleNet tagger top vs QCD discriminator |
| **FatJet_particleNet_WvsQCD** | Float_t| ParticleNet tagger W vs QCD discriminator |
| **FatJet_particleNet_ZvsQCD** | Float_t| ParticleNet tagger Z vs QCD discriminator |
| **FatJet_particleNet_mass** | Float_t| ParticleNet mass regression |
| **FatJet_phi** | Float_t| phi |
| **FatJet_pt** | Float_t| pt |
| **FatJet_rawFactor** | Float_t| 1 - Factor to get back to raw pT |
| **FatJet_subJetIdx1** | Int_t(index to Subjet)| index of first subjet |
| **FatJet_subJetIdx2** | Int_t(index to Subjet)| index of second subjet |
| **FatJet_tau1** | Float_t| Nsubjettiness (1 axis) |
| **FatJet_tau2** | Float_t| Nsubjettiness (2 axis) |
| **FatJet_tau3** | Float_t| Nsubjettiness (3 axis) |
| **FatJet_tau4** | Float_t| Nsubjettiness (4 axis) |
| **nFatJet** | UInt_t| slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis |

## <a id='flag'></a>Flag [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Flag_BadChargedCandidateFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadChargedCandidateFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_BadChargedCandidateSummer16Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadChargedCandidateSummer16Filter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_BadPFMuonDzFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadPFMuonFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadPFMuonFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_BadPFMuonSummer16Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadPFMuonSummer16Filter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_CSCTightHalo2015Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_CSCTightHalo2015Filter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_CSCTightHaloFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_CSCTightHaloFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_CSCTightHaloTrkMuUnvetoFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_EcalDeadCellBoundaryEnergyFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_EcalDeadCellBoundaryEnergyFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_EcalDeadCellTriggerPrimitiveFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_HBHENoiseFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_HBHENoiseFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_HBHENoiseIsoFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_HBHENoiseIsoFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_HcalStripHaloFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_HcalStripHaloFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_METFilters** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_METFilters_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_chargedHadronTrackResolutionFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_chargedHadronTrackResolutionFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_ecalBadCalibFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_ecalBadCalibFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_ecalLaserCorrFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_ecalLaserCorrFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_eeBadScFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_eeBadScFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_globalSuperTightHalo2016Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_globalSuperTightHalo2016Filter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_globalTightHalo2016Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_globalTightHalo2016Filter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_goodVertices** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_goodVertices_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_hcalLaserEventFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_hcalLaserEventFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_hfNoisyHitsFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_muonBadTrackFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_muonBadTrackFilter_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_trkPOGFilters** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOGFilters_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_trkPOG_logErrorTooManyClusters** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOG_logErrorTooManyClusters_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_trkPOG_manystripclus53X** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOG_manystripclus53X_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |
| **Flag_trkPOG_toomanystripclus53X** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOG_toomanystripclus53X_pRECO** | Bool_t| Trigger/flag bit (process: RECO) |

## <a id='fsrphoton'></a>FsrPhoton [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **FsrPhoton_dROverEt2** | Float_t| deltaR to associated muon divided by photon et2 |
| **FsrPhoton_eta** | Float_t| eta |
| **FsrPhoton_muonIdx** | Int_t(index to Muon)| index of associated muon |
| **FsrPhoton_phi** | Float_t| phi |
| **FsrPhoton_pt** | Float_t| pt |
| **FsrPhoton_relIso03** | Float_t| relative isolation in a 0.3 cone without CHS |
| **nFsrPhoton** | UInt_t| Final state radiation photons emitted by muons |

## <a id='hlt'></a>HLT [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HLT_AK4CaloJet100** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4CaloJet120** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4CaloJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4CaloJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4CaloJet50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4CaloJet80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4PFJet100** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4PFJet120** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4PFJet50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK4PFJet80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFHT750_TrimMass50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFHT800_TrimMass50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFHT850_TrimMass50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFHT900_TrimMass50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet140** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet25** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet260** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet320** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet360_TrimMass30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet380_TrimMass30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet400** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet400_TrimMass30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet420_TrimMass30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet500** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet550** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJet80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd140** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd25** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd260** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd320** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd400** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd500** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_AK8PFJetFwd80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK4DiJet110_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK4DiJet170_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK4DiJet20_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK4DiJet40_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK4DiJet70_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK4Jet300_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK8DiJet170_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK8Jet170_DoubleMu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_BTagMu_AK8Jet300_Mu5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CDC_L2cosmic_5_er1p0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CDC_L2cosmic_5p5_er1p0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloJet500_NoJetID** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloJet550_NoJetID** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET100_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET100_NotCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET110_NotCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET250_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET250_NotCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET300_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET350_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET70_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET80_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET80_NotCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET90_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMET90_NotCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_CaloMHT90** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiEle27_WPTightCaloOnly_L1DoubleEG** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiJet110_35_Mjj650_PFMET110** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiJet110_35_Mjj650_PFMET120** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiJet110_35_Mjj650_PFMET130** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiMu9_Ele9_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve100_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve140** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve15_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve160_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve220_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve25_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve260** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve300_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve320** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve400** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve500** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve60_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiPFJetAve80_HFJEC** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiSC30_18_EIso_AND_HE_Mass70** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi3p5_Muon2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi_L1_4R_0er1p5R** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi_L1_NoOS** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi_NoVertexing** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Jpsi_NoVertexing_NoOS** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_LowMass** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_LowMass_L1_0er1p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_LowMass_L1_0er1p5R** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_LowMass_L1_4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_LowMass_L1_4R** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_LowMass_L1_TM530** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_L1_4p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_L1_4p5NoOS** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_L1_4p5er2p0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_L1_4p5er2p0M** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_L1_5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_L1_5M** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_Muon_L1_TM0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_Muon_NoL1Mass** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon0_Upsilon_NoVertexing** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon10_PsiPrime_Barrel_Seagulls** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon12_Upsilon_y1p4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon14_Phi_Barrel_Seagulls** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon18_PsiPrime** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon18_PsiPrime_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon20_Jpsi_Barrel_Seagulls** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon24_Phi_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon24_Upsilon_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon25_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Dimuon25_Jpsi_noCorrL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle24_eta2p1_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle25_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle27_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle33_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleIsoMu20_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu23NoVtx_2Cha** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu25NoVtx_2Cha** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleL2Mu50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu20_7_Mass0to30_L1_DM4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu20_7_Mass0to30_L1_DM4EG** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu20_7_Mass0to30_Photon23** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu33NoFiltersNoVtxDisplaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_DCA_PFMET50_PFMHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_DZ_PFMET50_PFMHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_DZ_PFMET70_PFMHT70** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_DZ_PFMET90_PFMHT90** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_TkMu_DsTau3Mu** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_Trk_Tau3mu** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu40NoFiltersNoVtxDisplaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu43NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu48NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_3_Bs** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_3_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_JpsiTrkTrk_Displaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_JpsiTrk_Displaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_Jpsi_Displaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_Jpsi_NoVertexing** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_LowMassNonResonantTrk_Displaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_Mass3p8_DZ_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu4_PsiPrimeTrk_Displaced** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePFJets100_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePFJets200_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePFJets350_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePFJets40_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePhoton33_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePhoton70** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoublePhoton85** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleTrkMu_16_6_NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ECALHT800** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_EcalCalibration** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele115_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele135_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele145_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT450_PFMET50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT600** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_WPLoose_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele17_CaloIdM_TrackIdM_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele17_WPLoose_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele200_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele20_WPLoose_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele20_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele20_eta2p1_WPLoose_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_CaloIdM_TrackIdM_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele250_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele27_Ele37_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele27_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele28_HighEta_SC20_Mass55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele28_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele28_eta2p1_WPTight_Gsf_HT150** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele300_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele30_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele32_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele32_WPTight_Gsf_L1DoubleEG** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele35_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele35_WPTight_Gsf_L1EGMT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele38_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele40_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele50_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele8_CaloIdM_TrackIdM_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT300_Beamspot** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT400_DisplacedDijet40_DisplacedTrack** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT425** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT430_DisplacedDijet40_DisplacedTrack** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT430_DisplacedDijet60_DisplacedTrack** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT450_Beamspot** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT500_DisplacedDijet40_DisplacedTrack** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT550_DisplacedDijet60_Inclusive** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HT650_DisplacedDijet60_Inclusive** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HcalCalibration** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HcalIsolatedbunch** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HcalNZS** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_HcalPhiSym** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_TwoProngs35** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_MET90** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoTrackHB** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoTrackHE** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1ETMHadSeeds** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1NotBptxOR** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1SingleMu18** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1SingleMu25** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1UnpairedBunchBptxMinus** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1UnpairedBunchBptxPlus** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu10** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu10_NoVertex_NoBPTX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu10_NoVertex_NoBPTX3BX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu23NoVtx_2Cha** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu23NoVtx_2Cha_CosmicSeed** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_L2Mu50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MET105_IsoTrk50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MET120_IsoTrk50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePhoton20** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_IP6_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_IP6_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_IP6_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_IP6_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_IP6_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT450_PFMET50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT600** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_Photon30_IsoCaloId** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu18_Mu9** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu18_Mu9_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu18_Mu9_SameSign** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu18_Mu9_SameSign_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20_Mu10** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20_Mu10_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20_Mu10_SameSign** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20_Mu10_SameSign_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20_TkMu0_Phi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_Mu12** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_Mu12_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_Mu12_SameSign** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_Mu12_SameSign_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu25_TkMu0_Onia** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu25_TkMu0_Phi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu27** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu27_Ele37_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu30_TkMu0_Psi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu30_TkMu0_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu37_Ele27_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu37_TkMu27** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3_PFJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu50_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7_IP4_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7_IP4_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7_IP4_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7_IP4_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7_IP4_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_L2Mu2_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_L2Mu2_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_Track2_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_Track2_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_Track3p5_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_Track3p5_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_Track7_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_Track7_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_DiEle12_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP3_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP3_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP3_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP3_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP3_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP5_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP5_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP5_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP5_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP5_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP6_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP6_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP6_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP6_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_IP6_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP4_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP4_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP4_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP4_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP4_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP5_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP5_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP5_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP5_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP5_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP6_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP6_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP6_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP6_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu9_IP6_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_OldMu100** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT1050** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT180** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT250** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT330PT30_QuadPFJet_75_60_45_40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT350MinPFJet15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT370** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT400_SixPFJet32** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT430** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT450_SixPFJet36** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT500_PFMET100_PFMHT100_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT500_PFMET110_PFMHT110_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT510** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT590** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT680** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT700_PFMET85_PFMHT85_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT700_PFMET95_PFMHT95_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT780** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT800_PFMET75_PFMHT75_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT800_PFMET85_PFMHT85_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFHT890** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet140** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet25** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet260** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet320** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet400** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet500** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet550** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJet80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd140** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd25** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd260** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd320** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd400** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd500** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFJetFwd80** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET100_PFMHT100_IDTight_PFHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET110_PFMHT110_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET120_PFMHT120_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET120_PFMHT120_IDTight_PFHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET130_PFMHT130_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET140_PFMHT140_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET200_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET200_HBHE_BeamHaloCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET200_NotCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET250_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMET300_HBHECleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETNoMu110_PFMHTNoMu110_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETNoMu120_PFMHTNoMu120_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETNoMu130_PFMHTNoMu130_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETNoMu140_PFMHTNoMu140_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne110_PFMHT110_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne120_PFMHT120_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne130_PFMHT130_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne140_PFMHT140_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon100EBHE10** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon100EB_TightID_TightIso** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon100EEHE10** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon100EE_TightID_TightIso** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon110EB_TightID_TightIso** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon120** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon120EB_TightID_TightIso** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon120_R9Id90_HE10_IsoM** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon150** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon165_R9Id90_HE10_IsoM** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon175** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon20** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon20_HoverELoose** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon300_NoHE** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon30_HoverELoose** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon33** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon35_TwoProngs35** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon50_R9Id90_HE10_IsoM** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon60_R9Id90_CaloIdL_IsoL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon75** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon75_R9Id90_HE10_IsoM** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon90** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Photon90_R9Id90_HE10_IsoM** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Physics_part7** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet103_88_75_15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet105_88_76_15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet111_90_80_15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet98_83_71_15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Random** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Rsq0p35** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Rsq0p40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_RsqMR300_Rsq0p09_MR200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_RsqMR300_Rsq0p09_MR200_4jet** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_RsqMR320_Rsq0p09_MR200** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_RsqMR320_Rsq0p09_MR200_4jet** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_SingleJet30_Mu12_SinglePFJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_SinglePhoton10_Eta3p1ForPPRef** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_SinglePhoton20_Eta3p1ForPPRef** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_SinglePhoton30_Eta3p1ForPPRef** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TkMu100** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Trimuon5_3p5_2_Upsilon_Muon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleJet110_35_35_Mjj650_PFMET110** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleJet110_35_35_Mjj650_PFMET120** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleJet110_35_35_Mjj650_PFMET130** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_10_5_5_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_12_10_5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_5_3_3_Mass3p8_DCA** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_5_3_3_Mass3p8_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TriplePhoton_20_20_20_CaloIdLV2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TriplePhoton_30_30_10_CaloIdLV2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TrkMu16NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TrkMu6NoFiltersNoVtx** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_UncorrectedJetE30_NoBPTX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_UncorrectedJetE30_NoBPTX3BX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_UncorrectedJetE60_NoBPTX3BX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_UncorrectedJetE70_NoBPTX3BX** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_Alignment** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_Beamspot** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_FirstBXAfterTrain** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_FirstCollisionAfterAbortGap** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_FirstCollisionInTrain** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_IsolatedBunches** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_LastCollisionInTrain** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part0** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part2** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part4** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_ZeroBias_part7** | Bool_t| Trigger/flag bit (process: HLT) |

## <a id='hltriggerfinalpath'></a>HLTriggerFinalPath [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HLTriggerFinalPath** | Bool_t| Trigger/flag bit (process: HLT) |

## <a id='hltriggerfirstpath'></a>HLTriggerFirstPath [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **HLTriggerFirstPath** | Bool_t| Trigger/flag bit (process: HLT) |

## <a id='isotrack'></a>IsoTrack [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **IsoTrack_charge** | Int_t| electric charge |
| **IsoTrack_dxy** | Float_t| dxy (with sign) wrt first PV, in cm |
| **IsoTrack_dz** | Float_t| dz (with sign) wrt first PV, in cm |
| **IsoTrack_eta** | Float_t| eta |
| **IsoTrack_fromPV** | Int_t| isolated track comes from PV |
| **IsoTrack_isFromLostTrack** | Bool_t| if isolated track comes from a lost track |
| **IsoTrack_isHighPurityTrack** | Bool_t| track is high purity |
| **IsoTrack_isPFcand** | Bool_t| if isolated track is a PF candidate |
| **IsoTrack_miniPFRelIso_all** | Float_t| mini PF relative isolation, total (with scaled rho*EA PU corrections) |
| **IsoTrack_miniPFRelIso_chg** | Float_t| mini PF relative isolation, charged component |
| **IsoTrack_pdgId** | Int_t| PDG id of PF cand |
| **IsoTrack_pfRelIso03_all** | Float_t| PF relative isolation dR=0.3, total (deltaBeta corrections) |
| **IsoTrack_pfRelIso03_chg** | Float_t| PF relative isolation dR=0.3, charged component |
| **IsoTrack_phi** | Float_t| phi |
| **IsoTrack_pt** | Float_t| pt |
| **nIsoTrack** | UInt_t| isolated tracks after basic selection (((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && ((abs(dxy) < 0.2 && abs(dz) < 0.1) || pt>15) && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) and lepton veto |

## <a id='jet'></a>Jet [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Jet_area** | Float_t| jet catchment area, for JECs |
| **Jet_bRegCorr** | Float_t| pt correction for b-jet energy regression |
| **Jet_bRegRes** | Float_t| res on pt corrected with b-jet regression |
| **Jet_btagCSVV2** | Float_t|  pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2) |
| **Jet_btagDeepB** | Float_t| DeepCSV b+bb tag discriminator |
| **Jet_btagDeepCvB** | Float_t| DeepCSV c vs b+bb discriminator |
| **Jet_btagDeepCvL** | Float_t| DeepCSV c vs udsg discriminator |
| **Jet_btagDeepFlavB** | Float_t| DeepJet b+bb+lepb tag discriminator |
| **Jet_btagDeepFlavCvB** | Float_t| DeepJet c vs b+bb+lepb discriminator |
| **Jet_btagDeepFlavCvL** | Float_t| DeepJet c vs uds+g discriminator |
| **Jet_btagDeepFlavQG** | Float_t| DeepJet g vs uds discriminator |
| **Jet_cRegCorr** | Float_t| pt correction for c-jet energy regression |
| **Jet_cRegRes** | Float_t| res on pt corrected with c-jet regression |
| **Jet_chEmEF** | Float_t| charged Electromagnetic Energy Fraction |
| **Jet_chFPV0EF** | Float_t| charged fromPV==0 Energy Fraction (energy excluded from CHS jets). Previously called betastar. |
| **Jet_chHEF** | Float_t| charged Hadron Energy Fraction |
| **Jet_cleanmask** | UChar_t| simple cleaning mask with priority to leptons |
| **Jet_electronIdx1** | Int_t(index to Electron)| index of first matching electron |
| **Jet_electronIdx2** | Int_t(index to Electron)| index of second matching electron |
| **Jet_eta** | Float_t| eta |
| **Jet_hfadjacentEtaStripsSize** | Int_t| eta size of the strips next to the central tower strip in HF (noise discriminating variable)  |
| **Jet_hfcentralEtaStripSize** | Int_t| eta size of the central tower strip in HF (noise discriminating variable)  |
| **Jet_hfsigmaEtaEta** | Float_t| sigmaEtaEta for HF jets (noise discriminating variable) |
| **Jet_hfsigmaPhiPhi** | Float_t| sigmaPhiPhi for HF jets (noise discriminating variable) |
| **Jet_jetId** | Int_t| Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto |
| **Jet_mass** | Float_t| mass |
| **Jet_muEF** | Float_t| muon Energy Fraction |
| **Jet_muonIdx1** | Int_t(index to Muon)| index of first matching muon |
| **Jet_muonIdx2** | Int_t(index to Muon)| index of second matching muon |
| **Jet_muonSubtrFactor** | Float_t| 1-(muon-subtracted raw pt)/(raw pt) |
| **Jet_nConstituents** | UChar_t| Number of particles in the jet |
| **Jet_nElectrons** | Int_t| number of electrons in the jet |
| **Jet_nMuons** | Int_t| number of muons in the jet |
| **Jet_neEmEF** | Float_t| neutral Electromagnetic Energy Fraction |
| **Jet_neHEF** | Float_t| neutral Hadron Energy Fraction |
| **Jet_phi** | Float_t| phi |
| **Jet_pt** | Float_t| pt |
| **Jet_puId** | Int_t| Pileup ID flags with 106X (2018) training |
| **Jet_puIdDisc** | Float_t| Pileup ID discriminant with 106X (2018) training |
| **Jet_qgl** | Float_t| Quark vs Gluon likelihood discriminator |
| **Jet_rawFactor** | Float_t| 1 - Factor to get back to raw pT |
| **nJet** | UInt_t| slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15) |

## <a id='l1'></a>L1 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **L1_AlwaysTrue** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_AND_Ref1_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_AND_Ref3_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_AND_Ref4_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_BeamGas_B1_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_BeamGas_B2_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_BeamGas_Ref1_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_BeamGas_Ref2_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_NotOR_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_OR_Ref3_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_OR_Ref4_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BPTX_RefAND_VME** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BptxMinus** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BptxOR** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BptxPlus** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_BptxXOR** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG8er2p5_HTT260er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG8er2p5_HTT280er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG8er2p5_HTT300er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG8er2p5_HTT320er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG8er2p5_HTT340er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_15_10_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_20_10_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_22_10_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_25_12_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_25_14_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_27_14_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_LooseIso20_10_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_LooseIso22_10_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_LooseIso22_12_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleEG_LooseIso25_12_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleIsoTau28er2p1_Mass_Max80** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleIsoTau28er2p1_Mass_Max90** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleIsoTau32er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleIsoTau34er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleIsoTau36er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet100er2p3_dEta_Max1p6** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet100er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet112er2p3_dEta_Max1p6** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet120er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet150er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet40er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_100_30_DoubleJet30_Mass_Min620** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_110_35_DoubleJet35_Mass_Min620** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_115_40_DoubleJet40_Mass_Min620** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_120_45_DoubleJet45_Mass_Min620** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_80_30_Mass_Min420_Mu8** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleJet_90_30_DoubleJet30_Mass_Min620** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleLooseIsoEG22er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleLooseIsoEG24er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0_Mass_Min1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0_OQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0_SQ_OS** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er1p5_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er1p5_SQ_OS** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er1p5_SQ_dR_Max1p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu0er2p0_SQ_dR_Max1p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu10_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu18er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_OS_DoubleEG7p5Upsilon** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_ETMHF50_HTT60er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_HTT220er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_HTT240er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_SQ_HTT260er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4_SQ_EG9er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4_SQ_OS** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4_SQ_OS_dR_Max1p2** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4p5_SQ_OS** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4p5_SQ_OS_dR_Max1p2** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4p5er2p0_SQ_OS** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu5Upsilon_OS_DoubleEG3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu5_SQ_EG9er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu8_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu9_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu_12_5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu_15_5_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu_15_7** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu_15_7_Mass_Min1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleMu_15_7_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_DoubleTau70er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETM120** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETM150** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF100** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF100_HTT60er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF110** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF110_HTT60er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF110_HTT60er_NotSecondBunchInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF120** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF120_HTT60er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF120_NotSecondBunchInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF130** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF130_HTT60er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF140** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF150** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETMHF90_HTT60er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETT1200** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETT1600** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ETT2000** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_FirstBunchAfterTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_FirstBunchBeforeTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_FirstBunchInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_FirstCollisionInOrbit** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_FirstCollisionInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HCAL_LaserMon_Trig** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HCAL_LaserMon_Veto** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT120er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT160er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT200er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT255er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT280er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT280er_QuadJet_70_55_40_35_er2p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT320er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT320er_QuadJet_70_55_40_40_er2p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT360er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT400er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_HTT450er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoEG32er2p5_Mt40** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoEG32er2p5_Mt44** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoEG32er2p5_Mt48** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoTau40er2p1_ETMHF100** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoTau40er2p1_ETMHF110** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoTau40er2p1_ETMHF120** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoTau40er2p1_ETMHF80** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsoTau40er2p1_ETMHF90** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_IsolatedBunch** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LastBunchInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LastCollisionInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG24er2p1_HTT100er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG26er2p1_HTT100er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG28er2p1_HTT100er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG30er2p1_HTT100er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_MinimumBiasHF0_AND_BptxAND** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu18er2p1_Tau24er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu18er2p1_Tau26er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu20_EG10er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu22er2p1_IsoTau32er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu22er2p1_IsoTau34er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu22er2p1_IsoTau36er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu22er2p1_IsoTau40er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu22er2p1_Tau70er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet120er2p5_dR_Max0p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet120er2p5_dR_Max0p8** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet16er2p5_dR_Max0p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet30er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet35er2p5_dR_Max0p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet60er2p5_dR_Max0p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3_Jet80er2p5_dR_Max0p4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3er1p5_Jet100er2p5_ETMHF40** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu3er1p5_Jet100er2p5_ETMHF50** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu5_EG23er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu5_LooseIsoEG20er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu6_DoubleEG10er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu6_DoubleEG12er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu6_DoubleEG15er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu6_DoubleEG17er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu6_HTT240er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu6_HTT250er** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu7_EG20er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu7_EG23er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu7_LooseIsoEG20er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_Mu7_LooseIsoEG23er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_NotBptxOR** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_QuadJet36er2p5_IsoTau52er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_QuadJet60er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_QuadMu0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_QuadMu0_OQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_QuadMu0_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SecondBunchInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SecondLastBunchInTrain** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG10er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG15er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG26er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG28_FWD2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG28er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG28er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG28er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG34er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG36er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG38er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG40er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG42er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG45er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG50** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG60** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleEG8er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG24er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG24er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG26er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG26er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG26er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG28_FWD2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG28er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG28er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG28er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG30er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG30er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG32er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG32er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleIsoEG34er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet10erHE** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet120** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet120_FWD3p0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet120er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet12erHE** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet140er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet140er2p5_ETMHF70** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet140er2p5_ETMHF80** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet140er2p5_ETMHF90** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet160er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet180** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet180er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet200** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet20er2p5_NotBptxOR** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet20er2p5_NotBptxOR_3BX** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet35** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet35_FWD3p0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet35er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet43er2p5_NotBptxOR_3BX** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet46er2p5_NotBptxOR_3BX** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet60** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet60_FWD3p0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet60er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet8erHE** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet90** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet90_FWD3p0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleJet90er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG26er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG26er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG28_FWD2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG28er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG28er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG28er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG30er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleLooseIsoEG30er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu0_BMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu0_DQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu0_EMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu0_OMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu10er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu12_DQ_BMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu12_DQ_EMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu12_DQ_OMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu12er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu14er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu15_DQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu16er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu18** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu18er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu20** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu22** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu22_BMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu22_EMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu22_OMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu25** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu6er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu7** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu7_DQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu7er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu8er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMu9er1p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuCosmics** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuCosmics_BMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuCosmics_EMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuCosmics_OMTF** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuOpen** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuOpen_NotBptxOR** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuOpen_er1p1_NotBptxOR_3BX** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleMuOpen_er1p4_NotBptxOR_3BX** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleTau120er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_SingleTau130er2p1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TOTEM_1** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TOTEM_2** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TOTEM_3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TOTEM_4** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleEG16er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleEG_16_12_8_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleEG_16_15_8_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleEG_18_17_8_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleEG_18_18_12_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu0** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu0_OQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu0_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu3_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5SQ_3SQ_0OQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_3_3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_3_3_SQ** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_3p5_2p5** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_TripleMu_5_5_3** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_UnpairedBunchBptxMinus** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_UnpairedBunchBptxPlus** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_UnprefireableEvent** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ZeroBias** | Bool_t| Trigger/flag bit (process: NANO) |
| **L1_ZeroBias_copy** | Bool_t| Trigger/flag bit (process: NANO) |

## <a id='l1prefiringweight'></a>L1PreFiringWeight [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **L1PreFiringWeight_Dn** | Float_t| L1 pre-firing event correction weight (1-probability), down var. |
| **L1PreFiringWeight_ECAL_Dn** | Float_t| ECAL L1 pre-firing event correction weight (1-probability), down var. |
| **L1PreFiringWeight_ECAL_Nom** | Float_t| ECAL L1 pre-firing event correction weight (1-probability) |
| **L1PreFiringWeight_ECAL_Up** | Float_t| ECAL L1 pre-firing event correction weight (1-probability), up var. |
| **L1PreFiringWeight_Muon_Nom** | Float_t| Muon L1 pre-firing event correction weight (1-probability) |
| **L1PreFiringWeight_Muon_StatDn** | Float_t| Muon L1 pre-firing event correction weight (1-probability), down var. stat. |
| **L1PreFiringWeight_Muon_StatUp** | Float_t| Muon L1 pre-firing event correction weight (1-probability), up var. stat. |
| **L1PreFiringWeight_Muon_SystDn** | Float_t| Muon L1 pre-firing event correction weight (1-probability), down var. syst. |
| **L1PreFiringWeight_Muon_SystUp** | Float_t| Muon L1 pre-firing event correction weight (1-probability), up var. syst. |
| **L1PreFiringWeight_Nom** | Float_t| L1 pre-firing event correction weight (1-probability) |
| **L1PreFiringWeight_Up** | Float_t| L1 pre-firing event correction weight (1-probability), up var. |

## <a id='l1reco'></a>L1Reco [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **L1Reco_step** | Bool_t| Trigger/flag bit (process: RECO) |

## <a id='lowptelectron'></a>LowPtElectron [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **LowPtElectron_ID** | Float_t| New ID, BDT (raw) score |
| **LowPtElectron_charge** | Int_t| electric charge |
| **LowPtElectron_convVeto** | Bool_t| pass conversion veto |
| **LowPtElectron_convVtxRadius** | Float_t| conversion vertex radius (cm) |
| **LowPtElectron_convWP** | Int_t| conversion flag bit map: 1=Veto, 2=Loose, 3=Tight |
| **LowPtElectron_deltaEtaSC** | Float_t| delta eta (SC,ele) with sign |
| **LowPtElectron_dxy** | Float_t| dxy (with sign) wrt first PV, in cm |
| **LowPtElectron_dxyErr** | Float_t| dxy uncertainty, in cm |
| **LowPtElectron_dz** | Float_t| dz (with sign) wrt first PV, in cm |
| **LowPtElectron_dzErr** | Float_t| dz uncertainty, in cm |
| **LowPtElectron_eInvMinusPInv** | Float_t| 1/E_SC - 1/p_trk |
| **LowPtElectron_embeddedID** | Float_t| ID, BDT (raw) score |
| **LowPtElectron_energyErr** | Float_t| energy error of the cluster-track combination |
| **LowPtElectron_eta** | Float_t| eta |
| **LowPtElectron_hoe** | Float_t| H over E |
| **LowPtElectron_lostHits** | UChar_t| number of missing inner hits |
| **LowPtElectron_mass** | Float_t| mass |
| **LowPtElectron_miniPFRelIso_all** | Float_t| mini PF relative isolation, total (with scaled rho*EA PU corrections) |
| **LowPtElectron_miniPFRelIso_chg** | Float_t| mini PF relative isolation, charged component |
| **LowPtElectron_pdgId** | Int_t| PDG code assigned by the event reconstruction (not by MC truth) |
| **LowPtElectron_phi** | Float_t| phi |
| **LowPtElectron_pt** | Float_t| pt |
| **LowPtElectron_ptbiased** | Float_t| ElectronSeed, pT- and dxy- dependent BDT (raw) score |
| **LowPtElectron_r9** | Float_t| R9 of the SC, calculated with full 5x5 region |
| **LowPtElectron_scEtOverPt** | Float_t| (SC energy)/pt-1 |
| **LowPtElectron_sieie** | Float_t| sigma_IetaIeta of the SC, calculated with full 5x5 region |
| **LowPtElectron_unbiased** | Float_t| ElectronSeed, pT- and dxy- agnostic BDT (raw) score |
| **nLowPtElectron** | UInt_t| slimmedLowPtElectrons after basic selection (pt > 1. && userFloat('ID') > -0.25) |

## <a id='met'></a>MET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **MET_MetUnclustEnUpDeltaX** | Float_t| Delta (METx_mod-METx) Unclustered Energy Up |
| **MET_MetUnclustEnUpDeltaY** | Float_t| Delta (METy_mod-METy) Unclustered Energy Up |
| **MET_covXX** | Float_t| xx element of met covariance matrix |
| **MET_covXY** | Float_t| xy element of met covariance matrix |
| **MET_covYY** | Float_t| yy element of met covariance matrix |
| **MET_phi** | Float_t| phi |
| **MET_pt** | Float_t| pt |
| **MET_significance** | Float_t| MET significance |
| **MET_sumEt** | Float_t| scalar sum of Et |
| **MET_sumPtUnclustered** | Float_t| sumPt used for MET significance |

## <a id='muon'></a>Muon [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Muon_charge** | Int_t| electric charge |
| **Muon_cleanmask** | UChar_t| simple cleaning mask with priority to leptons |
| **Muon_dxy** | Float_t| dxy (with sign) wrt first PV, in cm |
| **Muon_dxyErr** | Float_t| dxy uncertainty, in cm |
| **Muon_dxybs** | Float_t| dxy (with sign) wrt the beam spot, in cm |
| **Muon_dz** | Float_t| dz (with sign) wrt first PV, in cm |
| **Muon_dzErr** | Float_t| dz uncertainty, in cm |
| **Muon_eta** | Float_t| eta |
| **Muon_fsrPhotonIdx** | Int_t(index to Fsrphoton)| Index of the associated FSR photon |
| **Muon_highPtId** | UChar_t| high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT) |
| **Muon_highPurity** | Bool_t| inner track is high purity |
| **Muon_inTimeMuon** | Bool_t| inTimeMuon ID |
| **Muon_ip3d** | Float_t| 3D impact parameter wrt first PV, in cm |
| **Muon_isGlobal** | Bool_t| muon is global muon |
| **Muon_isPFcand** | Bool_t| muon is PF candidate |
| **Muon_isStandalone** | Bool_t| muon is a standalone muon |
| **Muon_isTracker** | Bool_t| muon is tracker muon |
| **Muon_jetIdx** | Int_t(index to Jet)| index of the associated jet (-1 if none) |
| **Muon_jetNDauCharged** | UChar_t| number of charged daughters of the closest jet |
| **Muon_jetPtRelv2** | Float_t| Relative momentum of the lepton with respect to the closest jet after subtracting the lepton |
| **Muon_jetRelIso** | Float_t| Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet) |
| **Muon_looseId** | Bool_t| muon is loose muon |
| **Muon_mass** | Float_t| mass |
| **Muon_mediumId** | Bool_t| cut-based ID, medium WP |
| **Muon_mediumPromptId** | Bool_t| cut-based ID, medium prompt WP |
| **Muon_miniIsoId** | UChar_t| MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight) |
| **Muon_miniPFRelIso_all** | Float_t| mini PF relative isolation, total (with scaled rho*EA PU corrections) |
| **Muon_miniPFRelIso_chg** | Float_t| mini PF relative isolation, charged component |
| **Muon_multiIsoId** | UChar_t| MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium) |
| **Muon_mvaId** | UChar_t| Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight, 4=MvaVTight, 5=MvaVVTight) |
| **Muon_mvaLowPt** | Float_t| Low pt muon ID score |
| **Muon_mvaLowPtId** | UChar_t| Low Pt Mva ID from miniAOD selector (1=LowPtMvaLoose, 2=LowPtMvaMedium) |
| **Muon_mvaTTH** | Float_t| TTH MVA lepton ID score |
| **Muon_nStations** | Int_t| number of matched stations with default arbitration (segment & track) |
| **Muon_nTrackerLayers** | Int_t| number of layers in the tracker |
| **Muon_pdgId** | Int_t| PDG code assigned by the event reconstruction (not by MC truth) |
| **Muon_pfIsoId** | UChar_t| PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight) |
| **Muon_pfRelIso03_all** | Float_t| PF relative isolation dR=0.3, total (deltaBeta corrections) |
| **Muon_pfRelIso03_chg** | Float_t| PF relative isolation dR=0.3, charged component |
| **Muon_pfRelIso04_all** | Float_t| PF relative isolation dR=0.4, total (deltaBeta corrections) |
| **Muon_phi** | Float_t| phi |
| **Muon_pt** | Float_t| pt |
| **Muon_ptErr** | Float_t| ptError of the muon track |
| **Muon_puppiIsoId** | UChar_t| PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight) |
| **Muon_segmentComp** | Float_t| muon segment compatibility |
| **Muon_sip3d** | Float_t| 3D impact parameter significance wrt first PV |
| **Muon_softId** | Bool_t| soft cut-based ID |
| **Muon_softMva** | Float_t| soft MVA ID score |
| **Muon_softMvaId** | Bool_t| soft MVA ID |
| **Muon_tightCharge** | Int_t| Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass) |
| **Muon_tightId** | Bool_t| cut-based ID, tight WP |
| **Muon_tkIsoId** | UChar_t| TkIso ID (1=TkIsoLoose, 2=TkIsoTight) |
| **Muon_tkRelIso** | Float_t| Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt |
| **Muon_triggerIdLoose** | Bool_t| TriggerIdLoose ID |
| **Muon_tunepRelPt** | Float_t| TuneP relative pt, tunePpt/pt |
| **nMuon** | UInt_t| slimmedMuons after basic selection (pt > 15 || (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))) |

## <a id='otherpv'></a>OtherPV [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **OtherPV_z** | Float_t| Z position of other primary vertices, excluding the main PV |
| **nOtherPV** | UInt_t|  |

## <a id='ppslocaltrack'></a>PPSLocalTrack [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **PPSLocalTrack_decRPId** | Int_t| local track detector dec id |
| **PPSLocalTrack_multiRPProtonIdx** | Int_t(index to Multirpproton)| local track - proton correspondence |
| **PPSLocalTrack_rpType** | Int_t| strip=3, pixel=4, diamond=5, timing=6 |
| **PPSLocalTrack_singleRPProtonIdx** | Int_t(index to Singlerpproton)| local track - proton correspondence |
| **PPSLocalTrack_time** | Float_t| local track time |
| **PPSLocalTrack_timeUnc** | Float_t| local track time uncertainty |
| **PPSLocalTrack_x** | Float_t| local track x |
| **PPSLocalTrack_y** | Float_t| local track y |
| **nPPSLocalTrack** | UInt_t| ppsLocalTrack variables |

## <a id='pv'></a>PV [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **PV_chi2** | Float_t| main primary vertex reduced chi2 |
| **PV_ndof** | Float_t| main primary vertex number of degree of freedom |
| **PV_npvs** | Int_t| total number of reconstructed primary vertices |
| **PV_npvsGood** | Int_t| number of good reconstructed primary vertices. selection:!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2 |
| **PV_score** | Float_t| main primary vertex score, i.e. sum pt2 of clustered objects |
| **PV_x** | Float_t| main primary vertex position x coordinate |
| **PV_y** | Float_t| main primary vertex position y coordinate |
| **PV_z** | Float_t| main primary vertex position z coordinate |

## <a id='photon'></a>Photon [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Photon_charge** | Int_t| electric charge |
| **Photon_cleanmask** | UChar_t| simple cleaning mask with priority to leptons |
| **Photon_cutBased** | Int_t| cut-based ID bitmap, Fall17V2, (0:fail, 1:loose, 2:medium, 3:tight) |
| **Photon_cutBased_Fall17V1Bitmap** | Int_t| cut-based ID bitmap, Fall17V1, 2^(0:loose, 1:medium, 2:tight). |
| **Photon_dEscaleDown** | Float_t| ecal energy scale shifted 1 sigma down (adding gain/stat/syst in quadrature) |
| **Photon_dEscaleUp** | Float_t| ecal energy scale shifted 1 sigma up (adding gain/stat/syst in quadrature) |
| **Photon_dEsigmaDown** | Float_t| ecal energy smearing value shifted 1 sigma up |
| **Photon_dEsigmaUp** | Float_t| ecal energy smearing value shifted 1 sigma up |
| **Photon_eCorr** | Float_t| ratio of the calibrated energy/miniaod energy |
| **Photon_electronIdx** | Int_t(index to Electron)| index of the associated electron (-1 if none) |
| **Photon_electronVeto** | Bool_t| pass electron veto |
| **Photon_energyErr** | Float_t| energy error of the cluster from regression |
| **Photon_eta** | Float_t| eta |
| **Photon_hoe** | Float_t| H over E |
| **Photon_isScEtaEB** | Bool_t| is supercluster eta within barrel acceptance |
| **Photon_isScEtaEE** | Bool_t| is supercluster eta within endcap acceptance |
| **Photon_jetIdx** | Int_t(index to Jet)| index of the associated jet (-1 if none) |
| **Photon_mass** | Float_t| mass |
| **Photon_mvaID** | Float_t| MVA ID score, Fall17V2 |
| **Photon_mvaID_Fall17V1p1** | Float_t| MVA ID score, Fall17V1p1 |
| **Photon_mvaID_WP80** | Bool_t| MVA ID WP80, Fall17V2 |
| **Photon_mvaID_WP90** | Bool_t| MVA ID WP90, Fall17V2 |
| **Photon_pdgId** | Int_t| PDG code assigned by the event reconstruction (not by MC truth) |
| **Photon_pfRelIso03_all** | Float_t| PF relative isolation dR=0.3, total (with rho*EA PU corrections) |
| **Photon_pfRelIso03_chg** | Float_t| PF relative isolation dR=0.3, charged component (with rho*EA PU corrections) |
| **Photon_phi** | Float_t| phi |
| **Photon_pixelSeed** | Bool_t| has pixel seed |
| **Photon_pt** | Float_t| p_{T} |
| **Photon_r9** | Float_t| R9 of the supercluster, calculated with full 5x5 region |
| **Photon_seedGain** | UChar_t| Gain of the seed crystal |
| **Photon_sieie** | Float_t| sigma_IetaIeta of the supercluster, calculated with full 5x5 region |
| **Photon_vidNestedWPBitmap** | Int_t| Fall17V2 VID compressed bitmap (MinPtCut,PhoSCEtaMultiRangeCut,PhoSingleTowerHadOverEmCut,PhoFull5x5SigmaIEtaIEtaCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut), 2 bits per cut |
| **nPhoton** | UInt_t| slimmedPhotons after basic selection (pt > 5 ) |

## <a id='proton'></a>Proton [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Proton_multiRP_arm** | Int_t| 0 = sector45, 1 = sector56 |
| **Proton_multiRP_t** | Float_t| Mandelstam variable t |
| **Proton_multiRP_thetaX** | Float_t| theta x |
| **Proton_multiRP_thetaY** | Float_t| theta y |
| **Proton_multiRP_time** | Float_t| time |
| **Proton_multiRP_timeUnc** | Float_t| time uncertainty |
| **Proton_multiRP_xi** | Float_t| xi or dp/p |
| **Proton_singleRP_decRPId** | Int_t| Detector ID |
| **Proton_singleRP_thetaY** | Float_t| th y |
| **Proton_singleRP_xi** | Float_t| xi or dp/p |
| **nProton_multiRP** | UInt_t| bon |
| **nProton_singleRP** | UInt_t| bon |

## <a id='puppimet'></a>PuppiMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **PuppiMET_phi** | Float_t| phi |
| **PuppiMET_phiJERDown** | Float_t| JER down phi |
| **PuppiMET_phiJERUp** | Float_t| JER up phi |
| **PuppiMET_phiJESDown** | Float_t| JES down phi |
| **PuppiMET_phiJESUp** | Float_t| JES up phi |
| **PuppiMET_phiUnclusteredDown** | Float_t| Unclustered down phi |
| **PuppiMET_phiUnclusteredUp** | Float_t| Unclustered up phi |
| **PuppiMET_pt** | Float_t| pt |
| **PuppiMET_ptJERDown** | Float_t| JER down pt |
| **PuppiMET_ptJERUp** | Float_t| JER up pt |
| **PuppiMET_ptJESDown** | Float_t| JES down pt |
| **PuppiMET_ptJESUp** | Float_t| JES up pt |
| **PuppiMET_ptUnclusteredDown** | Float_t| Unclustered down pt |
| **PuppiMET_ptUnclusteredUp** | Float_t| Unclustered up pt |
| **PuppiMET_sumEt** | Float_t| scalar sum of Et |

## <a id='rawmet'></a>RawMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **RawMET_phi** | Float_t| phi |
| **RawMET_pt** | Float_t| pt |
| **RawMET_sumEt** | Float_t| scalar sum of Et |

## <a id='rawpuppimet'></a>RawPuppiMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **RawPuppiMET_phi** | Float_t| phi |
| **RawPuppiMET_pt** | Float_t| pt |
| **RawPuppiMET_sumEt** | Float_t| scalar sum of Et |

## <a id='sv'></a>SV [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SV_charge** | Int_t| sum of the charge of the SV tracks |
| **SV_chi2** | Float_t| reduced chi2, i.e. chi/ndof |
| **SV_dlen** | Float_t| decay length in cm |
| **SV_dlenSig** | Float_t| decay length significance |
| **SV_dxy** | Float_t| 2D decay length in cm |
| **SV_dxySig** | Float_t| 2D decay length significance |
| **SV_eta** | Float_t| eta |
| **SV_mass** | Float_t| mass |
| **SV_ndof** | Float_t| number of degrees of freedom |
| **SV_ntracks** | UChar_t| number of tracks |
| **SV_pAngle** | Float_t| pointing angle, i.e. acos(p_SV * (SV - PV))  |
| **SV_phi** | Float_t| phi |
| **SV_pt** | Float_t| pt |
| **SV_x** | Float_t| secondary vertex X position, in cm |
| **SV_y** | Float_t| secondary vertex Y position, in cm |
| **SV_z** | Float_t| secondary vertex Z position, in cm |
| **nSV** | UInt_t|  |

## <a id='softactivityjet'></a>SoftActivityJet [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJet_eta** | Float_t| eta |
| **SoftActivityJet_phi** | Float_t| phi |
| **SoftActivityJet_pt** | Float_t| pt |
| **nSoftActivityJet** | UInt_t| jets clustered from charged candidates compatible with primary vertex (charge()!=0 && pvAssociationQuality()>=5 && vertexRef().key()==0) |

## <a id='softactivityjetht'></a>SoftActivityJetHT [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetHT** | Float_t| scalar sum of soft activity jet pt, pt>1 |

## <a id='softactivityjetht10'></a>SoftActivityJetHT10 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetHT10** | Float_t| scalar sum of soft activity jet pt , pt >10 |

## <a id='softactivityjetht2'></a>SoftActivityJetHT2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetHT2** | Float_t| scalar sum of soft activity jet pt, pt >2 |

## <a id='softactivityjetht5'></a>SoftActivityJetHT5 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetHT5** | Float_t| scalar sum of soft activity jet pt, pt>5 |

## <a id='softactivityjetnjets10'></a>SoftActivityJetNjets10 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetNjets10** | Int_t| number of soft activity jet pt, pt >2 |

## <a id='softactivityjetnjets2'></a>SoftActivityJetNjets2 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetNjets2** | Int_t| number of soft activity jet pt, pt >10 |

## <a id='softactivityjetnjets5'></a>SoftActivityJetNjets5 [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SoftActivityJetNjets5** | Int_t| number of soft activity jet pt, pt >5 |

## <a id='subjet'></a>SubJet [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **SubJet_btagCSVV2** | Float_t|  pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2) |
| **SubJet_btagDeepB** | Float_t| DeepCSV b+bb tag discriminator |
| **SubJet_eta** | Float_t| eta |
| **SubJet_mass** | Float_t| mass |
| **SubJet_n2b1** | Float_t| N2 with beta=1 |
| **SubJet_n3b1** | Float_t| N3 with beta=1 |
| **SubJet_phi** | Float_t| phi |
| **SubJet_pt** | Float_t| pt |
| **SubJet_rawFactor** | Float_t| 1 - Factor to get back to raw pT |
| **SubJet_tau1** | Float_t| Nsubjettiness (1 axis) |
| **SubJet_tau2** | Float_t| Nsubjettiness (2 axis) |
| **SubJet_tau3** | Float_t| Nsubjettiness (3 axis) |
| **SubJet_tau4** | Float_t| Nsubjettiness (4 axis) |
| **nSubJet** | UInt_t| slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis |

## <a id='tau'></a>Tau [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **Tau_charge** | Int_t| electric charge |
| **Tau_chargedIso** | Float_t| charged isolation |
| **Tau_cleanmask** | UChar_t| simple cleaning mask with priority to leptons |
| **Tau_decayMode** | Int_t| decayMode() |
| **Tau_dxy** | Float_t| d_{xy} of lead track with respect to PV, in cm (with sign) |
| **Tau_dz** | Float_t| d_{z} of lead track with respect to PV, in cm (with sign) |
| **Tau_eta** | Float_t| eta |
| **Tau_idAntiEleDeadECal** | Bool_t| Anti-electron dead-ECal discriminator |
| **Tau_idAntiMu** | UChar_t| Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight |
| **Tau_idDecayModeOldDMs** | Bool_t| tauID('decayModeFinding') |
| **Tau_idDeepTau2017v2p1VSe** | UChar_t| byDeepTau2017v2p1VSe ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight |
| **Tau_idDeepTau2017v2p1VSjet** | UChar_t| byDeepTau2017v2p1VSjet ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight |
| **Tau_idDeepTau2017v2p1VSmu** | UChar_t| byDeepTau2017v2p1VSmu ID working points (deepTau2017v2p1): bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight |
| **Tau_jetIdx** | Int_t(index to Jet)| index of the associated jet (-1 if none) |
| **Tau_leadTkDeltaEta** | Float_t| eta of the leading track, minus tau eta |
| **Tau_leadTkDeltaPhi** | Float_t| phi of the leading track, minus tau phi |
| **Tau_leadTkPtOverTauPt** | Float_t| pt of the leading track divided by tau pt |
| **Tau_mass** | Float_t| mass |
| **Tau_neutralIso** | Float_t| neutral (photon) isolation |
| **Tau_phi** | Float_t| phi |
| **Tau_photonsOutsideSignalCone** | Float_t| sum of photons outside signal cone |
| **Tau_pt** | Float_t| pt |
| **Tau_puCorr** | Float_t| pileup correction |
| **Tau_rawDeepTau2017v2p1VSe** | Float_t| byDeepTau2017v2p1VSe raw output discriminator (deepTau2017v2p1) |
| **Tau_rawDeepTau2017v2p1VSjet** | Float_t| byDeepTau2017v2p1VSjet raw output discriminator (deepTau2017v2p1) |
| **Tau_rawDeepTau2017v2p1VSmu** | Float_t| byDeepTau2017v2p1VSmu raw output discriminator (deepTau2017v2p1) |
| **Tau_rawIso** | Float_t| combined isolation (deltaBeta corrections) |
| **Tau_rawIsodR03** | Float_t| combined isolation (deltaBeta corrections, dR=0.3) |
| **nTau** | UInt_t| slimmedTaus after basic selection (pt > 18 && tauID('decayModeFindingNewDMs') && (tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') || (tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum'))<2.5) || tauID('byVVVLooseDeepTau2017v2p1VSjet'))) |

## <a id='tkmet'></a>TkMET [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **TkMET_phi** | Float_t| raw track MET phi |
| **TkMET_pt** | Float_t| raw track MET pt |
| **TkMET_sumEt** | Float_t| raw track scalar sum of Et |

## <a id='trigobj'></a>TrigObj [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **TrigObj_eta** | Float_t| eta |
| **TrigObj_filterBits** | Int_t| extra bits of associated information: 1 = CaloIdL_TrackIdL_IsoVL, 2 = 1e (WPTight), 4 = 1e (WPLoose), 8 = OverlapFilter PFTau, 16 = 2e, 32 = 1e-1mu, 64 = 1e-1tau, 128 = 3e, 256 = 2e-1mu, 512 = 1e-2mu, 1024 = 1e (32_L1DoubleEG_AND_L1SingleEGOr), 2048 = 1e (CaloIdVT_GsfTrkIdT), 4096 = 1e (PFJet), 8192 = 1e (Photon175_OR_Photon200) for Electron (PixelMatched e/gamma); 1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = 1mu, 16 = 2mu, 32 = 1mu-1e, 64 = 1mu-1tau, 128 = 3mu, 256 = 2mu-1e, 512 = 1mu-2e, 1024 = 1mu (Mu50), 2048 = 1mu (Mu100) for Muon; 1 = LooseChargedIso, 2 = MediumChargedIso, 4 = TightChargedIso, 8 = TightID OOSC photons, 16 = HPS, 32 = single-tau + tau+MET, 64 = di-tau, 128 = e-tau, 256 = mu-tau, 512 = VBF+di-tau for Tau; Jet bits: bit 0 for VBF cross-cleaned from loose iso PFTau, bit 1 for hltBTagCaloCSVp087Triple, bit 2 for hltDoubleCentralJet90, bit 3 for hltDoublePFCentralJetLooseID90, bit 4 for hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet, bit 5 for hltQuadCentralJet30, bit 6 for hltQuadPFCentralJetLooseID30, bit 7 for hltL1sQuadJetC50IorQuadJetC60IorHTT280IorHTT300IorHTT320IorTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBF or hltL1sQuadJetCIorTripleJetVBFIorHTT, bit 8 for hltQuadCentralJet45, bit 9 for hltQuadPFCentralJetLooseID45, bit 10  for hltL1sQuadJetC60IorHTT380IorHTT280QuadJetIorHTT300QuadJet or hltL1sQuadJetC50to60IorHTT280to500IorHTT250to340QuadJet bit 11 for hltBTagCaloCSVp05Double or hltBTagCaloDeepCSVp17Double, bit 12 for hltPFCentralJetLooseIDQuad30, bit 13 for hlt1PFCentralJetLooseID75, bit 14 for hlt2PFCentralJetLooseID60, bit 15 for hlt3PFCentralJetLooseID45, bit 16 for hlt4PFCentralJetLooseID40, bit 17 for hltBTagPFCSVp070Triple or hltBTagPFDeepCSVp24Triple or hltBTagPFDeepCSV4p5Triple  for Jet; HT bits: bit 0 for hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet, bit 1 for hltL1sQuadJetC50IorQuadJetC60IorHTT280IorHTT300IorHTT320IorTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBF or hltL1sQuadJetCIorTripleJetVBFIorHTT, bit 2 for hltL1sQuadJetC60IorHTT380IorHTT280QuadJetIorHTT300QuadJet or hltL1sQuadJetC50to60IorHTT280to500IorHTT250to340QuadJet, bit 3 for hltCaloQuadJet30HT300 or hltCaloQuadJet30HT320, bit 4 for hltPFCentralJetsLooseIDQuad30HT300 or hltPFCentralJetsLooseIDQuad30HT330 for HT; MHT bits: bit 0 for hltCaloQuadJet30HT300 or hltCaloQuadJet30HT320, bit 1 for hltPFCentralJetsLooseIDQuad30HT300 or hltPFCentralJetsLooseIDQuad30HT330 for MHT;  |
| **TrigObj_id** | Int_t| ID of the object: 11 = Electron (PixelMatched e/gamma), 22 = Photon (PixelMatch-vetoed e/gamma), 13 = Muon, 15 = Tau, 1 = Jet, 6 = FatJet, 2 = MET, 3 = HT, 4 = MHT |
| **TrigObj_l1charge** | Int_t| charge of associated L1 seed |
| **TrigObj_l1iso** | Int_t| iso of associated L1 seed |
| **TrigObj_l1pt** | Float_t| pt of associated L1 seed |
| **TrigObj_l1pt_2** | Float_t| pt of associated secondary L1 seed |
| **TrigObj_l2pt** | Float_t| pt of associated 'L2' seed (i.e. HLT before tracking/PF) |
| **TrigObj_phi** | Float_t| phi |
| **TrigObj_pt** | Float_t| pt |
| **nTrigObj** | UInt_t|  |

## <a id='boostedtau'></a>boostedTau [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **boostedTau_charge** | Int_t| electric charge |
| **boostedTau_chargedIso** | Float_t| charged isolation |
| **boostedTau_decayMode** | Int_t| decayMode() |
| **boostedTau_eta** | Float_t| eta |
| **boostedTau_idAntiEle2018** | UChar_t| Anti-electron MVA discriminator V6 (2018): bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight |
| **boostedTau_idAntiMu** | UChar_t| Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight |
| **boostedTau_idMVAnewDM2017v2** | UChar_t| IsolationMVArun2017v2DBnewDMwLT ID working point (2017v2): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight |
| **boostedTau_idMVAoldDM2017v2** | UChar_t| IsolationMVArun2017v2DBoldDMwLT ID working point (2017v2): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight |
| **boostedTau_idMVAoldDMdR032017v2** | UChar_t| IsolationMVArun2017v2DBoldDMdR0p3wLT ID working point (2017v2): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight |
| **boostedTau_jetIdx** | Int_t(index to Jet)| index of the associated jet (-1 if none) |
| **boostedTau_leadTkDeltaEta** | Float_t| eta of the leading track, minus tau eta |
| **boostedTau_leadTkDeltaPhi** | Float_t| phi of the leading track, minus tau phi |
| **boostedTau_leadTkPtOverTauPt** | Float_t| pt of the leading track divided by tau pt |
| **boostedTau_mass** | Float_t| mass |
| **boostedTau_neutralIso** | Float_t| neutral (photon) isolation |
| **boostedTau_phi** | Float_t| phi |
| **boostedTau_photonsOutsideSignalCone** | Float_t| sum of photons outside signal cone |
| **boostedTau_pt** | Float_t| pt |
| **boostedTau_puCorr** | Float_t| pileup correction |
| **boostedTau_rawAntiEle2018** | Float_t| Anti-electron MVA discriminator V6 raw output discriminator (2018) |
| **boostedTau_rawAntiEleCat2018** | Int_t| Anti-electron MVA discriminator V6 category (2018) |
| **boostedTau_rawIso** | Float_t| combined isolation (deltaBeta corrections) |
| **boostedTau_rawIsodR03** | Float_t| combined isolation (deltaBeta corrections, dR=0.3) |
| **boostedTau_rawMVAnewDM2017v2** | Float_t| byIsolationMVArun2017v2DBnewDMwLT raw output discriminator (2017v2) |
| **boostedTau_rawMVAoldDM2017v2** | Float_t| byIsolationMVArun2017v2DBoldDMwLT raw output discriminator (2017v2) |
| **boostedTau_rawMVAoldDMdR032017v2** | Float_t| byIsolationMVArun2017v2DBoldDMdR0p3wLT raw output discriminator (2017v2) |
| **nboostedTau** | UInt_t| slimmedBoostedTaus after basic selection (pt > 40 && tauID('decayModeFindingNewDMs') && (tauID('byVVLooseIsolationMVArun2017v2DBoldDMwLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBnewDMwLT2017'))) |

## <a id='event'></a>event [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **event** | ULong64_t| event/l |

## <a id='fixedgridrhofastjetall'></a>fixedGridRhoFastjetAll [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **fixedGridRhoFastjetAll** | Float_t| rho from all PF Candidates, used e.g. for JECs |

## <a id='fixedgridrhofastjetcentral'></a>fixedGridRhoFastjetCentral [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **fixedGridRhoFastjetCentral** | Float_t| rho from all PF Candidates for central region, used e.g. for JECs |

## <a id='fixedgridrhofastjetcentralcalo'></a>fixedGridRhoFastjetCentralCalo [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **fixedGridRhoFastjetCentralCalo** | Float_t| rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation |

## <a id='fixedgridrhofastjetcentralchargedpileup'></a>fixedGridRhoFastjetCentralChargedPileUp [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **fixedGridRhoFastjetCentralChargedPileUp** | Float_t| rho from charged PF Candidates for central region, used e.g. for JECs |

## <a id='fixedgridrhofastjetcentralneutral'></a>fixedGridRhoFastjetCentralNeutral [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **fixedGridRhoFastjetCentralNeutral** | Float_t| rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations |

## <a id='luminosityblock'></a>luminosityBlock [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **luminosityBlock** | UInt_t| luminosityBlock/i |

## <a id='run'></a>run [<sup>[back to top]</sup>](#content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

