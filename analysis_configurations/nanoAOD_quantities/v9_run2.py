from code_generation.quantity import NanoAODQuantity

CaloMET_phi = NanoAODQuantity("CaloMET_phi")                                                                                                                                              
"""dtype: Float_t; description: phi"""
CaloMET_pt = NanoAODQuantity("CaloMET_pt")                                                                                                                                                
"""dtype: Float_t; description: pt"""
CaloMET_sumEt = NanoAODQuantity("CaloMET_sumEt")                                                                                                                                          
"""dtype: Float_t; description: scalar sum of Et"""

ChsMET_phi = NanoAODQuantity("ChsMET_phi")                                                                                                                                                
"""dtype: Float_t; description: raw chs PF MET phi"""
ChsMET_pt = NanoAODQuantity("ChsMET_pt")                                                                                                                                                  
"""dtype: Float_t; description: raw chs PF MET pt"""
ChsMET_sumEt = NanoAODQuantity("ChsMET_sumEt")                                                                                                                                            
"""dtype: Float_t; description: raw chs PF scalar sum of Et"""

nCorrT1METJet = NanoAODQuantity("nCorrT1METJet")                                                                                                                                          
"""dtype: UInt_t; description: Additional low-pt jets for Type-1 MET re-correction"""
CorrT1METJet_area = NanoAODQuantity("CorrT1METJet_area")                                                                                                                                  
"""dtype: Float_t; description: jet catchment area, for JECs"""
CorrT1METJet_eta = NanoAODQuantity("CorrT1METJet_eta")                                                                                                                                    
"""dtype: Float_t; description: eta"""
CorrT1METJet_muonSubtrFactor = NanoAODQuantity("CorrT1METJet_muonSubtrFactor")                                                                                                            
"""dtype: Float_t; description: 1-(muon-subtracted raw pt)/(raw pt)"""
CorrT1METJet_phi = NanoAODQuantity("CorrT1METJet_phi")                                                                                                                                    
"""dtype: Float_t; description: phi"""
CorrT1METJet_rawPt = NanoAODQuantity("CorrT1METJet_rawPt")                                                                                                                                
"""dtype: Float_t; description: pt()*jecFactor('Uncorrected')"""

DeepMETResolutionTune_phi = NanoAODQuantity("DeepMETResolutionTune_phi")                                                                                                                  
"""dtype: Float_t; description: DeepmET ResolutionTune phi"""
DeepMETResolutionTune_pt = NanoAODQuantity("DeepMETResolutionTune_pt")                                                                                                                    
"""dtype: Float_t; description: DeepMET ResolutionTune pt"""

DeepMETResponseTune_phi = NanoAODQuantity("DeepMETResponseTune_phi")                                                                                                                      
"""dtype: Float_t; description: DeepMET ResponseTune phi"""
DeepMETResponseTune_pt = NanoAODQuantity("DeepMETResponseTune_pt")                                                                                                                        
"""dtype: Float_t; description: DeepMET ResponseTune pt"""

nElectron = NanoAODQuantity("nElectron")                                                                                                                                                  
"""dtype: UInt_t; description: slimmedElectrons after basic selection (pt > 5 )"""
Electron_charge = NanoAODQuantity("Electron_charge")                                                                                                                                      
"""dtype: Int_t; description: electric charge"""
Electron_cleanmask = NanoAODQuantity("Electron_cleanmask")                                                                                                                                
"""dtype: UChar_t; description: simple cleaning mask with priority to leptons"""
Electron_convVeto = NanoAODQuantity("Electron_convVeto")                                                                                                                                  
"""dtype: Bool_t; description: pass conversion veto"""
Electron_cutBased = NanoAODQuantity("Electron_cutBased")                                                                                                                                  
"""dtype: Int_t; description: cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)"""
Electron_cutBased_HEEP = NanoAODQuantity("Electron_cutBased_HEEP")                                                                                                                        
"""dtype: Bool_t; description: cut-based HEEP ID"""
Electron_dEscaleDown = NanoAODQuantity("Electron_dEscaleDown")                                                                                                                            
"""dtype: Float_t; description: ecal energy scale shifted 1 sigma down (adding gain/stat/syst in quadrature)"""
Electron_dEscaleUp = NanoAODQuantity("Electron_dEscaleUp")                                                                                                                                
"""dtype: Float_t; description: ecal energy scale shifted 1 sigma up(adding gain/stat/syst in quadrature)"""
Electron_dEsigmaDown = NanoAODQuantity("Electron_dEsigmaDown")                                                                                                                            
"""dtype: Float_t; description: ecal energy smearing value shifted 1 sigma up"""
Electron_dEsigmaUp = NanoAODQuantity("Electron_dEsigmaUp")                                                                                                                                
"""dtype: Float_t; description: ecal energy smearing value shifted 1 sigma up"""
Electron_deltaEtaSC = NanoAODQuantity("Electron_deltaEtaSC")                                                                                                                              
"""dtype: Float_t; description: delta eta (SC,ele) with sign"""
Electron_dr03EcalRecHitSumEt = NanoAODQuantity("Electron_dr03EcalRecHitSumEt")                                                                                                            
"""dtype: Float_t; description: Non-PF Ecal isolation within a delta R cone of 0.3 with electron pt > 35 GeV"""
Electron_dr03HcalDepth1TowerSumEt = NanoAODQuantity("Electron_dr03HcalDepth1TowerSumEt")                                                                                                  
"""dtype: Float_t; description: Non-PF Hcal isolation within a delta R cone of 0.3 with electron pt > 35 GeV"""
Electron_dr03TkSumPt = NanoAODQuantity("Electron_dr03TkSumPt")                                                                                                                            
"""dtype: Float_t; description: Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV"""
Electron_dr03TkSumPtHEEP = NanoAODQuantity("Electron_dr03TkSumPtHEEP")                                                                                                                    
"""dtype: Float_t; description: Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV used in HEEP ID"""
Electron_dxy = NanoAODQuantity("Electron_dxy")                                                                                                                                            
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm"""
Electron_dxyErr = NanoAODQuantity("Electron_dxyErr")                                                                                                                                      
"""dtype: Float_t; description: dxy uncertainty, in cm"""
Electron_dz = NanoAODQuantity("Electron_dz")                                                                                                                                              
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm"""
Electron_dzErr = NanoAODQuantity("Electron_dzErr")                                                                                                                                        
"""dtype: Float_t; description: dz uncertainty, in cm"""
Electron_eCorr = NanoAODQuantity("Electron_eCorr")                                                                                                                                        
"""dtype: Float_t; description: ratio of the calibrated energy/miniaod energy"""
Electron_eInvMinusPInv = NanoAODQuantity("Electron_eInvMinusPInv")                                                                                                                        
"""dtype: Float_t; description: 1/E_SC - 1/p_trk"""
Electron_energyErr = NanoAODQuantity("Electron_energyErr")                                                                                                                                
"""dtype: Float_t; description: energy error of the cluster-track combination"""
Electron_eta = NanoAODQuantity("Electron_eta")                                                                                                                                            
"""dtype: Float_t; description: eta"""
Electron_genPartFlav = NanoAODQuantity("Electron_genPartFlav")                                                                                                                            
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 electrons or photons: 1 = prompt electron (including gamma*->mu mu), 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched"""
Electron_genPartIdx = NanoAODQuantity("Electron_genPartIdx")                                                                                                                              
"""dtype: Int_t; description: Index into genParticle list for MC matching to status==1 electrons or photons"""
Electron_hoe = NanoAODQuantity("Electron_hoe")                                                                                                                                            
"""dtype: Float_t; description: H over E"""
Electron_ip3d = NanoAODQuantity("Electron_ip3d")                                                                                                                                          
"""dtype: Float_t; description: 3D impact parameter wrt first PV, in cm"""
Electron_isPFcand = NanoAODQuantity("Electron_isPFcand")                                                                                                                                  
"""dtype: Bool_t; description: electron is PF candidate"""
Electron_jetIdx = NanoAODQuantity("Electron_jetIdx")                                                                                                                                      
"""dtype: Int_t; description: index of the associated jet (-1 if none)"""
Electron_jetNDauCharged = NanoAODQuantity("Electron_jetNDauCharged")                                                                                                                      
"""dtype: UChar_t; description: number of charged daughters of the closest jet"""
Electron_jetPtRelv2 = NanoAODQuantity("Electron_jetPtRelv2")                                                                                                                              
"""dtype: Float_t; description: Relative momentum of the lepton with respect to the closest jet after subtracting the lepton"""
Electron_jetRelIso = NanoAODQuantity("Electron_jetRelIso")                                                                                                                                
"""dtype: Float_t; description: Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)"""
Electron_lostHits = NanoAODQuantity("Electron_lostHits")                                                                                                                                  
"""dtype: UChar_t; description: number of missing inner hits"""
Electron_mass = NanoAODQuantity("Electron_mass")                                                                                                                                          
"""dtype: Float_t; description: mass"""
Electron_miniPFRelIso_all = NanoAODQuantity("Electron_miniPFRelIso_all")                                                                                                                  
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections)"""
Electron_miniPFRelIso_chg = NanoAODQuantity("Electron_miniPFRelIso_chg")                                                                                                                  
"""dtype: Float_t; description: mini PF relative isolation, charged component"""
Electron_mvaFall17V2Iso = NanoAODQuantity("Electron_mvaFall17V2Iso")                                                                                                                      
"""dtype: Float_t; description: MVA Iso ID V2 score"""
Electron_mvaFall17V2Iso_WP80 = NanoAODQuantity("Electron_mvaFall17V2Iso_WP80")                                                                                                            
"""dtype: Bool_t; description: MVA Iso ID V2 WP80"""
Electron_mvaFall17V2Iso_WP90 = NanoAODQuantity("Electron_mvaFall17V2Iso_WP90")                                                                                                            
"""dtype: Bool_t; description: MVA Iso ID V2 WP90"""
Electron_mvaFall17V2Iso_WPL = NanoAODQuantity("Electron_mvaFall17V2Iso_WPL")                                                                                                              
"""dtype: Bool_t; description: MVA Iso ID V2 loose WP"""
Electron_mvaFall17V2noIso = NanoAODQuantity("Electron_mvaFall17V2noIso")                                                                                                                  
"""dtype: Float_t; description: MVA noIso ID V2 score"""
Electron_mvaFall17V2noIso_WP80 = NanoAODQuantity("Electron_mvaFall17V2noIso_WP80")                                                                                                        
"""dtype: Bool_t; description: MVA noIso ID V2 WP80"""
Electron_mvaFall17V2noIso_WP90 = NanoAODQuantity("Electron_mvaFall17V2noIso_WP90")                                                                                                        
"""dtype: Bool_t; description: MVA noIso ID V2 WP90"""
Electron_mvaFall17V2noIso_WPL = NanoAODQuantity("Electron_mvaFall17V2noIso_WPL")                                                                                                          
"""dtype: Bool_t; description: MVA noIso ID V2 loose WP"""
Electron_mvaTTH = NanoAODQuantity("Electron_mvaTTH")                                                                                                                                      
"""dtype: Float_t; description: TTH MVA lepton ID score"""
Electron_pdgId = NanoAODQuantity("Electron_pdgId")                                                                                                                                        
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth)"""
Electron_pfRelIso03_all = NanoAODQuantity("Electron_pfRelIso03_all")                                                                                                                      
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (with rho*EA PU corrections)"""
Electron_pfRelIso03_chg = NanoAODQuantity("Electron_pfRelIso03_chg")                                                                                                                      
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component"""
Electron_phi = NanoAODQuantity("Electron_phi")                                                                                                                                            
"""dtype: Float_t; description: phi"""
Electron_photonIdx = NanoAODQuantity("Electron_photonIdx")                                                                                                                                
"""dtype: Int_t; description: index of the associated photon (-1 if none)"""
Electron_pt = NanoAODQuantity("Electron_pt")                                                                                                                                              
"""dtype: Float_t; description: p_{T}"""
Electron_r9 = NanoAODQuantity("Electron_r9")                                                                                                                                              
"""dtype: Float_t; description: R9 of the supercluster, calculated with full 5x5 region"""
Electron_scEtOverPt = NanoAODQuantity("Electron_scEtOverPt")                                                                                                                              
"""dtype: Float_t; description: (supercluster transverse energy)/pt-1"""
Electron_seedGain = NanoAODQuantity("Electron_seedGain")                                                                                                                                  
"""dtype: UChar_t; description: Gain of the seed crystal"""
Electron_sieie = NanoAODQuantity("Electron_sieie")                                                                                                                                        
"""dtype: Float_t; description: sigma_IetaIeta of the supercluster, calculated with full 5x5 region"""
Electron_sip3d = NanoAODQuantity("Electron_sip3d")                                                                                                                                        
"""dtype: Float_t; description: 3D impact parameter significance wrt first PV, in cm"""
Electron_tightCharge = NanoAODQuantity("Electron_tightCharge")                                                                                                                            
"""dtype: Int_t; description: Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent)"""
Electron_vidNestedWPBitmap = NanoAODQuantity("Electron_vidNestedWPBitmap")                                                                                                                
"""dtype: Int_t; description: VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaCut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEInverseMinusPInverseCut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 3 bits per cut"""
Electron_vidNestedWPBitmapHEEP = NanoAODQuantity("Electron_vidNestedWPBitmapHEEP")                                                                                                        
"""dtype: Int_t; description: VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaWithSatCut,GsfEleFull5x5E2x5OverE5x5WithSatCut,GsfEleHadronicOverEMLinearCut,GsfEleTrkPtIsoCut,GsfEleEmHadD1IsoRhoCut,GsfEleDxyCut,GsfEleMissingHitsCut,GsfEleEcalDrivenCut), 1 bits per cut"""

nFatJet = NanoAODQuantity("nFatJet")                                                                                                                                                      
"""dtype: UInt_t; description: slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"""
FatJet_area = NanoAODQuantity("FatJet_area")                                                                                                                                              
"""dtype: Float_t; description: jet catchment area, for JECs"""
FatJet_btagCSVV2 = NanoAODQuantity("FatJet_btagCSVV2")                                                                                                                                    
"""dtype: Float_t; description:  pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)"""
FatJet_btagDDBvLV2 = NanoAODQuantity("FatJet_btagDDBvLV2")                                                                                                                                
"""dtype: Float_t; description: DeepDoubleX V2(mass-decorrelated) discriminator for H(Z)->bb vs QCD"""
FatJet_btagDDCvBV2 = NanoAODQuantity("FatJet_btagDDCvBV2")                                                                                                                                
"""dtype: Float_t; description: DeepDoubleX V2 (mass-decorrelated) discriminator for H(Z)->cc vs H(Z)->bb"""
FatJet_btagDDCvLV2 = NanoAODQuantity("FatJet_btagDDCvLV2")                                                                                                                                
"""dtype: Float_t; description: DeepDoubleX V2 (mass-decorrelated) discriminator for H(Z)->cc vs QCD"""
FatJet_btagDeepB = NanoAODQuantity("FatJet_btagDeepB")                                                                                                                                    
"""dtype: Float_t; description: DeepCSV b+bb tag discriminator"""
FatJet_btagHbb = NanoAODQuantity("FatJet_btagHbb")                                                                                                                                        
"""dtype: Float_t; description: Higgs to BB tagger discriminator"""
FatJet_deepTagMD_H4qvsQCD = NanoAODQuantity("FatJet_deepTagMD_H4qvsQCD")                                                                                                                  
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger H->4q vs QCD discriminator"""
FatJet_deepTagMD_HbbvsQCD = NanoAODQuantity("FatJet_deepTagMD_HbbvsQCD")                                                                                                                  
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger H->bb vs QCD discriminator"""
FatJet_deepTagMD_TvsQCD = NanoAODQuantity("FatJet_deepTagMD_TvsQCD")                                                                                                                      
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger top vs QCD discriminator"""
FatJet_deepTagMD_WvsQCD = NanoAODQuantity("FatJet_deepTagMD_WvsQCD")                                                                                                                      
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger W vs QCD discriminator"""
FatJet_deepTagMD_ZHbbvsQCD = NanoAODQuantity("FatJet_deepTagMD_ZHbbvsQCD")                                                                                                                
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger Z/H->bb vs QCD discriminator"""
FatJet_deepTagMD_ZHccvsQCD = NanoAODQuantity("FatJet_deepTagMD_ZHccvsQCD")                                                                                                                
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger Z/H->cc vs QCD discriminator"""
FatJet_deepTagMD_ZbbvsQCD = NanoAODQuantity("FatJet_deepTagMD_ZbbvsQCD")                                                                                                                  
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger Z->bb vs QCD discriminator"""
FatJet_deepTagMD_ZvsQCD = NanoAODQuantity("FatJet_deepTagMD_ZvsQCD")                                                                                                                      
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger Z vs QCD discriminator"""
FatJet_deepTagMD_bbvsLight = NanoAODQuantity("FatJet_deepTagMD_bbvsLight")                                                                                                                
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->bb vs light flavour discriminator"""
FatJet_deepTagMD_ccvsLight = NanoAODQuantity("FatJet_deepTagMD_ccvsLight")                                                                                                                
"""dtype: Float_t; description: Mass-decorrelated DeepBoostedJet tagger Z/H/gluon->cc vs light flavour discriminator"""
FatJet_deepTag_H = NanoAODQuantity("FatJet_deepTag_H")                                                                                                                                    
"""dtype: Float_t; description: DeepBoostedJet tagger H(bb,cc,4q) sum"""
FatJet_deepTag_QCD = NanoAODQuantity("FatJet_deepTag_QCD")                                                                                                                                
"""dtype: Float_t; description: DeepBoostedJet tagger QCD(bb,cc,b,c,others) sum"""
FatJet_deepTag_QCDothers = NanoAODQuantity("FatJet_deepTag_QCDothers")                                                                                                                    
"""dtype: Float_t; description: DeepBoostedJet tagger QCDothers value"""
FatJet_deepTag_TvsQCD = NanoAODQuantity("FatJet_deepTag_TvsQCD")                                                                                                                          
"""dtype: Float_t; description: DeepBoostedJet tagger top vs QCD discriminator"""
FatJet_deepTag_WvsQCD = NanoAODQuantity("FatJet_deepTag_WvsQCD")                                                                                                                          
"""dtype: Float_t; description: DeepBoostedJet tagger W vs QCD discriminator"""
FatJet_deepTag_ZvsQCD = NanoAODQuantity("FatJet_deepTag_ZvsQCD")                                                                                                                          
"""dtype: Float_t; description: DeepBoostedJet tagger Z vs QCD discriminator"""
FatJet_electronIdx3SJ = NanoAODQuantity("FatJet_electronIdx3SJ")                                                                                                                          
"""dtype: Int_t; description: index of electron matched to jet"""
FatJet_eta = NanoAODQuantity("FatJet_eta")                                                                                                                                                
"""dtype: Float_t; description: eta"""
FatJet_genJetAK8Idx = NanoAODQuantity("FatJet_genJetAK8Idx")                                                                                                                              
"""dtype: Int_t; description: index of matched gen AK8 jet"""
FatJet_hadronFlavour = NanoAODQuantity("FatJet_hadronFlavour")                                                                                                                            
"""dtype: Int_t; description: flavour from hadron ghost clustering"""
FatJet_jetId = NanoAODQuantity("FatJet_jetId")                                                                                                                                            
"""dtype: Int_t; description: Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"""
FatJet_lsf3 = NanoAODQuantity("FatJet_lsf3")                                                                                                                                              
"""dtype: Float_t; description: Lepton Subjet Fraction (3 subjets)"""
FatJet_mass = NanoAODQuantity("FatJet_mass")                                                                                                                                              
"""dtype: Float_t; description: mass"""
FatJet_msoftdrop = NanoAODQuantity("FatJet_msoftdrop")                                                                                                                                    
"""dtype: Float_t; description: Corrected soft drop mass with PUPPI"""
FatJet_muonIdx3SJ = NanoAODQuantity("FatJet_muonIdx3SJ")                                                                                                                                  
"""dtype: Int_t; description: index of muon matched to jet"""
FatJet_n2b1 = NanoAODQuantity("FatJet_n2b1")                                                                                                                                              
"""dtype: Float_t; description: N2 with beta=1"""
FatJet_n3b1 = NanoAODQuantity("FatJet_n3b1")                                                                                                                                              
"""dtype: Float_t; description: N3 with beta=1"""
FatJet_nBHadrons = NanoAODQuantity("FatJet_nBHadrons")                                                                                                                                    
"""dtype: UChar_t; description: number of b-hadrons"""
FatJet_nCHadrons = NanoAODQuantity("FatJet_nCHadrons")                                                                                                                                    
"""dtype: UChar_t; description: number of c-hadrons"""
FatJet_nConstituents = NanoAODQuantity("FatJet_nConstituents")                                                                                                                            
"""dtype: UChar_t; description: Number of particles in the jet"""
FatJet_particleNetMD_QCD = NanoAODQuantity("FatJet_particleNetMD_QCD")                                                                                                                    
"""dtype: Float_t; description: Mass-decorrelated ParticleNet tagger raw QCD score"""
FatJet_particleNetMD_Xbb = NanoAODQuantity("FatJet_particleNetMD_Xbb")                                                                                                                    
"""dtype: Float_t; description: Mass-decorrelated ParticleNet tagger raw X->bb score. For X->bb vs QCD tagging, use Xbb/(Xbb+QCD)"""
FatJet_particleNetMD_Xcc = NanoAODQuantity("FatJet_particleNetMD_Xcc")                                                                                                                    
"""dtype: Float_t; description: Mass-decorrelated ParticleNet tagger raw X->cc score. For X->cc vs QCD tagging, use Xcc/(Xcc+QCD)"""
FatJet_particleNetMD_Xqq = NanoAODQuantity("FatJet_particleNetMD_Xqq")                                                                                                                    
"""dtype: Float_t; description: Mass-decorrelated ParticleNet tagger raw X->qq (uds) score. For X->qq vs QCD tagging, use Xqq/(Xqq+QCD). For W vs QCD tagging, use (Xcc+Xqq)/(Xcc+Xqq+QCD)"""
FatJet_particleNet_H4qvsQCD = NanoAODQuantity("FatJet_particleNet_H4qvsQCD")                                                                                                              
"""dtype: Float_t; description: ParticleNet tagger H(->VV->qqqq) vs QCD discriminator"""
FatJet_particleNet_HbbvsQCD = NanoAODQuantity("FatJet_particleNet_HbbvsQCD")                                                                                                              
"""dtype: Float_t; description: ParticleNet tagger H(->bb) vs QCD discriminator"""
FatJet_particleNet_HccvsQCD = NanoAODQuantity("FatJet_particleNet_HccvsQCD")                                                                                                              
"""dtype: Float_t; description: ParticleNet tagger H(->cc) vs QCD discriminator"""
FatJet_particleNet_QCD = NanoAODQuantity("FatJet_particleNet_QCD")                                                                                                                        
"""dtype: Float_t; description: ParticleNet tagger QCD(bb,cc,b,c,others) sum"""
FatJet_particleNet_TvsQCD = NanoAODQuantity("FatJet_particleNet_TvsQCD")                                                                                                                  
"""dtype: Float_t; description: ParticleNet tagger top vs QCD discriminator"""
FatJet_particleNet_WvsQCD = NanoAODQuantity("FatJet_particleNet_WvsQCD")                                                                                                                  
"""dtype: Float_t; description: ParticleNet tagger W vs QCD discriminator"""
FatJet_particleNet_ZvsQCD = NanoAODQuantity("FatJet_particleNet_ZvsQCD")                                                                                                                  
"""dtype: Float_t; description: ParticleNet tagger Z vs QCD discriminator"""
FatJet_particleNet_mass = NanoAODQuantity("FatJet_particleNet_mass")                                                                                                                      
"""dtype: Float_t; description: ParticleNet mass regression"""
FatJet_phi = NanoAODQuantity("FatJet_phi")                                                                                                                                                
"""dtype: Float_t; description: phi"""
FatJet_pt = NanoAODQuantity("FatJet_pt")                                                                                                                                                  
"""dtype: Float_t; description: pt"""
FatJet_rawFactor = NanoAODQuantity("FatJet_rawFactor")                                                                                                                                    
"""dtype: Float_t; description: 1 - Factor to get back to raw pT"""
FatJet_subJetIdx1 = NanoAODQuantity("FatJet_subJetIdx1")                                                                                                                                  
"""dtype: Int_t; description: index of first subjet"""
FatJet_subJetIdx2 = NanoAODQuantity("FatJet_subJetIdx2")                                                                                                                                  
"""dtype: Int_t; description: index of second subjet"""
FatJet_tau1 = NanoAODQuantity("FatJet_tau1")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (1 axis)"""
FatJet_tau2 = NanoAODQuantity("FatJet_tau2")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (2 axis)"""
FatJet_tau3 = NanoAODQuantity("FatJet_tau3")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (3 axis)"""
FatJet_tau4 = NanoAODQuantity("FatJet_tau4")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (4 axis)"""

Flag_BadChargedCandidateFilter = NanoAODQuantity("Flag_BadChargedCandidateFilter")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_BadChargedCandidateSummer16Filter = NanoAODQuantity("Flag_BadChargedCandidateSummer16Filter")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_BadPFMuonDzFilter = NanoAODQuantity("Flag_BadPFMuonDzFilter")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_BadPFMuonFilter = NanoAODQuantity("Flag_BadPFMuonFilter")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_BadPFMuonSummer16Filter = NanoAODQuantity("Flag_BadPFMuonSummer16Filter")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_CSCTightHalo2015Filter = NanoAODQuantity("Flag_CSCTightHalo2015Filter")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_CSCTightHaloFilter = NanoAODQuantity("Flag_CSCTightHaloFilter")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_CSCTightHaloTrkMuUnvetoFilter = NanoAODQuantity("Flag_CSCTightHaloTrkMuUnvetoFilter")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_EcalDeadCellBoundaryEnergyFilter = NanoAODQuantity("Flag_EcalDeadCellBoundaryEnergyFilter")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_EcalDeadCellTriggerPrimitiveFilter = NanoAODQuantity("Flag_EcalDeadCellTriggerPrimitiveFilter")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_HBHENoiseFilter = NanoAODQuantity("Flag_HBHENoiseFilter")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_HBHENoiseIsoFilter = NanoAODQuantity("Flag_HBHENoiseIsoFilter")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_HcalStripHaloFilter = NanoAODQuantity("Flag_HcalStripHaloFilter")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_METFilters = NanoAODQuantity("Flag_METFilters")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_chargedHadronTrackResolutionFilter = NanoAODQuantity("Flag_chargedHadronTrackResolutionFilter")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_ecalBadCalibFilter = NanoAODQuantity("Flag_ecalBadCalibFilter")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_ecalLaserCorrFilter = NanoAODQuantity("Flag_ecalLaserCorrFilter")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_eeBadScFilter = NanoAODQuantity("Flag_eeBadScFilter")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_globalSuperTightHalo2016Filter = NanoAODQuantity("Flag_globalSuperTightHalo2016Filter")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_globalTightHalo2016Filter = NanoAODQuantity("Flag_globalTightHalo2016Filter")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_goodVertices = NanoAODQuantity("Flag_goodVertices")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_hcalLaserEventFilter = NanoAODQuantity("Flag_hcalLaserEventFilter")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_hfNoisyHitsFilter = NanoAODQuantity("Flag_hfNoisyHitsFilter")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_muonBadTrackFilter = NanoAODQuantity("Flag_muonBadTrackFilter")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_trkPOGFilters = NanoAODQuantity("Flag_trkPOGFilters")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_trkPOG_logErrorTooManyClusters = NanoAODQuantity("Flag_trkPOG_logErrorTooManyClusters")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_trkPOG_manystripclus53X = NanoAODQuantity("Flag_trkPOG_manystripclus53X")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""
Flag_trkPOG_toomanystripclus53X = NanoAODQuantity("Flag_trkPOG_toomanystripclus53X")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: PAT)"""

nFsrPhoton = NanoAODQuantity("nFsrPhoton")                                                                                                                                                
"""dtype: UInt_t; description: Final state radiation photons emitted by muons"""
FsrPhoton_dROverEt2 = NanoAODQuantity("FsrPhoton_dROverEt2")                                                                                                                              
"""dtype: Float_t; description: deltaR to associated muon divided by photon et2"""
FsrPhoton_eta = NanoAODQuantity("FsrPhoton_eta")                                                                                                                                          
"""dtype: Float_t; description: eta"""
FsrPhoton_muonIdx = NanoAODQuantity("FsrPhoton_muonIdx")                                                                                                                                  
"""dtype: Int_t; description: index of associated muon"""
FsrPhoton_phi = NanoAODQuantity("FsrPhoton_phi")                                                                                                                                          
"""dtype: Float_t; description: phi"""
FsrPhoton_pt = NanoAODQuantity("FsrPhoton_pt")                                                                                                                                            
"""dtype: Float_t; description: pt"""
FsrPhoton_relIso03 = NanoAODQuantity("FsrPhoton_relIso03")                                                                                                                                
"""dtype: Float_t; description: relative isolation in a 0.3 cone without CHS"""

nGenDressedLepton = NanoAODQuantity("nGenDressedLepton")                                                                                                                                  
"""dtype: UInt_t; description: Dressed leptons from Rivet-based ParticleLevelProducer"""
GenDressedLepton_eta = NanoAODQuantity("GenDressedLepton_eta")                                                                                                                            
"""dtype: Float_t; description: eta"""
GenDressedLepton_hasTauAnc = NanoAODQuantity("GenDressedLepton_hasTauAnc")                                                                                                                
"""dtype: Bool_t; description: true if Dressed lepton has a tau as ancestor"""
GenDressedLepton_mass = NanoAODQuantity("GenDressedLepton_mass")                                                                                                                          
"""dtype: Float_t; description: mass"""
GenDressedLepton_pdgId = NanoAODQuantity("GenDressedLepton_pdgId")                                                                                                                        
"""dtype: Int_t; description: PDG id"""
GenDressedLepton_phi = NanoAODQuantity("GenDressedLepton_phi")                                                                                                                            
"""dtype: Float_t; description: phi"""
GenDressedLepton_pt = NanoAODQuantity("GenDressedLepton_pt")                                                                                                                              
"""dtype: Float_t; description: pt"""

nGenIsolatedPhoton = NanoAODQuantity("nGenIsolatedPhoton")                                                                                                                                
"""dtype: UInt_t; description: Isolated photons from Rivet-based ParticleLevelProducer"""
GenIsolatedPhoton_eta = NanoAODQuantity("GenIsolatedPhoton_eta")                                                                                                                          
"""dtype: Float_t; description: eta"""
GenIsolatedPhoton_mass = NanoAODQuantity("GenIsolatedPhoton_mass")                                                                                                                        
"""dtype: Float_t; description: mass"""
GenIsolatedPhoton_phi = NanoAODQuantity("GenIsolatedPhoton_phi")                                                                                                                          
"""dtype: Float_t; description: phi"""
GenIsolatedPhoton_pt = NanoAODQuantity("GenIsolatedPhoton_pt")                                                                                                                            
"""dtype: Float_t; description: pt"""

nGenJet = NanoAODQuantity("nGenJet")                                                                                                                                                      
"""dtype: UInt_t; description: slimmedGenJets, i.e. ak4 Jets made with visible genparticles"""
GenJet_eta = NanoAODQuantity("GenJet_eta")                                                                                                                                                
"""dtype: Float_t; description: eta"""
GenJet_hadronFlavour = NanoAODQuantity("GenJet_hadronFlavour")                                                                                                                            
"""dtype: UChar_t; description: flavour from hadron ghost clustering"""
GenJet_mass = NanoAODQuantity("GenJet_mass")                                                                                                                                              
"""dtype: Float_t; description: mass"""
GenJet_partonFlavour = NanoAODQuantity("GenJet_partonFlavour")                                                                                                                            
"""dtype: Int_t; description: flavour from parton matching"""
GenJet_phi = NanoAODQuantity("GenJet_phi")                                                                                                                                                
"""dtype: Float_t; description: phi"""
GenJet_pt = NanoAODQuantity("GenJet_pt")                                                                                                                                                  
"""dtype: Float_t; description: pt"""

nGenJetAK8 = NanoAODQuantity("nGenJetAK8")                                                                                                                                                
"""dtype: UInt_t; description: slimmedGenJetsAK8, i.e. ak8 Jets made with visible genparticles"""
GenJetAK8_eta = NanoAODQuantity("GenJetAK8_eta")                                                                                                                                          
"""dtype: Float_t; description: eta"""
GenJetAK8_hadronFlavour = NanoAODQuantity("GenJetAK8_hadronFlavour")                                                                                                                      
"""dtype: UChar_t; description: flavour from hadron ghost clustering"""
GenJetAK8_mass = NanoAODQuantity("GenJetAK8_mass")                                                                                                                                        
"""dtype: Float_t; description: mass"""
GenJetAK8_partonFlavour = NanoAODQuantity("GenJetAK8_partonFlavour")                                                                                                                      
"""dtype: Int_t; description: flavour from parton matching"""
GenJetAK8_phi = NanoAODQuantity("GenJetAK8_phi")                                                                                                                                          
"""dtype: Float_t; description: phi"""
GenJetAK8_pt = NanoAODQuantity("GenJetAK8_pt")                                                                                                                                            
"""dtype: Float_t; description: pt"""

GenMET_phi = NanoAODQuantity("GenMET_phi")                                                                                                                                                
"""dtype: Float_t; description: phi"""
GenMET_pt = NanoAODQuantity("GenMET_pt")                                                                                                                                                  
"""dtype: Float_t; description: pt"""

nGenPart = NanoAODQuantity("nGenPart")                                                                                                                                                    
"""dtype: UInt_t; description: interesting gen particles """
GenPart_eta = NanoAODQuantity("GenPart_eta")                                                                                                                                              
"""dtype: Float_t; description: eta"""
GenPart_genPartIdxMother = NanoAODQuantity("GenPart_genPartIdxMother")                                                                                                                    
"""dtype: Int_t; description: index of the mother particle"""
GenPart_mass = NanoAODQuantity("GenPart_mass")                                                                                                                                            
"""dtype: Float_t; description: Mass stored for all particles with the exception of quarks (except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411). For these particles, you can lookup the value from PDG."""
GenPart_pdgId = NanoAODQuantity("GenPart_pdgId")                                                                                                                                          
"""dtype: Int_t; description: PDG id"""
GenPart_phi = NanoAODQuantity("GenPart_phi")                                                                                                                                              
"""dtype: Float_t; description: phi"""
GenPart_pt = NanoAODQuantity("GenPart_pt")                                                                                                                                                
"""dtype: Float_t; description: pt"""
GenPart_status = NanoAODQuantity("GenPart_status")                                                                                                                                        
"""dtype: Int_t; description: Particle status. 1=stable"""
GenPart_statusFlags = NanoAODQuantity("GenPart_statusFlags")                                                                                                                              
"""dtype: Int_t; description: gen status flags stored bitwise, bits are: 0 : isPrompt, 1 : isDecayedLeptonHadron, 2 : isTauDecayProduct, 3 : isPromptTauDecayProduct, 4 : isDirectTauDecayProduct, 5 : isDirectPromptTauDecayProduct, 6 : isDirectHadronDecayProduct, 7 : isHardProcess, 8 : fromHardProcess, 9 : isHardProcessTauDecayProduct, 10 : isDirectHardProcessTauDecayProduct, 11 : fromHardProcessBeforeFSR, 12 : isFirstCopy, 13 : isLastCopy, 14 : isLastCopyBeforeFSR, """

nGenVisTau = NanoAODQuantity("nGenVisTau")                                                                                                                                                
"""dtype: UInt_t; description: gen hadronic taus """
GenVisTau_charge = NanoAODQuantity("GenVisTau_charge")                                                                                                                                    
"""dtype: Int_t; description: charge"""
GenVisTau_eta = NanoAODQuantity("GenVisTau_eta")                                                                                                                                          
"""dtype: Float_t; description: eta"""
GenVisTau_genPartIdxMother = NanoAODQuantity("GenVisTau_genPartIdxMother")                                                                                                                
"""dtype: Int_t; description: index of the mother particle"""
GenVisTau_mass = NanoAODQuantity("GenVisTau_mass")                                                                                                                                        
"""dtype: Float_t; description: mass"""
GenVisTau_phi = NanoAODQuantity("GenVisTau_phi")                                                                                                                                          
"""dtype: Float_t; description: phi"""
GenVisTau_pt = NanoAODQuantity("GenVisTau_pt")                                                                                                                                            
"""dtype: Float_t; description: pt"""
GenVisTau_status = NanoAODQuantity("GenVisTau_status")                                                                                                                                    
"""dtype: Int_t; description: Hadronic tau decay mode. 0=OneProng0PiZero, 1=OneProng1PiZero, 2=OneProng2PiZero, 10=ThreeProng0PiZero, 11=ThreeProng1PiZero, 15=Other"""

GenVtx_t0 = NanoAODQuantity("GenVtx_t0")                                                                                                                                                  
"""dtype: Float_t; description: gen vertex t0"""
GenVtx_x = NanoAODQuantity("GenVtx_x")                                                                                                                                                    
"""dtype: Float_t; description: gen vertex x"""
GenVtx_y = NanoAODQuantity("GenVtx_y")                                                                                                                                                    
"""dtype: Float_t; description: gen vertex y"""
GenVtx_z = NanoAODQuantity("GenVtx_z")                                                                                                                                                    
"""dtype: Float_t; description: gen vertex z"""

Generator_binvar = NanoAODQuantity("Generator_binvar")                                                                                                                                    
"""dtype: Float_t; description: MC generation binning value"""
Generator_id1 = NanoAODQuantity("Generator_id1")                                                                                                                                          
"""dtype: Int_t; description: id of first parton"""
Generator_id2 = NanoAODQuantity("Generator_id2")                                                                                                                                          
"""dtype: Int_t; description: id of second parton"""
Generator_scalePDF = NanoAODQuantity("Generator_scalePDF")                                                                                                                                
"""dtype: Float_t; description: Q2 scale for PDF"""
Generator_weight = NanoAODQuantity("Generator_weight")                                                                                                                                    
"""dtype: Float_t; description: MC generator weight"""
Generator_x1 = NanoAODQuantity("Generator_x1")                                                                                                                                            
"""dtype: Float_t; description: x1 fraction of proton momentum carried by the first parton"""
Generator_x2 = NanoAODQuantity("Generator_x2")                                                                                                                                            
"""dtype: Float_t; description: x2 fraction of proton momentum carried by the second parton"""
Generator_xpdf1 = NanoAODQuantity("Generator_xpdf1")                                                                                                                                      
"""dtype: Float_t; description: x*pdf(x) for the first parton"""
Generator_xpdf2 = NanoAODQuantity("Generator_xpdf2")                                                                                                                                      
"""dtype: Float_t; description: x*pdf(x) for the second parton"""

HLT_AK4CaloJet100 = NanoAODQuantity("HLT_AK4CaloJet100")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4CaloJet120 = NanoAODQuantity("HLT_AK4CaloJet120")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4CaloJet30 = NanoAODQuantity("HLT_AK4CaloJet30")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4CaloJet40 = NanoAODQuantity("HLT_AK4CaloJet40")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4CaloJet50 = NanoAODQuantity("HLT_AK4CaloJet50")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4CaloJet80 = NanoAODQuantity("HLT_AK4CaloJet80")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4PFJet100 = NanoAODQuantity("HLT_AK4PFJet100")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4PFJet120 = NanoAODQuantity("HLT_AK4PFJet120")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4PFJet30 = NanoAODQuantity("HLT_AK4PFJet30")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4PFJet50 = NanoAODQuantity("HLT_AK4PFJet50")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK4PFJet80 = NanoAODQuantity("HLT_AK4PFJet80")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFHT750_TrimMass50 = NanoAODQuantity("HLT_AK8PFHT750_TrimMass50")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFHT800_TrimMass50 = NanoAODQuantity("HLT_AK8PFHT800_TrimMass50")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFHT850_TrimMass50 = NanoAODQuantity("HLT_AK8PFHT850_TrimMass50")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFHT900_TrimMass50 = NanoAODQuantity("HLT_AK8PFHT900_TrimMass50")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet140 = NanoAODQuantity("HLT_AK8PFJet140")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet15 = NanoAODQuantity("HLT_AK8PFJet15")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet200 = NanoAODQuantity("HLT_AK8PFJet200")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet25 = NanoAODQuantity("HLT_AK8PFJet25")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet260 = NanoAODQuantity("HLT_AK8PFJet260")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet320 = NanoAODQuantity("HLT_AK8PFJet320")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1 = NanoAODQuantity("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1")                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17 = NanoAODQuantity("HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17")                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2 = NanoAODQuantity("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 = NanoAODQuantity("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02 = NanoAODQuantity("HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet360_TrimMass30 = NanoAODQuantity("HLT_AK8PFJet360_TrimMass30")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet380_TrimMass30 = NanoAODQuantity("HLT_AK8PFJet380_TrimMass30")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet40 = NanoAODQuantity("HLT_AK8PFJet40")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet400 = NanoAODQuantity("HLT_AK8PFJet400")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet400_TrimMass30 = NanoAODQuantity("HLT_AK8PFJet400_TrimMass30")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet420_TrimMass30 = NanoAODQuantity("HLT_AK8PFJet420_TrimMass30")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet450 = NanoAODQuantity("HLT_AK8PFJet450")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet500 = NanoAODQuantity("HLT_AK8PFJet500")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet550 = NanoAODQuantity("HLT_AK8PFJet550")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet60 = NanoAODQuantity("HLT_AK8PFJet60")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJet80 = NanoAODQuantity("HLT_AK8PFJet80")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd140 = NanoAODQuantity("HLT_AK8PFJetFwd140")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd15 = NanoAODQuantity("HLT_AK8PFJetFwd15")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd200 = NanoAODQuantity("HLT_AK8PFJetFwd200")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd25 = NanoAODQuantity("HLT_AK8PFJetFwd25")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd260 = NanoAODQuantity("HLT_AK8PFJetFwd260")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd320 = NanoAODQuantity("HLT_AK8PFJetFwd320")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd40 = NanoAODQuantity("HLT_AK8PFJetFwd40")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd400 = NanoAODQuantity("HLT_AK8PFJetFwd400")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd450 = NanoAODQuantity("HLT_AK8PFJetFwd450")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd500 = NanoAODQuantity("HLT_AK8PFJetFwd500")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd60 = NanoAODQuantity("HLT_AK8PFJetFwd60")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_AK8PFJetFwd80 = NanoAODQuantity("HLT_AK8PFJetFwd80")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet110_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet110_Mu5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet110_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK4DiJet110_Mu5_noalgo")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet170_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet170_Mu5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet170_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK4DiJet170_Mu5_noalgo")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet20_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet20_Mu5")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet20_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK4DiJet20_Mu5_noalgo")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet40_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet40_Mu5")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet40_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK4DiJet40_Mu5_noalgo")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet70_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4DiJet70_Mu5")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4DiJet70_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK4DiJet70_Mu5_noalgo")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4Jet300_Mu5 = NanoAODQuantity("HLT_BTagMu_AK4Jet300_Mu5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK4Jet300_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK4Jet300_Mu5_noalgo")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK8DiJet170_Mu5 = NanoAODQuantity("HLT_BTagMu_AK8DiJet170_Mu5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK8DiJet170_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK8DiJet170_Mu5_noalgo")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK8Jet170_DoubleMu5 = NanoAODQuantity("HLT_BTagMu_AK8Jet170_DoubleMu5")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK8Jet300_Mu5 = NanoAODQuantity("HLT_BTagMu_AK8Jet300_Mu5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_BTagMu_AK8Jet300_Mu5_noalgo = NanoAODQuantity("HLT_BTagMu_AK8Jet300_Mu5_noalgo")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CDC_L2cosmic_5_er1p0 = NanoAODQuantity("HLT_CDC_L2cosmic_5_er1p0")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CDC_L2cosmic_5p5_er1p0 = NanoAODQuantity("HLT_CDC_L2cosmic_5p5_er1p0")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloJet500_NoJetID = NanoAODQuantity("HLT_CaloJet500_NoJetID")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloJet550_NoJetID = NanoAODQuantity("HLT_CaloJet550_NoJetID")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET100_HBHECleaned = NanoAODQuantity("HLT_CaloMET100_HBHECleaned")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET100_NotCleaned = NanoAODQuantity("HLT_CaloMET100_NotCleaned")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET110_NotCleaned = NanoAODQuantity("HLT_CaloMET110_NotCleaned")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET250_HBHECleaned = NanoAODQuantity("HLT_CaloMET250_HBHECleaned")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET250_NotCleaned = NanoAODQuantity("HLT_CaloMET250_NotCleaned")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET300_HBHECleaned = NanoAODQuantity("HLT_CaloMET300_HBHECleaned")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET350_HBHECleaned = NanoAODQuantity("HLT_CaloMET350_HBHECleaned")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET70_HBHECleaned = NanoAODQuantity("HLT_CaloMET70_HBHECleaned")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET80_HBHECleaned = NanoAODQuantity("HLT_CaloMET80_HBHECleaned")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET80_NotCleaned = NanoAODQuantity("HLT_CaloMET80_NotCleaned")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET90_HBHECleaned = NanoAODQuantity("HLT_CaloMET90_HBHECleaned")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMET90_NotCleaned = NanoAODQuantity("HLT_CaloMET90_NotCleaned")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_CaloMHT90 = NanoAODQuantity("HLT_CaloMHT90")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiEle27_WPTightCaloOnly_L1DoubleEG = NanoAODQuantity("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiJet110_35_Mjj650_PFMET110 = NanoAODQuantity("HLT_DiJet110_35_Mjj650_PFMET110")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiJet110_35_Mjj650_PFMET120 = NanoAODQuantity("HLT_DiJet110_35_Mjj650_PFMET120")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiJet110_35_Mjj650_PFMET130 = NanoAODQuantity("HLT_DiJet110_35_Mjj650_PFMET130")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8 = NanoAODQuantity("HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8")                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiMu9_Ele9_CaloIdL_TrackIdL = NanoAODQuantity("HLT_DiMu9_Ele9_CaloIdL_TrackIdL")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ = NanoAODQuantity("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve100_HFJEC = NanoAODQuantity("HLT_DiPFJetAve100_HFJEC")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve140 = NanoAODQuantity("HLT_DiPFJetAve140")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve160_HFJEC = NanoAODQuantity("HLT_DiPFJetAve160_HFJEC")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve200 = NanoAODQuantity("HLT_DiPFJetAve200")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve220_HFJEC = NanoAODQuantity("HLT_DiPFJetAve220_HFJEC")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve260 = NanoAODQuantity("HLT_DiPFJetAve260")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve300_HFJEC = NanoAODQuantity("HLT_DiPFJetAve300_HFJEC")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve320 = NanoAODQuantity("HLT_DiPFJetAve320")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve40 = NanoAODQuantity("HLT_DiPFJetAve40")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve400 = NanoAODQuantity("HLT_DiPFJetAve400")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve500 = NanoAODQuantity("HLT_DiPFJetAve500")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve60 = NanoAODQuantity("HLT_DiPFJetAve60")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve60_HFJEC = NanoAODQuantity("HLT_DiPFJetAve60_HFJEC")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve80 = NanoAODQuantity("HLT_DiPFJetAve80")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiPFJetAve80_HFJEC = NanoAODQuantity("HLT_DiPFJetAve80_HFJEC")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DiSC30_18_EIso_AND_HE_Mass70 = NanoAODQuantity("HLT_DiSC30_18_EIso_AND_HE_Mass70")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi = NanoAODQuantity("HLT_Dimuon0_Jpsi")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi3p5_Muon2 = NanoAODQuantity("HLT_Dimuon0_Jpsi3p5_Muon2")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi_L1_4R_0er1p5R = NanoAODQuantity("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi_L1_NoOS = NanoAODQuantity("HLT_Dimuon0_Jpsi_L1_NoOS")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi_NoVertexing = NanoAODQuantity("HLT_Dimuon0_Jpsi_NoVertexing")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R = NanoAODQuantity("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R")                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Jpsi_NoVertexing_NoOS = NanoAODQuantity("HLT_Dimuon0_Jpsi_NoVertexing_NoOS")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_LowMass = NanoAODQuantity("HLT_Dimuon0_LowMass")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_LowMass_L1_0er1p5 = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_0er1p5")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_LowMass_L1_0er1p5R = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_0er1p5R")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_LowMass_L1_4 = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_4")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_LowMass_L1_4R = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_4R")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_LowMass_L1_TM530 = NanoAODQuantity("HLT_Dimuon0_LowMass_L1_TM530")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_L1_4p5 = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_L1_4p5NoOS = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5NoOS")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_L1_4p5er2p0 = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5er2p0")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_L1_4p5er2p0M = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_4p5er2p0M")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_L1_5 = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_L1_5M = NanoAODQuantity("HLT_Dimuon0_Upsilon_L1_5M")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_Muon_L1_TM0 = NanoAODQuantity("HLT_Dimuon0_Upsilon_Muon_L1_TM0")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_Muon_NoL1Mass = NanoAODQuantity("HLT_Dimuon0_Upsilon_Muon_NoL1Mass")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon0_Upsilon_NoVertexing = NanoAODQuantity("HLT_Dimuon0_Upsilon_NoVertexing")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon10_PsiPrime_Barrel_Seagulls = NanoAODQuantity("HLT_Dimuon10_PsiPrime_Barrel_Seagulls")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon12_Upsilon_y1p4 = NanoAODQuantity("HLT_Dimuon12_Upsilon_y1p4")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon14_Phi_Barrel_Seagulls = NanoAODQuantity("HLT_Dimuon14_Phi_Barrel_Seagulls")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon18_PsiPrime = NanoAODQuantity("HLT_Dimuon18_PsiPrime")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon18_PsiPrime_noCorrL1 = NanoAODQuantity("HLT_Dimuon18_PsiPrime_noCorrL1")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon20_Jpsi_Barrel_Seagulls = NanoAODQuantity("HLT_Dimuon20_Jpsi_Barrel_Seagulls")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon24_Phi_noCorrL1 = NanoAODQuantity("HLT_Dimuon24_Phi_noCorrL1")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon24_Upsilon_noCorrL1 = NanoAODQuantity("HLT_Dimuon24_Upsilon_noCorrL1")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon25_Jpsi = NanoAODQuantity("HLT_Dimuon25_Jpsi")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Dimuon25_Jpsi_noCorrL1 = NanoAODQuantity("HLT_Dimuon25_Jpsi_noCorrL1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55 = NanoAODQuantity("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55")                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 = NanoAODQuantity("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55")                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto = NanoAODQuantity("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto")                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55 = NanoAODQuantity("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_Mass55")                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 = NanoAODQuantity("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90")                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95 = NanoAODQuantity("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95")                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleEle24_eta2p1_WPTight_Gsf = NanoAODQuantity("HLT_DoubleEle24_eta2p1_WPTight_Gsf")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleEle25_CaloIdL_MW = NanoAODQuantity("HLT_DoubleEle25_CaloIdL_MW")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleEle27_CaloIdL_MW = NanoAODQuantity("HLT_DoubleEle27_CaloIdL_MW")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleEle33_CaloIdL_MW = NanoAODQuantity("HLT_DoubleEle33_CaloIdL_MW")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350 = NanoAODQuantity("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350 = NanoAODQuantity("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350")                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleIsoMu20_eta2p1 = NanoAODQuantity("HLT_DoubleIsoMu20_eta2p1")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu23NoVtx_2Cha = NanoAODQuantity("HLT_DoubleL2Mu23NoVtx_2Cha")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed = NanoAODQuantity("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched = NanoAODQuantity("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched = NanoAODQuantity("HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu25NoVtx_2Cha = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4 = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4 = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched = NanoAODQuantity("HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4 = NanoAODQuantity("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4 = NanoAODQuantity("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleL2Mu50 = NanoAODQuantity("HLT_DoubleL2Mu50")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg = NanoAODQuantity("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg")                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg = NanoAODQuantity("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg")                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg = NanoAODQuantity("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg")                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg = NanoAODQuantity("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg")                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu20_7_Mass0to30_L1_DM4 = NanoAODQuantity("HLT_DoubleMu20_7_Mass0to30_L1_DM4")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu20_7_Mass0to30_L1_DM4EG = NanoAODQuantity("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu20_7_Mass0to30_Photon23 = NanoAODQuantity("HLT_DoubleMu20_7_Mass0to30_Photon23")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi = NanoAODQuantity("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 = NanoAODQuantity("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu33NoFiltersNoVtxDisplaced = NanoAODQuantity("HLT_DoubleMu33NoFiltersNoVtxDisplaced")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_DCA_PFMET50_PFMHT60 = NanoAODQuantity("HLT_DoubleMu3_DCA_PFMET50_PFMHT60")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_DZ_PFMET50_PFMHT60 = NanoAODQuantity("HLT_DoubleMu3_DZ_PFMET50_PFMHT60")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_DZ_PFMET70_PFMHT70 = NanoAODQuantity("HLT_DoubleMu3_DZ_PFMET70_PFMHT70")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_DZ_PFMET90_PFMHT90 = NanoAODQuantity("HLT_DoubleMu3_DZ_PFMET90_PFMHT90")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon = NanoAODQuantity("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_TkMu_DsTau3Mu = NanoAODQuantity("HLT_DoubleMu3_TkMu_DsTau3Mu")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_Trk_Tau3mu = NanoAODQuantity("HLT_DoubleMu3_Trk_Tau3mu")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass = NanoAODQuantity("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu40NoFiltersNoVtxDisplaced = NanoAODQuantity("HLT_DoubleMu40NoFiltersNoVtxDisplaced")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu43NoFiltersNoVtx = NanoAODQuantity("HLT_DoubleMu43NoFiltersNoVtx")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu48NoFiltersNoVtx = NanoAODQuantity("HLT_DoubleMu48NoFiltersNoVtx")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_3_Bs = NanoAODQuantity("HLT_DoubleMu4_3_Bs")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_3_Jpsi = NanoAODQuantity("HLT_DoubleMu4_3_Jpsi")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_JpsiTrkTrk_Displaced = NanoAODQuantity("HLT_DoubleMu4_JpsiTrkTrk_Displaced")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_JpsiTrk_Displaced = NanoAODQuantity("HLT_DoubleMu4_JpsiTrk_Displaced")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_Jpsi_Displaced = NanoAODQuantity("HLT_DoubleMu4_Jpsi_Displaced")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_Jpsi_NoVertexing = NanoAODQuantity("HLT_DoubleMu4_Jpsi_NoVertexing")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_LowMassNonResonantTrk_Displaced = NanoAODQuantity("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced")                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_Mass3p8_DZ_PFHT350 = NanoAODQuantity("HLT_DoubleMu4_Mass3p8_DZ_PFHT350")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu4_PsiPrimeTrk_Displaced = NanoAODQuantity("HLT_DoubleMu4_PsiPrimeTrk_Displaced")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL = NanoAODQuantity("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePFJets100_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_DoublePFJets100_CaloBTagDeepCSV_p71")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71")                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71")                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePFJets200_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_DoublePFJets200_CaloBTagDeepCSV_p71")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePFJets350_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_DoublePFJets350_CaloBTagDeepCSV_p71")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePFJets40_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_DoublePFJets40_CaloBTagDeepCSV_p71")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePhoton33_CaloIdL = NanoAODQuantity("HLT_DoublePhoton33_CaloIdL")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePhoton70 = NanoAODQuantity("HLT_DoublePhoton70")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoublePhoton85 = NanoAODQuantity("HLT_DoublePhoton85")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg = NanoAODQuantity("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg = NanoAODQuantity("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg = NanoAODQuantity("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg = NanoAODQuantity("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_DoubleTrkMu_16_6_NoFiltersNoVtx = NanoAODQuantity("HLT_DoubleTrkMu_16_6_NoFiltersNoVtx")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ECALHT800 = NanoAODQuantity("HLT_ECALHT800")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_EcalCalibration = NanoAODQuantity("HLT_EcalCalibration")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele115_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele115_CaloIdVT_GsfTrkIdT")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30")                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele135_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele135_CaloIdVT_GsfTrkIdT")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele145_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele145_CaloIdVT_GsfTrkIdT")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity("HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30")                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT450")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5")                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_IsoVVVL_PFHT450_PFMET50 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT450_PFMET50")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_IsoVVVL_PFHT600 = NanoAODQuantity("HLT_Ele15_IsoVVVL_PFHT600")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele15_WPLoose_Gsf = NanoAODQuantity("HLT_Ele15_WPLoose_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL = NanoAODQuantity("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele17_CaloIdM_TrackIdM_PFJet30 = NanoAODQuantity("HLT_Ele17_CaloIdM_TrackIdM_PFJet30")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele17_WPLoose_Gsf = NanoAODQuantity("HLT_Ele17_WPLoose_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele200_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele200_CaloIdVT_GsfTrkIdT")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele20_WPLoose_Gsf = NanoAODQuantity("HLT_Ele20_WPLoose_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele20_WPTight_Gsf = NanoAODQuantity("HLT_Ele20_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele20_eta2p1_WPLoose_Gsf = NanoAODQuantity("HLT_Ele20_eta2p1_WPLoose_Gsf")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30")                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele23_CaloIdM_TrackIdM_PFJet30 = NanoAODQuantity("HLT_Ele23_CaloIdM_TrackIdM_PFJet30")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1 = NanoAODQuantity("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1")                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 = NanoAODQuantity("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1")          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1 = NanoAODQuantity("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1")                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 = NanoAODQuantity("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1")        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1 = NanoAODQuantity("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1")                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1 = NanoAODQuantity("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1")          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele250_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele250_CaloIdVT_GsfTrkIdT")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele27_Ele37_CaloIdL_MW = NanoAODQuantity("HLT_Ele27_Ele37_CaloIdL_MW")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele27_WPTight_Gsf = NanoAODQuantity("HLT_Ele27_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele28_HighEta_SC20_Mass55 = NanoAODQuantity("HLT_Ele28_HighEta_SC20_Mass55")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele28_WPTight_Gsf = NanoAODQuantity("HLT_Ele28_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele28_eta2p1_WPTight_Gsf_HT150 = NanoAODQuantity("HLT_Ele28_eta2p1_WPTight_Gsf_HT150")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele300_CaloIdVT_GsfTrkIdT = NanoAODQuantity("HLT_Ele300_CaloIdVT_GsfTrkIdT")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele30_WPTight_Gsf = NanoAODQuantity("HLT_Ele30_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned = NanoAODQuantity("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned")                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele32_WPTight_Gsf = NanoAODQuantity("HLT_Ele32_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele32_WPTight_Gsf_L1DoubleEG = NanoAODQuantity("HLT_Ele32_WPTight_Gsf_L1DoubleEG")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele35_WPTight_Gsf = NanoAODQuantity("HLT_Ele35_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele35_WPTight_Gsf_L1EGMT = NanoAODQuantity("HLT_Ele35_WPTight_Gsf_L1EGMT")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele38_WPTight_Gsf = NanoAODQuantity("HLT_Ele38_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele40_WPTight_Gsf = NanoAODQuantity("HLT_Ele40_WPTight_Gsf")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = NanoAODQuantity("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele50_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Ele50_IsoVVVL_PFHT450")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30 = NanoAODQuantity("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Ele8_CaloIdM_TrackIdM_PFJet30 = NanoAODQuantity("HLT_Ele8_CaloIdM_TrackIdM_PFJet30")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT300_Beamspot = NanoAODQuantity("HLT_HT300_Beamspot")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT400_DisplacedDijet40_DisplacedTrack = NanoAODQuantity("HLT_HT400_DisplacedDijet40_DisplacedTrack")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT425 = NanoAODQuantity("HLT_HT425")                                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT430_DisplacedDijet40_DisplacedTrack = NanoAODQuantity("HLT_HT430_DisplacedDijet40_DisplacedTrack")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT430_DisplacedDijet60_DisplacedTrack = NanoAODQuantity("HLT_HT430_DisplacedDijet60_DisplacedTrack")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT450_Beamspot = NanoAODQuantity("HLT_HT450_Beamspot")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT500_DisplacedDijet40_DisplacedTrack = NanoAODQuantity("HLT_HT500_DisplacedDijet40_DisplacedTrack")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT550_DisplacedDijet60_Inclusive = NanoAODQuantity("HLT_HT550_DisplacedDijet60_Inclusive")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HT650_DisplacedDijet60_Inclusive = NanoAODQuantity("HLT_HT650_DisplacedDijet60_Inclusive")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HcalCalibration = NanoAODQuantity("HLT_HcalCalibration")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HcalIsolatedbunch = NanoAODQuantity("HLT_HcalIsolatedbunch")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HcalNZS = NanoAODQuantity("HLT_HcalNZS")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_HcalPhiSym = NanoAODQuantity("HLT_HcalPhiSym")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20 = NanoAODQuantity("HLT_IsoMu20")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 = NanoAODQuantity("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 = NanoAODQuantity("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1")                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1 = NanoAODQuantity("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1")                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 = NanoAODQuantity("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1")                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1 = NanoAODQuantity("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1 = NanoAODQuantity("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1")                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24 = NanoAODQuantity("HLT_IsoMu24")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_TwoProngs35 = NanoAODQuantity("HLT_IsoMu24_TwoProngs35")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_eta2p1 = NanoAODQuantity("HLT_IsoMu24_eta2p1")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr = NanoAODQuantity("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1 = NanoAODQuantity("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1")          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 = NanoAODQuantity("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1")                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1 = NanoAODQuantity("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1")            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 = NanoAODQuantity("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1")                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu27 = NanoAODQuantity("HLT_IsoMu27")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1 = NanoAODQuantity("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu27_MET90 = NanoAODQuantity("HLT_IsoMu27_MET90")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1 = NanoAODQuantity("HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1 = NanoAODQuantity("HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoMu30 = NanoAODQuantity("HLT_IsoMu30")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoTrackHB = NanoAODQuantity("HLT_IsoTrackHB")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_IsoTrackHE = NanoAODQuantity("HLT_IsoTrackHE")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1ETMHadSeeds = NanoAODQuantity("HLT_L1ETMHadSeeds")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1NotBptxOR = NanoAODQuantity("HLT_L1NotBptxOR")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1SingleMu18 = NanoAODQuantity("HLT_L1SingleMu18")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1SingleMu25 = NanoAODQuantity("HLT_L1SingleMu25")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1UnpairedBunchBptxMinus = NanoAODQuantity("HLT_L1UnpairedBunchBptxMinus")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1UnpairedBunchBptxPlus = NanoAODQuantity("HLT_L1UnpairedBunchBptxPlus")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 = NanoAODQuantity("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu10 = NanoAODQuantity("HLT_L2Mu10")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu10_NoVertex_NoBPTX = NanoAODQuantity("HLT_L2Mu10_NoVertex_NoBPTX")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu10_NoVertex_NoBPTX3BX = NanoAODQuantity("HLT_L2Mu10_NoVertex_NoBPTX3BX")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu23NoVtx_2Cha = NanoAODQuantity("HLT_L2Mu23NoVtx_2Cha")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu23NoVtx_2Cha_CosmicSeed = NanoAODQuantity("HLT_L2Mu23NoVtx_2Cha_CosmicSeed")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX = NanoAODQuantity("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX = NanoAODQuantity("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_L2Mu50 = NanoAODQuantity("HLT_L2Mu50")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MET105_IsoTrk50 = NanoAODQuantity("HLT_MET105_IsoTrk50")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MET120_IsoTrk50 = NanoAODQuantity("HLT_MET120_IsoTrk50")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1 = NanoAODQuantity("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1")                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr = NanoAODQuantity("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr")                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1 = NanoAODQuantity("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1")                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1 = NanoAODQuantity("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1")                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110 = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120 = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130 = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140 = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90 = NanoAODQuantity("HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight = NanoAODQuantity("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight")                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight = NanoAODQuantity("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight")                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight = NanoAODQuantity("HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight")                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight = NanoAODQuantity("HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight")                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60 = NanoAODQuantity("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60")                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12 = NanoAODQuantity("HLT_Mu12")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 = NanoAODQuantity("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_DoublePhoton20 = NanoAODQuantity("HLT_Mu12_DoublePhoton20")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_IP6_part0 = NanoAODQuantity("HLT_Mu12_IP6_part0")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_IP6_part1 = NanoAODQuantity("HLT_Mu12_IP6_part1")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_IP6_part2 = NanoAODQuantity("HLT_Mu12_IP6_part2")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_IP6_part3 = NanoAODQuantity("HLT_Mu12_IP6_part3")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_IP6_part4 = NanoAODQuantity("HLT_Mu12_IP6_part4")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu15 = NanoAODQuantity("HLT_Mu15")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu15_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT450")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu15_IsoVVVL_PFHT450_PFMET50 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT450_PFMET50")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu15_IsoVVVL_PFHT600 = NanoAODQuantity("HLT_Mu15_IsoVVVL_PFHT600")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17 = NanoAODQuantity("HLT_Mu17")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17_Photon30_IsoCaloId = NanoAODQuantity("HLT_Mu17_Photon30_IsoCaloId")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17_TrkIsoVVL = NanoAODQuantity("HLT_Mu17_TrkIsoVVL")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = NanoAODQuantity("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = NanoAODQuantity("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 = NanoAODQuantity("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 = NanoAODQuantity("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu18_Mu9 = NanoAODQuantity("HLT_Mu18_Mu9")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu18_Mu9_DZ = NanoAODQuantity("HLT_Mu18_Mu9_DZ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu18_Mu9_SameSign = NanoAODQuantity("HLT_Mu18_Mu9_SameSign")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu18_Mu9_SameSign_DZ = NanoAODQuantity("HLT_Mu18_Mu9_SameSign_DZ")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu19 = NanoAODQuantity("HLT_Mu19")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu19_TrkIsoVVL = NanoAODQuantity("HLT_Mu19_TrkIsoVVL")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL = NanoAODQuantity("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ = NanoAODQuantity("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 = NanoAODQuantity("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 = NanoAODQuantity("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu20 = NanoAODQuantity("HLT_Mu20")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu20_Mu10 = NanoAODQuantity("HLT_Mu20_Mu10")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu20_Mu10_DZ = NanoAODQuantity("HLT_Mu20_Mu10_DZ")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu20_Mu10_SameSign = NanoAODQuantity("HLT_Mu20_Mu10_SameSign")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu20_Mu10_SameSign_DZ = NanoAODQuantity("HLT_Mu20_Mu10_SameSign_DZ")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu20_TkMu0_Phi = NanoAODQuantity("HLT_Mu20_TkMu0_Phi")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu23_Mu12 = NanoAODQuantity("HLT_Mu23_Mu12")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu23_Mu12_DZ = NanoAODQuantity("HLT_Mu23_Mu12_DZ")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu23_Mu12_SameSign = NanoAODQuantity("HLT_Mu23_Mu12_SameSign")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu23_Mu12_SameSign_DZ = NanoAODQuantity("HLT_Mu23_Mu12_SameSign_DZ")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu25_TkMu0_Onia = NanoAODQuantity("HLT_Mu25_TkMu0_Onia")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu25_TkMu0_Phi = NanoAODQuantity("HLT_Mu25_TkMu0_Phi")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu27 = NanoAODQuantity("HLT_Mu27")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu27_Ele37_CaloIdL_MW = NanoAODQuantity("HLT_Mu27_Ele37_CaloIdL_MW")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu30_TkMu0_Psi = NanoAODQuantity("HLT_Mu30_TkMu0_Psi")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu30_TkMu0_Upsilon = NanoAODQuantity("HLT_Mu30_TkMu0_Upsilon")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu37_Ele27_CaloIdL_MW = NanoAODQuantity("HLT_Mu37_Ele27_CaloIdL_MW")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu37_TkMu27 = NanoAODQuantity("HLT_Mu37_TkMu27")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL = NanoAODQuantity("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3_L1SingleMu5orSingleMu7 = NanoAODQuantity("HLT_Mu3_L1SingleMu5orSingleMu7")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3_PFJet40 = NanoAODQuantity("HLT_Mu3_PFJet40")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight")                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight")                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight = NanoAODQuantity("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL = NanoAODQuantity("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL = NanoAODQuantity("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL = NanoAODQuantity("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL")                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 = NanoAODQuantity("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60")                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu50 = NanoAODQuantity("HLT_Mu50")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu50_IsoVVVL_PFHT450 = NanoAODQuantity("HLT_Mu50_IsoVVVL_PFHT450")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu55 = NanoAODQuantity("HLT_Mu55")                                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7_IP4_part0 = NanoAODQuantity("HLT_Mu7_IP4_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7_IP4_part1 = NanoAODQuantity("HLT_Mu7_IP4_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7_IP4_part2 = NanoAODQuantity("HLT_Mu7_IP4_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7_IP4_part3 = NanoAODQuantity("HLT_Mu7_IP4_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7_IP4_part4 = NanoAODQuantity("HLT_Mu7_IP4_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_L2Mu2_Jpsi = NanoAODQuantity("HLT_Mu7p5_L2Mu2_Jpsi")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_L2Mu2_Upsilon = NanoAODQuantity("HLT_Mu7p5_L2Mu2_Upsilon")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_Track2_Jpsi = NanoAODQuantity("HLT_Mu7p5_Track2_Jpsi")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_Track2_Upsilon = NanoAODQuantity("HLT_Mu7p5_Track2_Upsilon")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_Track3p5_Jpsi = NanoAODQuantity("HLT_Mu7p5_Track3p5_Jpsi")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_Track3p5_Upsilon = NanoAODQuantity("HLT_Mu7p5_Track3p5_Upsilon")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_Track7_Jpsi = NanoAODQuantity("HLT_Mu7p5_Track7_Jpsi")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu7p5_Track7_Upsilon = NanoAODQuantity("HLT_Mu7p5_Track7_Upsilon")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8 = NanoAODQuantity("HLT_Mu8")                                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_DiEle12_CaloIdL_TrackIdL = NanoAODQuantity("HLT_Mu8_DiEle12_CaloIdL_TrackIdL")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ = NanoAODQuantity("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350 = NanoAODQuantity("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ = NanoAODQuantity("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ")                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP3_part0 = NanoAODQuantity("HLT_Mu8_IP3_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP3_part1 = NanoAODQuantity("HLT_Mu8_IP3_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP3_part2 = NanoAODQuantity("HLT_Mu8_IP3_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP3_part3 = NanoAODQuantity("HLT_Mu8_IP3_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP3_part4 = NanoAODQuantity("HLT_Mu8_IP3_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP5_part0 = NanoAODQuantity("HLT_Mu8_IP5_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP5_part1 = NanoAODQuantity("HLT_Mu8_IP5_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP5_part2 = NanoAODQuantity("HLT_Mu8_IP5_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP5_part3 = NanoAODQuantity("HLT_Mu8_IP5_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP5_part4 = NanoAODQuantity("HLT_Mu8_IP5_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP6_part0 = NanoAODQuantity("HLT_Mu8_IP6_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP6_part1 = NanoAODQuantity("HLT_Mu8_IP6_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP6_part2 = NanoAODQuantity("HLT_Mu8_IP6_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP6_part3 = NanoAODQuantity("HLT_Mu8_IP6_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_IP6_part4 = NanoAODQuantity("HLT_Mu8_IP6_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL = NanoAODQuantity("HLT_Mu8_TrkIsoVVL")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60 = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60")                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30 = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30")                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5 = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5")  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30 = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30")                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5 = NanoAODQuantity("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5")          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP4_part0 = NanoAODQuantity("HLT_Mu9_IP4_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP4_part1 = NanoAODQuantity("HLT_Mu9_IP4_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP4_part2 = NanoAODQuantity("HLT_Mu9_IP4_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP4_part3 = NanoAODQuantity("HLT_Mu9_IP4_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP4_part4 = NanoAODQuantity("HLT_Mu9_IP4_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP5_part0 = NanoAODQuantity("HLT_Mu9_IP5_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP5_part1 = NanoAODQuantity("HLT_Mu9_IP5_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP5_part2 = NanoAODQuantity("HLT_Mu9_IP5_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP5_part3 = NanoAODQuantity("HLT_Mu9_IP5_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP5_part4 = NanoAODQuantity("HLT_Mu9_IP5_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP6_part0 = NanoAODQuantity("HLT_Mu9_IP6_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP6_part1 = NanoAODQuantity("HLT_Mu9_IP6_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP6_part2 = NanoAODQuantity("HLT_Mu9_IP6_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP6_part3 = NanoAODQuantity("HLT_Mu9_IP6_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Mu9_IP6_part4 = NanoAODQuantity("HLT_Mu9_IP6_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_OldMu100 = NanoAODQuantity("HLT_OldMu100")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT1050 = NanoAODQuantity("HLT_PFHT1050")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT180 = NanoAODQuantity("HLT_PFHT180")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT250 = NanoAODQuantity("HLT_PFHT250")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT330PT30_QuadPFJet_75_60_45_40 = NanoAODQuantity("HLT_PFHT330PT30_QuadPFJet_75_60_45_40")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 = NanoAODQuantity("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5")                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT350 = NanoAODQuantity("HLT_PFHT350")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT350MinPFJet15 = NanoAODQuantity("HLT_PFHT350MinPFJet15")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT370 = NanoAODQuantity("HLT_PFHT370")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT400_SixPFJet32 = NanoAODQuantity("HLT_PFHT400_SixPFJet32")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94 = NanoAODQuantity("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94")                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT430 = NanoAODQuantity("HLT_PFHT430")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT450_SixPFJet36 = NanoAODQuantity("HLT_PFHT450_SixPFJet36")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59 = NanoAODQuantity("HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT500_PFMET100_PFMHT100_IDTight = NanoAODQuantity("HLT_PFHT500_PFMET100_PFMHT100_IDTight")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT500_PFMET110_PFMHT110_IDTight = NanoAODQuantity("HLT_PFHT500_PFMET110_PFMHT110_IDTight")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT510 = NanoAODQuantity("HLT_PFHT510")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT590 = NanoAODQuantity("HLT_PFHT590")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT680 = NanoAODQuantity("HLT_PFHT680")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT700_PFMET85_PFMHT85_IDTight = NanoAODQuantity("HLT_PFHT700_PFMET85_PFMHT85_IDTight")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT700_PFMET95_PFMHT95_IDTight = NanoAODQuantity("HLT_PFHT700_PFMET95_PFMHT95_IDTight")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT780 = NanoAODQuantity("HLT_PFHT780")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT800_PFMET75_PFMHT75_IDTight = NanoAODQuantity("HLT_PFHT800_PFMET75_PFMHT75_IDTight")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT800_PFMET85_PFMHT85_IDTight = NanoAODQuantity("HLT_PFHT800_PFMET85_PFMHT85_IDTight")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFHT890 = NanoAODQuantity("HLT_PFHT890")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet140 = NanoAODQuantity("HLT_PFJet140")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet15 = NanoAODQuantity("HLT_PFJet15")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet200 = NanoAODQuantity("HLT_PFJet200")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet25 = NanoAODQuantity("HLT_PFJet25")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet260 = NanoAODQuantity("HLT_PFJet260")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet320 = NanoAODQuantity("HLT_PFJet320")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet40 = NanoAODQuantity("HLT_PFJet40")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet400 = NanoAODQuantity("HLT_PFJet400")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet450 = NanoAODQuantity("HLT_PFJet450")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet500 = NanoAODQuantity("HLT_PFJet500")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet550 = NanoAODQuantity("HLT_PFJet550")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet60 = NanoAODQuantity("HLT_PFJet60")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJet80 = NanoAODQuantity("HLT_PFJet80")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd140 = NanoAODQuantity("HLT_PFJetFwd140")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd15 = NanoAODQuantity("HLT_PFJetFwd15")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd200 = NanoAODQuantity("HLT_PFJetFwd200")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd25 = NanoAODQuantity("HLT_PFJetFwd25")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd260 = NanoAODQuantity("HLT_PFJetFwd260")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd320 = NanoAODQuantity("HLT_PFJetFwd320")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd40 = NanoAODQuantity("HLT_PFJetFwd40")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd400 = NanoAODQuantity("HLT_PFJetFwd400")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd450 = NanoAODQuantity("HLT_PFJetFwd450")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd500 = NanoAODQuantity("HLT_PFJetFwd500")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd60 = NanoAODQuantity("HLT_PFJetFwd60")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFJetFwd80 = NanoAODQuantity("HLT_PFJetFwd80")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1 = NanoAODQuantity("HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET100_PFMHT100_IDTight_PFHT60 = NanoAODQuantity("HLT_PFMET100_PFMHT100_IDTight_PFHT60")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET110_PFMHT110_IDTight = NanoAODQuantity("HLT_PFMET110_PFMHT110_IDTight")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1 = NanoAODQuantity("HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET120_PFMHT120_IDTight = NanoAODQuantity("HLT_PFMET120_PFMHT120_IDTight")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1 = NanoAODQuantity("HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET120_PFMHT120_IDTight_PFHT60 = NanoAODQuantity("HLT_PFMET120_PFMHT120_IDTight_PFHT60")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET130_PFMHT130_IDTight = NanoAODQuantity("HLT_PFMET130_PFMHT130_IDTight")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1 = NanoAODQuantity("HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET140_PFMHT140_IDTight = NanoAODQuantity("HLT_PFMET140_PFMHT140_IDTight")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1 = NanoAODQuantity("HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1")                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET200_HBHECleaned = NanoAODQuantity("HLT_PFMET200_HBHECleaned")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET200_HBHE_BeamHaloCleaned = NanoAODQuantity("HLT_PFMET200_HBHE_BeamHaloCleaned")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET200_NotCleaned = NanoAODQuantity("HLT_PFMET200_NotCleaned")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET250_HBHECleaned = NanoAODQuantity("HLT_PFMET250_HBHECleaned")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMET300_HBHECleaned = NanoAODQuantity("HLT_PFMET300_HBHECleaned")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 = NanoAODQuantity("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETNoMu110_PFMHTNoMu110_IDTight = NanoAODQuantity("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = NanoAODQuantity("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = NanoAODQuantity("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETNoMu130_PFMHTNoMu130_IDTight = NanoAODQuantity("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETNoMu140_PFMHTNoMu140_IDTight = NanoAODQuantity("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60 = NanoAODQuantity("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne110_PFMHT110_IDTight = NanoAODQuantity("HLT_PFMETTypeOne110_PFMHT110_IDTight")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne120_PFMHT120_IDTight = NanoAODQuantity("HLT_PFMETTypeOne120_PFMHT120_IDTight")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60 = NanoAODQuantity("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne130_PFMHT130_IDTight = NanoAODQuantity("HLT_PFMETTypeOne130_PFMHT130_IDTight")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne140_PFMHT140_IDTight = NanoAODQuantity("HLT_PFMETTypeOne140_PFMHT140_IDTight")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned = NanoAODQuantity("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned")                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon100EBHE10 = NanoAODQuantity("HLT_Photon100EBHE10")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon100EB_TightID_TightIso = NanoAODQuantity("HLT_Photon100EB_TightID_TightIso")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon100EEHE10 = NanoAODQuantity("HLT_Photon100EEHE10")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon100EE_TightID_TightIso = NanoAODQuantity("HLT_Photon100EE_TightID_TightIso")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon110EB_TightID_TightIso = NanoAODQuantity("HLT_Photon110EB_TightID_TightIso")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon120 = NanoAODQuantity("HLT_Photon120")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon120EB_TightID_TightIso = NanoAODQuantity("HLT_Photon120EB_TightID_TightIso")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon120_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon120_R9Id90_HE10_IsoM")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon150 = NanoAODQuantity("HLT_Photon150")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon165_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon165_R9Id90_HE10_IsoM")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon175 = NanoAODQuantity("HLT_Photon175")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon20 = NanoAODQuantity("HLT_Photon20")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon200 = NanoAODQuantity("HLT_Photon200")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon20_HoverELoose = NanoAODQuantity("HLT_Photon20_HoverELoose")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon300_NoHE = NanoAODQuantity("HLT_Photon300_NoHE")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon30_HoverELoose = NanoAODQuantity("HLT_Photon30_HoverELoose")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon33 = NanoAODQuantity("HLT_Photon33")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon35_TwoProngs35 = NanoAODQuantity("HLT_Photon35_TwoProngs35")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon50 = NanoAODQuantity("HLT_Photon50")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon50_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon50_R9Id90_HE10_IsoM")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50 = NanoAODQuantity("HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50")                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon60_R9Id90_CaloIdL_IsoL = NanoAODQuantity("HLT_Photon60_R9Id90_CaloIdL_IsoL")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL = NanoAODQuantity("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL")                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15 = NanoAODQuantity("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15")                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon75 = NanoAODQuantity("HLT_Photon75")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon75_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon75_R9Id90_HE10_IsoM")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3 = NanoAODQuantity("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3")                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3 = NanoAODQuantity("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3")                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3 = NanoAODQuantity("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3")                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3 = NanoAODQuantity("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3")                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon90 = NanoAODQuantity("HLT_Photon90")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon90_CaloIdL_PFHT700 = NanoAODQuantity("HLT_Photon90_CaloIdL_PFHT700")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Photon90_R9Id90_HE10_IsoM = NanoAODQuantity("HLT_Photon90_R9Id90_HE10_IsoM")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics = NanoAODQuantity("HLT_Physics")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part0 = NanoAODQuantity("HLT_Physics_part0")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part1 = NanoAODQuantity("HLT_Physics_part1")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part2 = NanoAODQuantity("HLT_Physics_part2")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part3 = NanoAODQuantity("HLT_Physics_part3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part4 = NanoAODQuantity("HLT_Physics_part4")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part5 = NanoAODQuantity("HLT_Physics_part5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part6 = NanoAODQuantity("HLT_Physics_part6")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Physics_part7 = NanoAODQuantity("HLT_Physics_part7")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet103_88_75_15 = NanoAODQuantity("HLT_QuadPFJet103_88_75_15")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = NanoAODQuantity("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2 = NanoAODQuantity("HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet105_88_76_15 = NanoAODQuantity("HLT_QuadPFJet105_88_76_15")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = NanoAODQuantity("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2 = NanoAODQuantity("HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet111_90_80_15 = NanoAODQuantity("HLT_QuadPFJet111_90_80_15")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = NanoAODQuantity("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1")                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2 = NanoAODQuantity("HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet98_83_71_15 = NanoAODQuantity("HLT_QuadPFJet98_83_71_15")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 = NanoAODQuantity("HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1")                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2 = NanoAODQuantity("HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2")                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Random = NanoAODQuantity("HLT_Random")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Rsq0p35 = NanoAODQuantity("HLT_Rsq0p35")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Rsq0p40 = NanoAODQuantity("HLT_Rsq0p40")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_RsqMR300_Rsq0p09_MR200 = NanoAODQuantity("HLT_RsqMR300_Rsq0p09_MR200")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_RsqMR300_Rsq0p09_MR200_4jet = NanoAODQuantity("HLT_RsqMR300_Rsq0p09_MR200_4jet")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_RsqMR320_Rsq0p09_MR200 = NanoAODQuantity("HLT_RsqMR320_Rsq0p09_MR200")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_RsqMR320_Rsq0p09_MR200_4jet = NanoAODQuantity("HLT_RsqMR320_Rsq0p09_MR200_4jet")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_SingleJet30_Mu12_SinglePFJet40 = NanoAODQuantity("HLT_SingleJet30_Mu12_SinglePFJet40")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_SinglePhoton10_Eta3p1ForPPRef = NanoAODQuantity("HLT_SinglePhoton10_Eta3p1ForPPRef")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_SinglePhoton20_Eta3p1ForPPRef = NanoAODQuantity("HLT_SinglePhoton20_Eta3p1ForPPRef")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_SinglePhoton30_Eta3p1ForPPRef = NanoAODQuantity("HLT_SinglePhoton30_Eta3p1ForPPRef")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = NanoAODQuantity("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = NanoAODQuantity("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15 = NanoAODQuantity("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1 = NanoAODQuantity("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TkMu100 = NanoAODQuantity("HLT_TkMu100")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_Trimuon5_3p5_2_Upsilon_Muon = NanoAODQuantity("HLT_Trimuon5_3p5_2_Upsilon_Muon")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon = NanoAODQuantity("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon")                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleJet110_35_35_Mjj650_PFMET110 = NanoAODQuantity("HLT_TripleJet110_35_35_Mjj650_PFMET110")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleJet110_35_35_Mjj650_PFMET120 = NanoAODQuantity("HLT_TripleJet110_35_35_Mjj650_PFMET120")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleJet110_35_35_Mjj650_PFMET130 = NanoAODQuantity("HLT_TripleJet110_35_35_Mjj650_PFMET130")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleMu_10_5_5_DZ = NanoAODQuantity("HLT_TripleMu_10_5_5_DZ")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleMu_12_10_5 = NanoAODQuantity("HLT_TripleMu_12_10_5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleMu_5_3_3_Mass3p8_DCA = NanoAODQuantity("HLT_TripleMu_5_3_3_Mass3p8_DCA")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TripleMu_5_3_3_Mass3p8_DZ = NanoAODQuantity("HLT_TripleMu_5_3_3_Mass3p8_DZ")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TriplePhoton_20_20_20_CaloIdLV2 = NanoAODQuantity("HLT_TriplePhoton_20_20_20_CaloIdLV2")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL = NanoAODQuantity("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL")                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TriplePhoton_30_30_10_CaloIdLV2 = NanoAODQuantity("HLT_TriplePhoton_30_30_10_CaloIdLV2")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL = NanoAODQuantity("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL")                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL = NanoAODQuantity("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx = NanoAODQuantity("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TrkMu16NoFiltersNoVtx = NanoAODQuantity("HLT_TrkMu16NoFiltersNoVtx")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx = NanoAODQuantity("HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx = NanoAODQuantity("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx")                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_TrkMu6NoFiltersNoVtx = NanoAODQuantity("HLT_TrkMu6NoFiltersNoVtx")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_UncorrectedJetE30_NoBPTX = NanoAODQuantity("HLT_UncorrectedJetE30_NoBPTX")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_UncorrectedJetE30_NoBPTX3BX = NanoAODQuantity("HLT_UncorrectedJetE30_NoBPTX3BX")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_UncorrectedJetE60_NoBPTX3BX = NanoAODQuantity("HLT_UncorrectedJetE60_NoBPTX3BX")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_UncorrectedJetE70_NoBPTX3BX = NanoAODQuantity("HLT_UncorrectedJetE70_NoBPTX3BX")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1 = NanoAODQuantity("HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1 = NanoAODQuantity("HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1")                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1 = NanoAODQuantity("HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1")                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias = NanoAODQuantity("HLT_ZeroBias")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_Alignment = NanoAODQuantity("HLT_ZeroBias_Alignment")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_Beamspot = NanoAODQuantity("HLT_ZeroBias_Beamspot")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_FirstBXAfterTrain = NanoAODQuantity("HLT_ZeroBias_FirstBXAfterTrain")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_FirstCollisionAfterAbortGap = NanoAODQuantity("HLT_ZeroBias_FirstCollisionAfterAbortGap")                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_FirstCollisionInTrain = NanoAODQuantity("HLT_ZeroBias_FirstCollisionInTrain")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_IsolatedBunches = NanoAODQuantity("HLT_ZeroBias_IsolatedBunches")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_LastCollisionInTrain = NanoAODQuantity("HLT_ZeroBias_LastCollisionInTrain")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part0 = NanoAODQuantity("HLT_ZeroBias_part0")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part1 = NanoAODQuantity("HLT_ZeroBias_part1")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part2 = NanoAODQuantity("HLT_ZeroBias_part2")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part3 = NanoAODQuantity("HLT_ZeroBias_part3")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part4 = NanoAODQuantity("HLT_ZeroBias_part4")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part5 = NanoAODQuantity("HLT_ZeroBias_part5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part6 = NanoAODQuantity("HLT_ZeroBias_part6")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLT_ZeroBias_part7 = NanoAODQuantity("HLT_ZeroBias_part7")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""

HTXS_Higgs_pt = NanoAODQuantity("HTXS_Higgs_pt")                                                                                                                                          
"""dtype: Float_t; description: pt of the Higgs boson as identified in HTXS"""
HTXS_Higgs_y = NanoAODQuantity("HTXS_Higgs_y")                                                                                                                                            
"""dtype: Float_t; description: rapidity of the Higgs boson as identified in HTXS"""
HTXS_njets25 = NanoAODQuantity("HTXS_njets25")                                                                                                                                            
"""dtype: UChar_t; description: number of jets with pt>25 GeV as identified in HTXS"""
HTXS_njets30 = NanoAODQuantity("HTXS_njets30")                                                                                                                                            
"""dtype: UChar_t; description: number of jets with pt>30 GeV as identified in HTXS"""
HTXS_stage1_1_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_1_cat_pTjet25GeV")                                                                                                            
"""dtype: Int_t; description: HTXS stage-1.1 category(jet pt>25 GeV)"""
HTXS_stage1_1_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_1_cat_pTjet30GeV")                                                                                                            
"""dtype: Int_t; description: HTXS stage-1.1 category(jet pt>30 GeV)"""
HTXS_stage1_1_fine_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_1_fine_cat_pTjet25GeV")                                                                                                  
"""dtype: Int_t; description: HTXS stage-1.1-fine category(jet pt>25 GeV)"""
HTXS_stage1_1_fine_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_1_fine_cat_pTjet30GeV")                                                                                                  
"""dtype: Int_t; description: HTXS stage-1.1-fine category(jet pt>30 GeV)"""
HTXS_stage1_2_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_2_cat_pTjet25GeV")                                                                                                            
"""dtype: Int_t; description: HTXS stage-1.2 category(jet pt>25 GeV)"""
HTXS_stage1_2_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_2_cat_pTjet30GeV")                                                                                                            
"""dtype: Int_t; description: HTXS stage-1.2 category(jet pt>30 GeV)"""
HTXS_stage1_2_fine_cat_pTjet25GeV = NanoAODQuantity("HTXS_stage1_2_fine_cat_pTjet25GeV")                                                                                                  
"""dtype: Int_t; description: HTXS stage-1.2-fine category(jet pt>25 GeV)"""
HTXS_stage1_2_fine_cat_pTjet30GeV = NanoAODQuantity("HTXS_stage1_2_fine_cat_pTjet30GeV")                                                                                                  
"""dtype: Int_t; description: HTXS stage-1.2-fine category(jet pt>30 GeV)"""
HTXS_stage_0 = NanoAODQuantity("HTXS_stage_0")                                                                                                                                            
"""dtype: Int_t; description: HTXS stage-0 category"""
HTXS_stage_1_pTjet25 = NanoAODQuantity("HTXS_stage_1_pTjet25")                                                                                                                            
"""dtype: Int_t; description: HTXS stage-1 category (jet pt>25 GeV)"""
HTXS_stage_1_pTjet30 = NanoAODQuantity("HTXS_stage_1_pTjet30")                                                                                                                            
"""dtype: Int_t; description: HTXS stage-1 category (jet pt>30 GeV)"""

nIsoTrack = NanoAODQuantity("nIsoTrack")                                                                                                                                                  
"""dtype: UInt_t; description: isolated tracks after basic selection (((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && ((abs(dxy) < 0.2 && abs(dz) < 0.1) || pt>15) && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) and lepton veto"""
IsoTrack_charge = NanoAODQuantity("IsoTrack_charge")                                                                                                                                      
"""dtype: Int_t; description: electric charge"""
IsoTrack_dxy = NanoAODQuantity("IsoTrack_dxy")                                                                                                                                            
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm"""
IsoTrack_dz = NanoAODQuantity("IsoTrack_dz")                                                                                                                                              
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm"""
IsoTrack_eta = NanoAODQuantity("IsoTrack_eta")                                                                                                                                            
"""dtype: Float_t; description: eta"""
IsoTrack_fromPV = NanoAODQuantity("IsoTrack_fromPV")                                                                                                                                      
"""dtype: Int_t; description: isolated track comes from PV"""
IsoTrack_isFromLostTrack = NanoAODQuantity("IsoTrack_isFromLostTrack")                                                                                                                    
"""dtype: Bool_t; description: if isolated track comes from a lost track"""
IsoTrack_isHighPurityTrack = NanoAODQuantity("IsoTrack_isHighPurityTrack")                                                                                                                
"""dtype: Bool_t; description: track is high purity"""
IsoTrack_isPFcand = NanoAODQuantity("IsoTrack_isPFcand")                                                                                                                                  
"""dtype: Bool_t; description: if isolated track is a PF candidate"""
IsoTrack_miniPFRelIso_all = NanoAODQuantity("IsoTrack_miniPFRelIso_all")                                                                                                                  
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections)"""
IsoTrack_miniPFRelIso_chg = NanoAODQuantity("IsoTrack_miniPFRelIso_chg")                                                                                                                  
"""dtype: Float_t; description: mini PF relative isolation, charged component"""
IsoTrack_pdgId = NanoAODQuantity("IsoTrack_pdgId")                                                                                                                                        
"""dtype: Int_t; description: PDG id of PF cand"""
IsoTrack_pfRelIso03_all = NanoAODQuantity("IsoTrack_pfRelIso03_all")                                                                                                                      
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (deltaBeta corrections)"""
IsoTrack_pfRelIso03_chg = NanoAODQuantity("IsoTrack_pfRelIso03_chg")                                                                                                                      
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component"""
IsoTrack_phi = NanoAODQuantity("IsoTrack_phi")                                                                                                                                            
"""dtype: Float_t; description: phi"""
IsoTrack_pt = NanoAODQuantity("IsoTrack_pt")                                                                                                                                              
"""dtype: Float_t; description: pt"""

nJet = NanoAODQuantity("nJet")                                                                                                                                                            
"""dtype: UInt_t; description: slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15)"""
Jet_area = NanoAODQuantity("Jet_area")                                                                                                                                                    
"""dtype: Float_t; description: jet catchment area, for JECs"""
Jet_bRegCorr = NanoAODQuantity("Jet_bRegCorr")                                                                                                                                            
"""dtype: Float_t; description: pt correction for b-jet energy regression"""
Jet_bRegRes = NanoAODQuantity("Jet_bRegRes")                                                                                                                                              
"""dtype: Float_t; description: res on pt corrected with b-jet regression"""
Jet_btagCSVV2 = NanoAODQuantity("Jet_btagCSVV2")                                                                                                                                          
"""dtype: Float_t; description:  pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)"""
Jet_btagDeepB = NanoAODQuantity("Jet_btagDeepB")                                                                                                                                          
"""dtype: Float_t; description: DeepCSV b+bb tag discriminator"""
Jet_btagDeepCvB = NanoAODQuantity("Jet_btagDeepCvB")                                                                                                                                      
"""dtype: Float_t; description: DeepCSV c vs b+bb discriminator"""
Jet_btagDeepCvL = NanoAODQuantity("Jet_btagDeepCvL")                                                                                                                                      
"""dtype: Float_t; description: DeepCSV c vs udsg discriminator"""
Jet_btagDeepFlavB = NanoAODQuantity("Jet_btagDeepFlavB")                                                                                                                                  
"""dtype: Float_t; description: DeepJet b+bb+lepb tag discriminator"""
Jet_btagDeepFlavCvB = NanoAODQuantity("Jet_btagDeepFlavCvB")                                                                                                                              
"""dtype: Float_t; description: DeepJet c vs b+bb+lepb discriminator"""
Jet_btagDeepFlavCvL = NanoAODQuantity("Jet_btagDeepFlavCvL")                                                                                                                              
"""dtype: Float_t; description: DeepJet c vs uds+g discriminator"""
Jet_btagDeepFlavQG = NanoAODQuantity("Jet_btagDeepFlavQG")                                                                                                                                
"""dtype: Float_t; description: DeepJet g vs uds discriminator"""
Jet_cRegCorr = NanoAODQuantity("Jet_cRegCorr")                                                                                                                                            
"""dtype: Float_t; description: pt correction for c-jet energy regression"""
Jet_cRegRes = NanoAODQuantity("Jet_cRegRes")                                                                                                                                              
"""dtype: Float_t; description: res on pt corrected with c-jet regression"""
Jet_chEmEF = NanoAODQuantity("Jet_chEmEF")                                                                                                                                                
"""dtype: Float_t; description: charged Electromagnetic Energy Fraction"""
Jet_chFPV0EF = NanoAODQuantity("Jet_chFPV0EF")                                                                                                                                            
"""dtype: Float_t; description: charged fromPV==0 Energy Fraction (energy excluded from CHS jets). Previously called betastar."""
Jet_chHEF = NanoAODQuantity("Jet_chHEF")                                                                                                                                                  
"""dtype: Float_t; description: charged Hadron Energy Fraction"""
Jet_cleanmask = NanoAODQuantity("Jet_cleanmask")                                                                                                                                          
"""dtype: UChar_t; description: simple cleaning mask with priority to leptons"""
Jet_electronIdx1 = NanoAODQuantity("Jet_electronIdx1")                                                                                                                                    
"""dtype: Int_t; description: index of first matching electron"""
Jet_electronIdx2 = NanoAODQuantity("Jet_electronIdx2")                                                                                                                                    
"""dtype: Int_t; description: index of second matching electron"""
Jet_eta = NanoAODQuantity("Jet_eta")                                                                                                                                                      
"""dtype: Float_t; description: eta"""
Jet_genJetIdx = NanoAODQuantity("Jet_genJetIdx")                                                                                                                                          
"""dtype: Int_t; description: index of matched gen jet"""
Jet_hadronFlavour = NanoAODQuantity("Jet_hadronFlavour")                                                                                                                                  
"""dtype: Int_t; description: flavour from hadron ghost clustering"""
Jet_hfadjacentEtaStripsSize = NanoAODQuantity("Jet_hfadjacentEtaStripsSize")                                                                                                              
"""dtype: Int_t; description: eta size of the strips next to the central tower strip in HF (noise discriminating variable) """
Jet_hfcentralEtaStripSize = NanoAODQuantity("Jet_hfcentralEtaStripSize")                                                                                                                  
"""dtype: Int_t; description: eta size of the central tower strip in HF (noise discriminating variable) """
Jet_hfsigmaEtaEta = NanoAODQuantity("Jet_hfsigmaEtaEta")                                                                                                                                  
"""dtype: Float_t; description: sigmaEtaEta for HF jets (noise discriminating variable)"""
Jet_hfsigmaPhiPhi = NanoAODQuantity("Jet_hfsigmaPhiPhi")                                                                                                                                  
"""dtype: Float_t; description: sigmaPhiPhi for HF jets (noise discriminating variable)"""
Jet_jetId = NanoAODQuantity("Jet_jetId")                                                                                                                                                  
"""dtype: Int_t; description: Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"""
Jet_mass = NanoAODQuantity("Jet_mass")                                                                                                                                                    
"""dtype: Float_t; description: mass"""
Jet_muEF = NanoAODQuantity("Jet_muEF")                                                                                                                                                    
"""dtype: Float_t; description: muon Energy Fraction"""
Jet_muonIdx1 = NanoAODQuantity("Jet_muonIdx1")                                                                                                                                            
"""dtype: Int_t; description: index of first matching muon"""
Jet_muonIdx2 = NanoAODQuantity("Jet_muonIdx2")                                                                                                                                            
"""dtype: Int_t; description: index of second matching muon"""
Jet_muonSubtrFactor = NanoAODQuantity("Jet_muonSubtrFactor")                                                                                                                              
"""dtype: Float_t; description: 1-(muon-subtracted raw pt)/(raw pt)"""
Jet_nConstituents = NanoAODQuantity("Jet_nConstituents")                                                                                                                                  
"""dtype: UChar_t; description: Number of particles in the jet"""
Jet_nElectrons = NanoAODQuantity("Jet_nElectrons")                                                                                                                                        
"""dtype: Int_t; description: number of electrons in the jet"""
Jet_nMuons = NanoAODQuantity("Jet_nMuons")                                                                                                                                                
"""dtype: Int_t; description: number of muons in the jet"""
Jet_neEmEF = NanoAODQuantity("Jet_neEmEF")                                                                                                                                                
"""dtype: Float_t; description: neutral Electromagnetic Energy Fraction"""
Jet_neHEF = NanoAODQuantity("Jet_neHEF")                                                                                                                                                  
"""dtype: Float_t; description: neutral Hadron Energy Fraction"""
Jet_partonFlavour = NanoAODQuantity("Jet_partonFlavour")                                                                                                                                  
"""dtype: Int_t; description: flavour from parton matching"""
Jet_phi = NanoAODQuantity("Jet_phi")                                                                                                                                                      
"""dtype: Float_t; description: phi"""
Jet_pt = NanoAODQuantity("Jet_pt")                                                                                                                                                        
"""dtype: Float_t; description: pt"""
Jet_puId = NanoAODQuantity("Jet_puId")                                                                                                                                                    
"""dtype: Int_t; description: Pileup ID flags with 106X (2018) training"""
Jet_puIdDisc = NanoAODQuantity("Jet_puIdDisc")                                                                                                                                            
"""dtype: Float_t; description: Pileup ID discriminant with 106X (2018) training"""
Jet_qgl = NanoAODQuantity("Jet_qgl")                                                                                                                                                      
"""dtype: Float_t; description: Quark vs Gluon likelihood discriminator"""
Jet_rawFactor = NanoAODQuantity("Jet_rawFactor")                                                                                                                                          
"""dtype: Float_t; description: 1 - Factor to get back to raw pT"""

L1_AlwaysTrue = NanoAODQuantity("L1_AlwaysTrue")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_AND_Ref1_VME = NanoAODQuantity("L1_BPTX_AND_Ref1_VME")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_AND_Ref3_VME = NanoAODQuantity("L1_BPTX_AND_Ref3_VME")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_AND_Ref4_VME = NanoAODQuantity("L1_BPTX_AND_Ref4_VME")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_BeamGas_B1_VME = NanoAODQuantity("L1_BPTX_BeamGas_B1_VME")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_BeamGas_B2_VME = NanoAODQuantity("L1_BPTX_BeamGas_B2_VME")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_BeamGas_Ref1_VME = NanoAODQuantity("L1_BPTX_BeamGas_Ref1_VME")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_BeamGas_Ref2_VME = NanoAODQuantity("L1_BPTX_BeamGas_Ref2_VME")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_NotOR_VME = NanoAODQuantity("L1_BPTX_NotOR_VME")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_OR_Ref3_VME = NanoAODQuantity("L1_BPTX_OR_Ref3_VME")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_OR_Ref4_VME = NanoAODQuantity("L1_BPTX_OR_Ref4_VME")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BPTX_RefAND_VME = NanoAODQuantity("L1_BPTX_RefAND_VME")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BptxMinus = NanoAODQuantity("L1_BptxMinus")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BptxOR = NanoAODQuantity("L1_BptxOR")                                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BptxPlus = NanoAODQuantity("L1_BptxPlus")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_BptxXOR = NanoAODQuantity("L1_BptxXOR")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142 = NanoAODQuantity("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142")                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG8er2p5_HTT260er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT260er")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG8er2p5_HTT280er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT280er")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG8er2p5_HTT300er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT300er")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG8er2p5_HTT320er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT320er")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG8er2p5_HTT340er = NanoAODQuantity("L1_DoubleEG8er2p5_HTT340er")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_15_10_er2p5 = NanoAODQuantity("L1_DoubleEG_15_10_er2p5")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_20_10_er2p5 = NanoAODQuantity("L1_DoubleEG_20_10_er2p5")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_22_10_er2p5 = NanoAODQuantity("L1_DoubleEG_22_10_er2p5")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_25_12_er2p5 = NanoAODQuantity("L1_DoubleEG_25_12_er2p5")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_25_14_er2p5 = NanoAODQuantity("L1_DoubleEG_25_14_er2p5")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_27_14_er2p5 = NanoAODQuantity("L1_DoubleEG_27_14_er2p5")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_LooseIso20_10_er2p5 = NanoAODQuantity("L1_DoubleEG_LooseIso20_10_er2p5")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_LooseIso22_10_er2p5 = NanoAODQuantity("L1_DoubleEG_LooseIso22_10_er2p5")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_LooseIso22_12_er2p5 = NanoAODQuantity("L1_DoubleEG_LooseIso22_12_er2p5")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleEG_LooseIso25_12_er2p5 = NanoAODQuantity("L1_DoubleEG_LooseIso25_12_er2p5")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleIsoTau32er2p1 = NanoAODQuantity("L1_DoubleIsoTau32er2p1")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleIsoTau34er2p1 = NanoAODQuantity("L1_DoubleIsoTau34er2p1")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleIsoTau36er2p1 = NanoAODQuantity("L1_DoubleIsoTau36er2p1")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet100er2p3_dEta_Max1p6 = NanoAODQuantity("L1_DoubleJet100er2p3_dEta_Max1p6")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet100er2p5 = NanoAODQuantity("L1_DoubleJet100er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet112er2p3_dEta_Max1p6 = NanoAODQuantity("L1_DoubleJet112er2p3_dEta_Max1p6")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet120er2p5 = NanoAODQuantity("L1_DoubleJet120er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet150er2p5 = NanoAODQuantity("L1_DoubleJet150er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5 = NanoAODQuantity("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5 = NanoAODQuantity("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 = NanoAODQuantity("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 = NanoAODQuantity("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5 = NanoAODQuantity("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5 = NanoAODQuantity("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp = NanoAODQuantity("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp")                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet40er2p5 = NanoAODQuantity("L1_DoubleJet40er2p5")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_100_30_DoubleJet30_Mass_Min620 = NanoAODQuantity("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_110_35_DoubleJet35_Mass_Min620 = NanoAODQuantity("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620 = NanoAODQuantity("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28 = NanoAODQuantity("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28")                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620 = NanoAODQuantity("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28 = NanoAODQuantity("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28")                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ = NanoAODQuantity("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp = NanoAODQuantity("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp")                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_80_30_Mass_Min420_Mu8 = NanoAODQuantity("L1_DoubleJet_80_30_Mass_Min420_Mu8")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleJet_90_30_DoubleJet30_Mass_Min620 = NanoAODQuantity("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620")                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleLooseIsoEG22er2p1 = NanoAODQuantity("L1_DoubleLooseIsoEG22er2p1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleLooseIsoEG24er2p1 = NanoAODQuantity("L1_DoubleLooseIsoEG24er2p1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0 = NanoAODQuantity("L1_DoubleMu0")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0_Mass_Min1 = NanoAODQuantity("L1_DoubleMu0_Mass_Min1")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0_OQ = NanoAODQuantity("L1_DoubleMu0_OQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0_SQ = NanoAODQuantity("L1_DoubleMu0_SQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0_SQ_OS = NanoAODQuantity("L1_DoubleMu0_SQ_OS")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8 = NanoAODQuantity("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er1p5_SQ = NanoAODQuantity("L1_DoubleMu0er1p5_SQ")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er1p5_SQ_OS = NanoAODQuantity("L1_DoubleMu0er1p5_SQ_OS")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er1p5_SQ_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er1p5_SQ_dR_Max1p4")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu0er2p0_SQ_dR_Max1p4 = NanoAODQuantity("L1_DoubleMu0er2p0_SQ_dR_Max1p4")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu10_SQ = NanoAODQuantity("L1_DoubleMu10_SQ")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu18er2p1 = NanoAODQuantity("L1_DoubleMu18er2p1")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_OS_DoubleEG7p5Upsilon = NanoAODQuantity("L1_DoubleMu3_OS_DoubleEG7p5Upsilon")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_ETMHF50_HTT60er = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF50_HTT60er")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5 = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5 = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5")                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5 = NanoAODQuantity("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5")                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_HTT220er = NanoAODQuantity("L1_DoubleMu3_SQ_HTT220er")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_HTT240er = NanoAODQuantity("L1_DoubleMu3_SQ_HTT240er")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_SQ_HTT260er = NanoAODQuantity("L1_DoubleMu3_SQ_HTT260er")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8 = NanoAODQuantity("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4_SQ_EG9er2p5 = NanoAODQuantity("L1_DoubleMu4_SQ_EG9er2p5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4_SQ_OS = NanoAODQuantity("L1_DoubleMu4_SQ_OS")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4_SQ_OS_dR_Max1p2 = NanoAODQuantity("L1_DoubleMu4_SQ_OS_dR_Max1p2")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4p5_SQ_OS = NanoAODQuantity("L1_DoubleMu4p5_SQ_OS")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4p5_SQ_OS_dR_Max1p2 = NanoAODQuantity("L1_DoubleMu4p5_SQ_OS_dR_Max1p2")                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4p5er2p0_SQ_OS = NanoAODQuantity("L1_DoubleMu4p5er2p0_SQ_OS")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18 = NanoAODQuantity("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu5Upsilon_OS_DoubleEG3 = NanoAODQuantity("L1_DoubleMu5Upsilon_OS_DoubleEG3")                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu5_SQ_EG9er2p5 = NanoAODQuantity("L1_DoubleMu5_SQ_EG9er2p5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu9_SQ = NanoAODQuantity("L1_DoubleMu9_SQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu_12_5 = NanoAODQuantity("L1_DoubleMu_12_5")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu_15_5_SQ = NanoAODQuantity("L1_DoubleMu_15_5_SQ")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu_15_7 = NanoAODQuantity("L1_DoubleMu_15_7")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu_15_7_Mass_Min1 = NanoAODQuantity("L1_DoubleMu_15_7_Mass_Min1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleMu_15_7_SQ = NanoAODQuantity("L1_DoubleMu_15_7_SQ")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_DoubleTau70er2p1 = NanoAODQuantity("L1_DoubleTau70er2p1")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETM120 = NanoAODQuantity("L1_ETM120")                                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETM150 = NanoAODQuantity("L1_ETM150")                                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF100 = NanoAODQuantity("L1_ETMHF100")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF100_HTT60er = NanoAODQuantity("L1_ETMHF100_HTT60er")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF110 = NanoAODQuantity("L1_ETMHF110")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF110_HTT60er = NanoAODQuantity("L1_ETMHF110_HTT60er")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF110_HTT60er_NotSecondBunchInTrain = NanoAODQuantity("L1_ETMHF110_HTT60er_NotSecondBunchInTrain")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF120 = NanoAODQuantity("L1_ETMHF120")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF120_HTT60er = NanoAODQuantity("L1_ETMHF120_HTT60er")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF120_NotSecondBunchInTrain = NanoAODQuantity("L1_ETMHF120_NotSecondBunchInTrain")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF130 = NanoAODQuantity("L1_ETMHF130")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF130_HTT60er = NanoAODQuantity("L1_ETMHF130_HTT60er")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF140 = NanoAODQuantity("L1_ETMHF140")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF150 = NanoAODQuantity("L1_ETMHF150")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETMHF90_HTT60er = NanoAODQuantity("L1_ETMHF90_HTT60er")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETT1200 = NanoAODQuantity("L1_ETT1200")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETT1600 = NanoAODQuantity("L1_ETT1600")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ETT2000 = NanoAODQuantity("L1_ETT2000")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_FirstBunchAfterTrain = NanoAODQuantity("L1_FirstBunchAfterTrain")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_FirstBunchBeforeTrain = NanoAODQuantity("L1_FirstBunchBeforeTrain")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_FirstBunchInTrain = NanoAODQuantity("L1_FirstBunchInTrain")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_FirstCollisionInOrbit = NanoAODQuantity("L1_FirstCollisionInOrbit")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_FirstCollisionInTrain = NanoAODQuantity("L1_FirstCollisionInTrain")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HCAL_LaserMon_Trig = NanoAODQuantity("L1_HCAL_LaserMon_Trig")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HCAL_LaserMon_Veto = NanoAODQuantity("L1_HCAL_LaserMon_Veto")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT120er = NanoAODQuantity("L1_HTT120er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT160er = NanoAODQuantity("L1_HTT160er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT200er = NanoAODQuantity("L1_HTT200er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT255er = NanoAODQuantity("L1_HTT255er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT280er = NanoAODQuantity("L1_HTT280er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT280er_QuadJet_70_55_40_35_er2p4 = NanoAODQuantity("L1_HTT280er_QuadJet_70_55_40_35_er2p4")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT320er = NanoAODQuantity("L1_HTT320er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT320er_QuadJet_70_55_40_40_er2p4 = NanoAODQuantity("L1_HTT320er_QuadJet_70_55_40_40_er2p4")                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 = NanoAODQuantity("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = NanoAODQuantity("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT360er = NanoAODQuantity("L1_HTT360er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT400er = NanoAODQuantity("L1_HTT400er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_HTT450er = NanoAODQuantity("L1_HTT450er")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoEG32er2p5_Mt40 = NanoAODQuantity("L1_IsoEG32er2p5_Mt40")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoEG32er2p5_Mt44 = NanoAODQuantity("L1_IsoEG32er2p5_Mt44")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoEG32er2p5_Mt48 = NanoAODQuantity("L1_IsoEG32er2p5_Mt48")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoTau40er2p1_ETMHF100 = NanoAODQuantity("L1_IsoTau40er2p1_ETMHF100")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoTau40er2p1_ETMHF110 = NanoAODQuantity("L1_IsoTau40er2p1_ETMHF110")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoTau40er2p1_ETMHF120 = NanoAODQuantity("L1_IsoTau40er2p1_ETMHF120")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsoTau40er2p1_ETMHF90 = NanoAODQuantity("L1_IsoTau40er2p1_ETMHF90")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_IsolatedBunch = NanoAODQuantity("L1_IsolatedBunch")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LastBunchInTrain = NanoAODQuantity("L1_LastBunchInTrain")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LastCollisionInTrain = NanoAODQuantity("L1_LastCollisionInTrain")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3 = NanoAODQuantity("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3 = NanoAODQuantity("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG24er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG24er2p1_HTT100er")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3 = NanoAODQuantity("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG26er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG26er2p1_HTT100er")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3 = NanoAODQuantity("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG28er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG28er2p1_HTT100er")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3 = NanoAODQuantity("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG30er2p1_HTT100er = NanoAODQuantity("L1_LooseIsoEG30er2p1_HTT100er")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3 = NanoAODQuantity("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3")                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_MinimumBiasHF0_AND_BptxAND = NanoAODQuantity("L1_MinimumBiasHF0_AND_BptxAND")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6 = NanoAODQuantity("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6")                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6 = NanoAODQuantity("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6")                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6 = NanoAODQuantity("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6")                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu18er2p1_Tau24er2p1 = NanoAODQuantity("L1_Mu18er2p1_Tau24er2p1")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu18er2p1_Tau26er2p1 = NanoAODQuantity("L1_Mu18er2p1_Tau26er2p1")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu20_EG10er2p5 = NanoAODQuantity("L1_Mu20_EG10er2p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu22er2p1_IsoTau32er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau32er2p1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu22er2p1_IsoTau34er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau34er2p1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu22er2p1_IsoTau36er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau36er2p1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu22er2p1_IsoTau40er2p1 = NanoAODQuantity("L1_Mu22er2p1_IsoTau40er2p1")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu22er2p1_Tau70er2p1 = NanoAODQuantity("L1_Mu22er2p1_Tau70er2p1")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet120er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet120er2p5_dR_Max0p4")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet120er2p5_dR_Max0p8 = NanoAODQuantity("L1_Mu3_Jet120er2p5_dR_Max0p8")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet16er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet16er2p5_dR_Max0p4")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet30er2p5 = NanoAODQuantity("L1_Mu3_Jet30er2p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet35er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet35er2p5_dR_Max0p4")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet60er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet60er2p5_dR_Max0p4")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3_Jet80er2p5_dR_Max0p4 = NanoAODQuantity("L1_Mu3_Jet80er2p5_dR_Max0p4")                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3er1p5_Jet100er2p5_ETMHF40 = NanoAODQuantity("L1_Mu3er1p5_Jet100er2p5_ETMHF40")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu3er1p5_Jet100er2p5_ETMHF50 = NanoAODQuantity("L1_Mu3er1p5_Jet100er2p5_ETMHF50")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu5_EG23er2p5 = NanoAODQuantity("L1_Mu5_EG23er2p5")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu5_LooseIsoEG20er2p5 = NanoAODQuantity("L1_Mu5_LooseIsoEG20er2p5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu6_DoubleEG10er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG10er2p5")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu6_DoubleEG12er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG12er2p5")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu6_DoubleEG15er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG15er2p5")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu6_DoubleEG17er2p5 = NanoAODQuantity("L1_Mu6_DoubleEG17er2p5")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu6_HTT240er = NanoAODQuantity("L1_Mu6_HTT240er")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu6_HTT250er = NanoAODQuantity("L1_Mu6_HTT250er")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu7_EG23er2p5 = NanoAODQuantity("L1_Mu7_EG23er2p5")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu7_LooseIsoEG20er2p5 = NanoAODQuantity("L1_Mu7_LooseIsoEG20er2p5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_Mu7_LooseIsoEG23er2p5 = NanoAODQuantity("L1_Mu7_LooseIsoEG23er2p5")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_NotBptxOR = NanoAODQuantity("L1_NotBptxOR")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_QuadJet36er2p5_IsoTau52er2p1 = NanoAODQuantity("L1_QuadJet36er2p5_IsoTau52er2p1")                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_QuadJet60er2p5 = NanoAODQuantity("L1_QuadJet60er2p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0 = NanoAODQuantity("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0")                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_QuadMu0 = NanoAODQuantity("L1_QuadMu0")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_QuadMu0_OQ = NanoAODQuantity("L1_QuadMu0_OQ")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_QuadMu0_SQ = NanoAODQuantity("L1_QuadMu0_SQ")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SecondBunchInTrain = NanoAODQuantity("L1_SecondBunchInTrain")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SecondLastBunchInTrain = NanoAODQuantity("L1_SecondLastBunchInTrain")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG10er2p5 = NanoAODQuantity("L1_SingleEG10er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG15er2p5 = NanoAODQuantity("L1_SingleEG15er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG26er2p5 = NanoAODQuantity("L1_SingleEG26er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG34er2p5 = NanoAODQuantity("L1_SingleEG34er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG36er2p5 = NanoAODQuantity("L1_SingleEG36er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG38er2p5 = NanoAODQuantity("L1_SingleEG38er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG40er2p5 = NanoAODQuantity("L1_SingleEG40er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG42er2p5 = NanoAODQuantity("L1_SingleEG42er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG45er2p5 = NanoAODQuantity("L1_SingleEG45er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG50 = NanoAODQuantity("L1_SingleEG50")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG60 = NanoAODQuantity("L1_SingleEG60")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleEG8er2p5 = NanoAODQuantity("L1_SingleEG8er2p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG24er1p5 = NanoAODQuantity("L1_SingleIsoEG24er1p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG24er2p1 = NanoAODQuantity("L1_SingleIsoEG24er2p1")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG26er1p5 = NanoAODQuantity("L1_SingleIsoEG26er1p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG26er2p1 = NanoAODQuantity("L1_SingleIsoEG26er2p1")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG26er2p5 = NanoAODQuantity("L1_SingleIsoEG26er2p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG28er1p5 = NanoAODQuantity("L1_SingleIsoEG28er1p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG28er2p1 = NanoAODQuantity("L1_SingleIsoEG28er2p1")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG28er2p5 = NanoAODQuantity("L1_SingleIsoEG28er2p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG30er2p1 = NanoAODQuantity("L1_SingleIsoEG30er2p1")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG30er2p5 = NanoAODQuantity("L1_SingleIsoEG30er2p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG32er2p1 = NanoAODQuantity("L1_SingleIsoEG32er2p1")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG32er2p5 = NanoAODQuantity("L1_SingleIsoEG32er2p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleIsoEG34er2p5 = NanoAODQuantity("L1_SingleIsoEG34er2p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet10erHE = NanoAODQuantity("L1_SingleJet10erHE")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet120 = NanoAODQuantity("L1_SingleJet120")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet120_FWD3p0 = NanoAODQuantity("L1_SingleJet120_FWD3p0")                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet120er2p5 = NanoAODQuantity("L1_SingleJet120er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet12erHE = NanoAODQuantity("L1_SingleJet12erHE")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet140er2p5 = NanoAODQuantity("L1_SingleJet140er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet140er2p5_ETMHF80 = NanoAODQuantity("L1_SingleJet140er2p5_ETMHF80")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet140er2p5_ETMHF90 = NanoAODQuantity("L1_SingleJet140er2p5_ETMHF90")                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet160er2p5 = NanoAODQuantity("L1_SingleJet160er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet180 = NanoAODQuantity("L1_SingleJet180")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet180er2p5 = NanoAODQuantity("L1_SingleJet180er2p5")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet200 = NanoAODQuantity("L1_SingleJet200")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet20er2p5_NotBptxOR = NanoAODQuantity("L1_SingleJet20er2p5_NotBptxOR")                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet20er2p5_NotBptxOR_3BX = NanoAODQuantity("L1_SingleJet20er2p5_NotBptxOR_3BX")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet35 = NanoAODQuantity("L1_SingleJet35")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet35_FWD3p0 = NanoAODQuantity("L1_SingleJet35_FWD3p0")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet35er2p5 = NanoAODQuantity("L1_SingleJet35er2p5")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet43er2p5_NotBptxOR_3BX = NanoAODQuantity("L1_SingleJet43er2p5_NotBptxOR_3BX")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet46er2p5_NotBptxOR_3BX = NanoAODQuantity("L1_SingleJet46er2p5_NotBptxOR_3BX")                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet60 = NanoAODQuantity("L1_SingleJet60")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet60_FWD3p0 = NanoAODQuantity("L1_SingleJet60_FWD3p0")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet60er2p5 = NanoAODQuantity("L1_SingleJet60er2p5")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet8erHE = NanoAODQuantity("L1_SingleJet8erHE")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet90 = NanoAODQuantity("L1_SingleJet90")                                                                                                                                        
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet90_FWD3p0 = NanoAODQuantity("L1_SingleJet90_FWD3p0")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleJet90er2p5 = NanoAODQuantity("L1_SingleJet90er2p5")                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleLooseIsoEG28er1p5 = NanoAODQuantity("L1_SingleLooseIsoEG28er1p5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleLooseIsoEG30er1p5 = NanoAODQuantity("L1_SingleLooseIsoEG30er1p5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu0_BMTF = NanoAODQuantity("L1_SingleMu0_BMTF")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu0_DQ = NanoAODQuantity("L1_SingleMu0_DQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu0_EMTF = NanoAODQuantity("L1_SingleMu0_EMTF")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu0_OMTF = NanoAODQuantity("L1_SingleMu0_OMTF")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu10er1p5 = NanoAODQuantity("L1_SingleMu10er1p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu12_DQ_BMTF = NanoAODQuantity("L1_SingleMu12_DQ_BMTF")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu12_DQ_EMTF = NanoAODQuantity("L1_SingleMu12_DQ_EMTF")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu12_DQ_OMTF = NanoAODQuantity("L1_SingleMu12_DQ_OMTF")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu12er1p5 = NanoAODQuantity("L1_SingleMu12er1p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu14er1p5 = NanoAODQuantity("L1_SingleMu14er1p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu15_DQ = NanoAODQuantity("L1_SingleMu15_DQ")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu16er1p5 = NanoAODQuantity("L1_SingleMu16er1p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu18 = NanoAODQuantity("L1_SingleMu18")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu18er1p5 = NanoAODQuantity("L1_SingleMu18er1p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu20 = NanoAODQuantity("L1_SingleMu20")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu22 = NanoAODQuantity("L1_SingleMu22")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu22_BMTF = NanoAODQuantity("L1_SingleMu22_BMTF")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu22_EMTF = NanoAODQuantity("L1_SingleMu22_EMTF")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu22_OMTF = NanoAODQuantity("L1_SingleMu22_OMTF")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu25 = NanoAODQuantity("L1_SingleMu25")                                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu3 = NanoAODQuantity("L1_SingleMu3")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu5 = NanoAODQuantity("L1_SingleMu5")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu6er1p5 = NanoAODQuantity("L1_SingleMu6er1p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu7 = NanoAODQuantity("L1_SingleMu7")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu7_DQ = NanoAODQuantity("L1_SingleMu7_DQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu7er1p5 = NanoAODQuantity("L1_SingleMu7er1p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu8er1p5 = NanoAODQuantity("L1_SingleMu8er1p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMu9er1p5 = NanoAODQuantity("L1_SingleMu9er1p5")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuCosmics = NanoAODQuantity("L1_SingleMuCosmics")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuCosmics_BMTF = NanoAODQuantity("L1_SingleMuCosmics_BMTF")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuCosmics_EMTF = NanoAODQuantity("L1_SingleMuCosmics_EMTF")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuCosmics_OMTF = NanoAODQuantity("L1_SingleMuCosmics_OMTF")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuOpen = NanoAODQuantity("L1_SingleMuOpen")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuOpen_NotBptxOR = NanoAODQuantity("L1_SingleMuOpen_NotBptxOR")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuOpen_er1p1_NotBptxOR_3BX = NanoAODQuantity("L1_SingleMuOpen_er1p1_NotBptxOR_3BX")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleMuOpen_er1p4_NotBptxOR_3BX = NanoAODQuantity("L1_SingleMuOpen_er1p4_NotBptxOR_3BX")                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleTau120er2p1 = NanoAODQuantity("L1_SingleTau120er2p1")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_SingleTau130er2p1 = NanoAODQuantity("L1_SingleTau130er2p1")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TOTEM_1 = NanoAODQuantity("L1_TOTEM_1")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TOTEM_2 = NanoAODQuantity("L1_TOTEM_2")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TOTEM_3 = NanoAODQuantity("L1_TOTEM_3")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TOTEM_4 = NanoAODQuantity("L1_TOTEM_4")                                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleEG16er2p5 = NanoAODQuantity("L1_TripleEG16er2p5")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleEG_16_12_8_er2p5 = NanoAODQuantity("L1_TripleEG_16_12_8_er2p5")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleEG_16_15_8_er2p5 = NanoAODQuantity("L1_TripleEG_16_15_8_er2p5")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleEG_18_17_8_er2p5 = NanoAODQuantity("L1_TripleEG_18_17_8_er2p5")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleEG_18_18_12_er2p5 = NanoAODQuantity("L1_TripleEG_18_18_12_er2p5")                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5 = NanoAODQuantity("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5 = NanoAODQuantity("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5")                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 = NanoAODQuantity("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5")                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu0 = NanoAODQuantity("L1_TripleMu0")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu0_OQ = NanoAODQuantity("L1_TripleMu0_OQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu0_SQ = NanoAODQuantity("L1_TripleMu0_SQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu3 = NanoAODQuantity("L1_TripleMu3")                                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu3_SQ = NanoAODQuantity("L1_TripleMu3_SQ")                                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5SQ_3SQ_0OQ = NanoAODQuantity("L1_TripleMu_5SQ_3SQ_0OQ")                                                                                                                      
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9 = NanoAODQuantity("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9")                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9 = NanoAODQuantity("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_3_3 = NanoAODQuantity("L1_TripleMu_5_3_3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_3_3_SQ = NanoAODQuantity("L1_TripleMu_5_3_3_SQ")                                                                                                                            
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_3p5_2p5 = NanoAODQuantity("L1_TripleMu_5_3p5_2p5")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17 = NanoAODQuantity("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17")                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17 = NanoAODQuantity("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17")                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17 = NanoAODQuantity("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17")                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_TripleMu_5_5_3 = NanoAODQuantity("L1_TripleMu_5_5_3")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_UnpairedBunchBptxMinus = NanoAODQuantity("L1_UnpairedBunchBptxMinus")                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_UnpairedBunchBptxPlus = NanoAODQuantity("L1_UnpairedBunchBptxPlus")                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_UnprefireableEvent = NanoAODQuantity("L1_UnprefireableEvent")                                                                                                                          
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ZeroBias = NanoAODQuantity("L1_ZeroBias")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""
L1_ZeroBias_copy = NanoAODQuantity("L1_ZeroBias_copy")                                                                                                                                    
"""dtype: Bool_t; description: Trigger/flag bit (process: NANO)"""

L1PreFiringWeight_Dn = NanoAODQuantity("L1PreFiringWeight_Dn")                                                                                                                            
"""dtype: Float_t; description: L1 pre-firing event correction weight (1-probability), down var."""
L1PreFiringWeight_ECAL_Dn = NanoAODQuantity("L1PreFiringWeight_ECAL_Dn")                                                                                                                  
"""dtype: Float_t; description: ECAL L1 pre-firing event correction weight (1-probability), down var."""
L1PreFiringWeight_ECAL_Nom = NanoAODQuantity("L1PreFiringWeight_ECAL_Nom")                                                                                                                
"""dtype: Float_t; description: ECAL L1 pre-firing event correction weight (1-probability)"""
L1PreFiringWeight_ECAL_Up = NanoAODQuantity("L1PreFiringWeight_ECAL_Up")                                                                                                                  
"""dtype: Float_t; description: ECAL L1 pre-firing event correction weight (1-probability), up var."""
L1PreFiringWeight_Muon_Nom = NanoAODQuantity("L1PreFiringWeight_Muon_Nom")                                                                                                                
"""dtype: Float_t; description: Muon L1 pre-firing event correction weight (1-probability)"""
L1PreFiringWeight_Muon_StatDn = NanoAODQuantity("L1PreFiringWeight_Muon_StatDn")                                                                                                          
"""dtype: Float_t; description: Muon L1 pre-firing event correction weight (1-probability), down var. stat."""
L1PreFiringWeight_Muon_StatUp = NanoAODQuantity("L1PreFiringWeight_Muon_StatUp")                                                                                                          
"""dtype: Float_t; description: Muon L1 pre-firing event correction weight (1-probability), up var. stat."""
L1PreFiringWeight_Muon_SystDn = NanoAODQuantity("L1PreFiringWeight_Muon_SystDn")                                                                                                          
"""dtype: Float_t; description: Muon L1 pre-firing event correction weight (1-probability), down var. syst."""
L1PreFiringWeight_Muon_SystUp = NanoAODQuantity("L1PreFiringWeight_Muon_SystUp")                                                                                                          
"""dtype: Float_t; description: Muon L1 pre-firing event correction weight (1-probability), up var. syst."""
L1PreFiringWeight_Nom = NanoAODQuantity("L1PreFiringWeight_Nom")                                                                                                                          
"""dtype: Float_t; description: L1 pre-firing event correction weight (1-probability)"""
L1PreFiringWeight_Up = NanoAODQuantity("L1PreFiringWeight_Up")                                                                                                                            
"""dtype: Float_t; description: L1 pre-firing event correction weight (1-probability), up var."""

L1Reco_step = NanoAODQuantity("L1Reco_step")                                                                                                                                              
"""dtype: Bool_t; description: Trigger/flag bit (process: RECO)"""

L1simulation_step = NanoAODQuantity("L1simulation_step")                                                                                                                                  
"""dtype: Bool_t; description: Trigger/flag bit (process: DIGI2RAW)"""

LHE_AlphaS = NanoAODQuantity("LHE_AlphaS")                                                                                                                                                
"""dtype: Float_t; description: Per-event alphaS"""
LHE_HT = NanoAODQuantity("LHE_HT")                                                                                                                                                        
"""dtype: Float_t; description: HT, scalar sum of parton pTs at LHE step"""
LHE_HTIncoming = NanoAODQuantity("LHE_HTIncoming")                                                                                                                                        
"""dtype: Float_t; description: HT, scalar sum of parton pTs at LHE step, restricted to partons"""
LHE_Nb = NanoAODQuantity("LHE_Nb")                                                                                                                                                        
"""dtype: UChar_t; description: Number of b partons at LHE step"""
LHE_Nc = NanoAODQuantity("LHE_Nc")                                                                                                                                                        
"""dtype: UChar_t; description: Number of c partons at LHE step"""
LHE_Nglu = NanoAODQuantity("LHE_Nglu")                                                                                                                                                    
"""dtype: UChar_t; description: Number of gluon partons at LHE step"""
LHE_Njets = NanoAODQuantity("LHE_Njets")                                                                                                                                                  
"""dtype: UChar_t; description: Number of jets (partons) at LHE step"""
LHE_NpLO = NanoAODQuantity("LHE_NpLO")                                                                                                                                                    
"""dtype: UChar_t; description: number of partons at LO"""
LHE_NpNLO = NanoAODQuantity("LHE_NpNLO")                                                                                                                                                  
"""dtype: UChar_t; description: number of partons at NLO"""
LHE_Nuds = NanoAODQuantity("LHE_Nuds")                                                                                                                                                    
"""dtype: UChar_t; description: Number of u,d,s partons at LHE step"""
LHE_Vpt = NanoAODQuantity("LHE_Vpt")                                                                                                                                                      
"""dtype: Float_t; description: pT of the W or Z boson at LHE step"""

nLHEPart = NanoAODQuantity("nLHEPart")                                                                                                                                                    
"""dtype: UInt_t; description: """
LHEPart_eta = NanoAODQuantity("LHEPart_eta")                                                                                                                                              
"""dtype: Float_t; description: Pseodorapidity of LHE particles"""
LHEPart_incomingpz = NanoAODQuantity("LHEPart_incomingpz")                                                                                                                                
"""dtype: Float_t; description: Pz of incoming LHE particles"""
LHEPart_mass = NanoAODQuantity("LHEPart_mass")                                                                                                                                            
"""dtype: Float_t; description: Mass of LHE particles"""
LHEPart_pdgId = NanoAODQuantity("LHEPart_pdgId")                                                                                                                                          
"""dtype: Int_t; description: PDG ID of LHE particles"""
LHEPart_phi = NanoAODQuantity("LHEPart_phi")                                                                                                                                              
"""dtype: Float_t; description: Phi of LHE particles"""
LHEPart_pt = NanoAODQuantity("LHEPart_pt")                                                                                                                                                
"""dtype: Float_t; description: Pt of LHE particles"""
LHEPart_spin = NanoAODQuantity("LHEPart_spin")                                                                                                                                            
"""dtype: Int_t; description: Spin of LHE particles"""
LHEPart_status = NanoAODQuantity("LHEPart_status")                                                                                                                                        
"""dtype: Int_t; description: LHE particle status; -1:incoming, 1:outgoing"""

LHEWeight_originalXWGTUP = NanoAODQuantity("LHEWeight_originalXWGTUP")                                                                                                                    
"""dtype: Float_t; description: Nominal event weight in the LHE file"""

nLowPtElectron = NanoAODQuantity("nLowPtElectron")                                                                                                                                        
"""dtype: UInt_t; description: slimmedLowPtElectrons after basic selection (pt > 1. && userFloat('ID') > -0.25)"""
LowPtElectron_ID = NanoAODQuantity("LowPtElectron_ID")                                                                                                                                    
"""dtype: Float_t; description: New ID, BDT (raw) score"""
LowPtElectron_charge = NanoAODQuantity("LowPtElectron_charge")                                                                                                                            
"""dtype: Int_t; description: electric charge"""
LowPtElectron_convVeto = NanoAODQuantity("LowPtElectron_convVeto")                                                                                                                        
"""dtype: Bool_t; description: pass conversion veto"""
LowPtElectron_convVtxRadius = NanoAODQuantity("LowPtElectron_convVtxRadius")                                                                                                              
"""dtype: Float_t; description: conversion vertex radius (cm)"""
LowPtElectron_convWP = NanoAODQuantity("LowPtElectron_convWP")                                                                                                                            
"""dtype: Int_t; description: conversion flag bit map: 1=Veto, 2=Loose, 3=Tight"""
LowPtElectron_deltaEtaSC = NanoAODQuantity("LowPtElectron_deltaEtaSC")                                                                                                                    
"""dtype: Float_t; description: delta eta (SC,ele) with sign"""
LowPtElectron_dxy = NanoAODQuantity("LowPtElectron_dxy")                                                                                                                                  
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm"""
LowPtElectron_dxyErr = NanoAODQuantity("LowPtElectron_dxyErr")                                                                                                                            
"""dtype: Float_t; description: dxy uncertainty, in cm"""
LowPtElectron_dz = NanoAODQuantity("LowPtElectron_dz")                                                                                                                                    
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm"""
LowPtElectron_dzErr = NanoAODQuantity("LowPtElectron_dzErr")                                                                                                                              
"""dtype: Float_t; description: dz uncertainty, in cm"""
LowPtElectron_eInvMinusPInv = NanoAODQuantity("LowPtElectron_eInvMinusPInv")                                                                                                              
"""dtype: Float_t; description: 1/E_SC - 1/p_trk"""
LowPtElectron_embeddedID = NanoAODQuantity("LowPtElectron_embeddedID")                                                                                                                    
"""dtype: Float_t; description: ID, BDT (raw) score"""
LowPtElectron_energyErr = NanoAODQuantity("LowPtElectron_energyErr")                                                                                                                      
"""dtype: Float_t; description: energy error of the cluster-track combination"""
LowPtElectron_eta = NanoAODQuantity("LowPtElectron_eta")                                                                                                                                  
"""dtype: Float_t; description: eta"""
LowPtElectron_genPartFlav = NanoAODQuantity("LowPtElectron_genPartFlav")                                                                                                                  
"""dtype: UChar_t; description: Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 electrons or photons: 1 = prompt electron (including gamma*->mu mu), 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched"""
LowPtElectron_genPartIdx = NanoAODQuantity("LowPtElectron_genPartIdx")                                                                                                                    
"""dtype: Int_t; description: Index into genParticle list for MC matching to status==1 electrons or photons"""
LowPtElectron_hoe = NanoAODQuantity("LowPtElectron_hoe")                                                                                                                                  
"""dtype: Float_t; description: H over E"""
LowPtElectron_lostHits = NanoAODQuantity("LowPtElectron_lostHits")                                                                                                                        
"""dtype: UChar_t; description: number of missing inner hits"""
LowPtElectron_mass = NanoAODQuantity("LowPtElectron_mass")                                                                                                                                
"""dtype: Float_t; description: mass"""
LowPtElectron_miniPFRelIso_all = NanoAODQuantity("LowPtElectron_miniPFRelIso_all")                                                                                                        
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections)"""
LowPtElectron_miniPFRelIso_chg = NanoAODQuantity("LowPtElectron_miniPFRelIso_chg")                                                                                                        
"""dtype: Float_t; description: mini PF relative isolation, charged component"""
LowPtElectron_pdgId = NanoAODQuantity("LowPtElectron_pdgId")                                                                                                                              
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth)"""
LowPtElectron_phi = NanoAODQuantity("LowPtElectron_phi")                                                                                                                                  
"""dtype: Float_t; description: phi"""
LowPtElectron_pt = NanoAODQuantity("LowPtElectron_pt")                                                                                                                                    
"""dtype: Float_t; description: pt"""
LowPtElectron_ptbiased = NanoAODQuantity("LowPtElectron_ptbiased")                                                                                                                        
"""dtype: Float_t; description: ElectronSeed, pT- and dxy- dependent BDT (raw) score"""
LowPtElectron_r9 = NanoAODQuantity("LowPtElectron_r9")                                                                                                                                    
"""dtype: Float_t; description: R9 of the SC, calculated with full 5x5 region"""
LowPtElectron_scEtOverPt = NanoAODQuantity("LowPtElectron_scEtOverPt")                                                                                                                    
"""dtype: Float_t; description: (SC energy)/pt-1"""
LowPtElectron_sieie = NanoAODQuantity("LowPtElectron_sieie")                                                                                                                              
"""dtype: Float_t; description: sigma_IetaIeta of the SC, calculated with full 5x5 region"""
LowPtElectron_unbiased = NanoAODQuantity("LowPtElectron_unbiased")                                                                                                                        
"""dtype: Float_t; description: ElectronSeed, pT- and dxy- agnostic BDT (raw) score"""

MET_MetUnclustEnUpDeltaX = NanoAODQuantity("MET_MetUnclustEnUpDeltaX")                                                                                                                    
"""dtype: Float_t; description: Delta (METx_mod-METx) Unclustered Energy Up"""
MET_MetUnclustEnUpDeltaY = NanoAODQuantity("MET_MetUnclustEnUpDeltaY")                                                                                                                    
"""dtype: Float_t; description: Delta (METy_mod-METy) Unclustered Energy Up"""
MET_covXX = NanoAODQuantity("MET_covXX")                                                                                                                                                  
"""dtype: Float_t; description: xx element of met covariance matrix"""
MET_covXY = NanoAODQuantity("MET_covXY")                                                                                                                                                  
"""dtype: Float_t; description: xy element of met covariance matrix"""
MET_covYY = NanoAODQuantity("MET_covYY")                                                                                                                                                  
"""dtype: Float_t; description: yy element of met covariance matrix"""
MET_fiducialGenPhi = NanoAODQuantity("MET_fiducialGenPhi")                                                                                                                                
"""dtype: Float_t; description: phi"""
MET_fiducialGenPt = NanoAODQuantity("MET_fiducialGenPt")                                                                                                                                  
"""dtype: Float_t; description: pt"""
MET_phi = NanoAODQuantity("MET_phi")                                                                                                                                                      
"""dtype: Float_t; description: phi"""
MET_pt = NanoAODQuantity("MET_pt")                                                                                                                                                        
"""dtype: Float_t; description: pt"""
MET_significance = NanoAODQuantity("MET_significance")                                                                                                                                    
"""dtype: Float_t; description: MET significance"""
MET_sumEt = NanoAODQuantity("MET_sumEt")                                                                                                                                                  
"""dtype: Float_t; description: scalar sum of Et"""
MET_sumPtUnclustered = NanoAODQuantity("MET_sumPtUnclustered")                                                                                                                            
"""dtype: Float_t; description: sumPt used for MET significance"""

nMuon = NanoAODQuantity("nMuon")                                                                                                                                                          
"""dtype: UInt_t; description: slimmedMuons after basic selection (pt > 15 || (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt'))))"""
Muon_charge = NanoAODQuantity("Muon_charge")                                                                                                                                              
"""dtype: Int_t; description: electric charge"""
Muon_cleanmask = NanoAODQuantity("Muon_cleanmask")                                                                                                                                        
"""dtype: UChar_t; description: simple cleaning mask with priority to leptons"""
Muon_dxy = NanoAODQuantity("Muon_dxy")                                                                                                                                                    
"""dtype: Float_t; description: dxy (with sign) wrt first PV, in cm"""
Muon_dxyErr = NanoAODQuantity("Muon_dxyErr")                                                                                                                                              
"""dtype: Float_t; description: dxy uncertainty, in cm"""
Muon_dxybs = NanoAODQuantity("Muon_dxybs")                                                                                                                                                
"""dtype: Float_t; description: dxy (with sign) wrt the beam spot, in cm"""
Muon_dz = NanoAODQuantity("Muon_dz")                                                                                                                                                      
"""dtype: Float_t; description: dz (with sign) wrt first PV, in cm"""
Muon_dzErr = NanoAODQuantity("Muon_dzErr")                                                                                                                                                
"""dtype: Float_t; description: dz uncertainty, in cm"""
Muon_eta = NanoAODQuantity("Muon_eta")                                                                                                                                                    
"""dtype: Float_t; description: eta"""
Muon_fsrPhotonIdx = NanoAODQuantity("Muon_fsrPhotonIdx")                                                                                                                                  
"""dtype: Int_t; description: Index of the associated FSR photon"""
Muon_genPartFlav = NanoAODQuantity("Muon_genPartFlav")                                                                                                                                    
"""dtype: UChar_t; description: Flavour of genParticle for MC matching to status==1 muons: 1 = prompt muon (including gamma*->mu mu), 15 = muon from prompt tau, 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched"""
Muon_genPartIdx = NanoAODQuantity("Muon_genPartIdx")                                                                                                                                      
"""dtype: Int_t; description: Index into genParticle list for MC matching to status==1 muons"""
Muon_highPtId = NanoAODQuantity("Muon_highPtId")                                                                                                                                          
"""dtype: UChar_t; description: high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"""
Muon_highPurity = NanoAODQuantity("Muon_highPurity")                                                                                                                                      
"""dtype: Bool_t; description: inner track is high purity"""
Muon_inTimeMuon = NanoAODQuantity("Muon_inTimeMuon")                                                                                                                                      
"""dtype: Bool_t; description: inTimeMuon ID"""
Muon_ip3d = NanoAODQuantity("Muon_ip3d")                                                                                                                                                  
"""dtype: Float_t; description: 3D impact parameter wrt first PV, in cm"""
Muon_isGlobal = NanoAODQuantity("Muon_isGlobal")                                                                                                                                          
"""dtype: Bool_t; description: muon is global muon"""
Muon_isPFcand = NanoAODQuantity("Muon_isPFcand")                                                                                                                                          
"""dtype: Bool_t; description: muon is PF candidate"""
Muon_isStandalone = NanoAODQuantity("Muon_isStandalone")                                                                                                                                  
"""dtype: Bool_t; description: muon is a standalone muon"""
Muon_isTracker = NanoAODQuantity("Muon_isTracker")                                                                                                                                        
"""dtype: Bool_t; description: muon is tracker muon"""
Muon_jetIdx = NanoAODQuantity("Muon_jetIdx")                                                                                                                                              
"""dtype: Int_t; description: index of the associated jet (-1 if none)"""
Muon_jetNDauCharged = NanoAODQuantity("Muon_jetNDauCharged")                                                                                                                              
"""dtype: UChar_t; description: number of charged daughters of the closest jet"""
Muon_jetPtRelv2 = NanoAODQuantity("Muon_jetPtRelv2")                                                                                                                                      
"""dtype: Float_t; description: Relative momentum of the lepton with respect to the closest jet after subtracting the lepton"""
Muon_jetRelIso = NanoAODQuantity("Muon_jetRelIso")                                                                                                                                        
"""dtype: Float_t; description: Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet)"""
Muon_looseId = NanoAODQuantity("Muon_looseId")                                                                                                                                            
"""dtype: Bool_t; description: muon is loose muon"""
Muon_mass = NanoAODQuantity("Muon_mass")                                                                                                                                                  
"""dtype: Float_t; description: mass"""
Muon_mediumId = NanoAODQuantity("Muon_mediumId")                                                                                                                                          
"""dtype: Bool_t; description: cut-based ID, medium WP"""
Muon_mediumPromptId = NanoAODQuantity("Muon_mediumPromptId")                                                                                                                              
"""dtype: Bool_t; description: cut-based ID, medium prompt WP"""
Muon_miniIsoId = NanoAODQuantity("Muon_miniIsoId")                                                                                                                                        
"""dtype: UChar_t; description: MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"""
Muon_miniPFRelIso_all = NanoAODQuantity("Muon_miniPFRelIso_all")                                                                                                                          
"""dtype: Float_t; description: mini PF relative isolation, total (with scaled rho*EA PU corrections)"""
Muon_miniPFRelIso_chg = NanoAODQuantity("Muon_miniPFRelIso_chg")                                                                                                                          
"""dtype: Float_t; description: mini PF relative isolation, charged component"""
Muon_multiIsoId = NanoAODQuantity("Muon_multiIsoId")                                                                                                                                      
"""dtype: UChar_t; description: MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"""
Muon_mvaId = NanoAODQuantity("Muon_mvaId")                                                                                                                                                
"""dtype: UChar_t; description: Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight, 4=MvaVTight, 5=MvaVVTight)"""
Muon_mvaLowPt = NanoAODQuantity("Muon_mvaLowPt")                                                                                                                                          
"""dtype: Float_t; description: Low pt muon ID score"""
Muon_mvaLowPtId = NanoAODQuantity("Muon_mvaLowPtId")                                                                                                                                      
"""dtype: UChar_t; description: Low Pt Mva ID from miniAOD selector (1=LowPtMvaLoose, 2=LowPtMvaMedium)"""
Muon_mvaTTH = NanoAODQuantity("Muon_mvaTTH")                                                                                                                                              
"""dtype: Float_t; description: TTH MVA lepton ID score"""
Muon_nStations = NanoAODQuantity("Muon_nStations")                                                                                                                                        
"""dtype: Int_t; description: number of matched stations with default arbitration (segment & track)"""
Muon_nTrackerLayers = NanoAODQuantity("Muon_nTrackerLayers")                                                                                                                              
"""dtype: Int_t; description: number of layers in the tracker"""
Muon_pdgId = NanoAODQuantity("Muon_pdgId")                                                                                                                                                
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth)"""
Muon_pfIsoId = NanoAODQuantity("Muon_pfIsoId")                                                                                                                                            
"""dtype: UChar_t; description: PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"""
Muon_pfRelIso03_all = NanoAODQuantity("Muon_pfRelIso03_all")                                                                                                                              
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (deltaBeta corrections)"""
Muon_pfRelIso03_chg = NanoAODQuantity("Muon_pfRelIso03_chg")                                                                                                                              
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component"""
Muon_pfRelIso04_all = NanoAODQuantity("Muon_pfRelIso04_all")                                                                                                                              
"""dtype: Float_t; description: PF relative isolation dR=0.4, total (deltaBeta corrections)"""
Muon_phi = NanoAODQuantity("Muon_phi")                                                                                                                                                    
"""dtype: Float_t; description: phi"""
Muon_pt = NanoAODQuantity("Muon_pt")                                                                                                                                                      
"""dtype: Float_t; description: pt"""
Muon_ptErr = NanoAODQuantity("Muon_ptErr")                                                                                                                                                
"""dtype: Float_t; description: ptError of the muon track"""
Muon_puppiIsoId = NanoAODQuantity("Muon_puppiIsoId")                                                                                                                                      
"""dtype: UChar_t; description: PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight)"""
Muon_segmentComp = NanoAODQuantity("Muon_segmentComp")                                                                                                                                    
"""dtype: Float_t; description: muon segment compatibility"""
Muon_sip3d = NanoAODQuantity("Muon_sip3d")                                                                                                                                                
"""dtype: Float_t; description: 3D impact parameter significance wrt first PV"""
Muon_softId = NanoAODQuantity("Muon_softId")                                                                                                                                              
"""dtype: Bool_t; description: soft cut-based ID"""
Muon_softMva = NanoAODQuantity("Muon_softMva")                                                                                                                                            
"""dtype: Float_t; description: soft MVA ID score"""
Muon_softMvaId = NanoAODQuantity("Muon_softMvaId")                                                                                                                                        
"""dtype: Bool_t; description: soft MVA ID"""
Muon_tightCharge = NanoAODQuantity("Muon_tightCharge")                                                                                                                                    
"""dtype: Int_t; description: Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"""
Muon_tightId = NanoAODQuantity("Muon_tightId")                                                                                                                                            
"""dtype: Bool_t; description: cut-based ID, tight WP"""
Muon_tkIsoId = NanoAODQuantity("Muon_tkIsoId")                                                                                                                                            
"""dtype: UChar_t; description: TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"""
Muon_tkRelIso = NanoAODQuantity("Muon_tkRelIso")                                                                                                                                          
"""dtype: Float_t; description: Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt"""
Muon_triggerIdLoose = NanoAODQuantity("Muon_triggerIdLoose")                                                                                                                              
"""dtype: Bool_t; description: TriggerIdLoose ID"""
Muon_tunepRelPt = NanoAODQuantity("Muon_tunepRelPt")                                                                                                                                      
"""dtype: Float_t; description: TuneP relative pt, tunePpt/pt"""

nOtherPV = NanoAODQuantity("nOtherPV")                                                                                                                                                    
"""dtype: UInt_t; description: """
OtherPV_z = NanoAODQuantity("OtherPV_z")                                                                                                                                                  
"""dtype: Float_t; description: Z position of other primary vertices, excluding the main PV"""

PV_chi2 = NanoAODQuantity("PV_chi2")                                                                                                                                                      
"""dtype: Float_t; description: main primary vertex reduced chi2"""
PV_ndof = NanoAODQuantity("PV_ndof")                                                                                                                                                      
"""dtype: Float_t; description: main primary vertex number of degree of freedom"""
PV_npvs = NanoAODQuantity("PV_npvs")                                                                                                                                                      
"""dtype: Int_t; description: total number of reconstructed primary vertices"""
PV_npvsGood = NanoAODQuantity("PV_npvsGood")                                                                                                                                              
"""dtype: Int_t; description: number of good reconstructed primary vertices. selection:!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"""
PV_score = NanoAODQuantity("PV_score")                                                                                                                                                    
"""dtype: Float_t; description: main primary vertex score, i.e. sum pt2 of clustered objects"""
PV_x = NanoAODQuantity("PV_x")                                                                                                                                                            
"""dtype: Float_t; description: main primary vertex position x coordinate"""
PV_y = NanoAODQuantity("PV_y")                                                                                                                                                            
"""dtype: Float_t; description: main primary vertex position y coordinate"""
PV_z = NanoAODQuantity("PV_z")                                                                                                                                                            
"""dtype: Float_t; description: main primary vertex position z coordinate"""

nPhoton = NanoAODQuantity("nPhoton")                                                                                                                                                      
"""dtype: UInt_t; description: slimmedPhotons after basic selection (pt > 5 )"""
Photon_charge = NanoAODQuantity("Photon_charge")                                                                                                                                          
"""dtype: Int_t; description: electric charge"""
Photon_cleanmask = NanoAODQuantity("Photon_cleanmask")                                                                                                                                    
"""dtype: UChar_t; description: simple cleaning mask with priority to leptons"""
Photon_cutBased = NanoAODQuantity("Photon_cutBased")                                                                                                                                      
"""dtype: Int_t; description: cut-based ID bitmap, Fall17V2, (0:fail, 1:loose, 2:medium, 3:tight)"""
Photon_cutBased_Fall17V1Bitmap = NanoAODQuantity("Photon_cutBased_Fall17V1Bitmap")                                                                                                        
"""dtype: Int_t; description: cut-based ID bitmap, Fall17V1, 2^(0:loose, 1:medium, 2:tight)."""
Photon_dEscaleDown = NanoAODQuantity("Photon_dEscaleDown")                                                                                                                                
"""dtype: Float_t; description: ecal energy scale shifted 1 sigma down (adding gain/stat/syst in quadrature)"""
Photon_dEscaleUp = NanoAODQuantity("Photon_dEscaleUp")                                                                                                                                    
"""dtype: Float_t; description: ecal energy scale shifted 1 sigma up (adding gain/stat/syst in quadrature)"""
Photon_dEsigmaDown = NanoAODQuantity("Photon_dEsigmaDown")                                                                                                                                
"""dtype: Float_t; description: ecal energy smearing value shifted 1 sigma up"""
Photon_dEsigmaUp = NanoAODQuantity("Photon_dEsigmaUp")                                                                                                                                    
"""dtype: Float_t; description: ecal energy smearing value shifted 1 sigma up"""
Photon_eCorr = NanoAODQuantity("Photon_eCorr")                                                                                                                                            
"""dtype: Float_t; description: ratio of the calibrated energy/miniaod energy"""
Photon_electronIdx = NanoAODQuantity("Photon_electronIdx")                                                                                                                                
"""dtype: Int_t; description: index of the associated electron (-1 if none)"""
Photon_electronVeto = NanoAODQuantity("Photon_electronVeto")                                                                                                                              
"""dtype: Bool_t; description: pass electron veto"""
Photon_energyErr = NanoAODQuantity("Photon_energyErr")                                                                                                                                    
"""dtype: Float_t; description: energy error of the cluster from regression"""
Photon_eta = NanoAODQuantity("Photon_eta")                                                                                                                                                
"""dtype: Float_t; description: eta"""
Photon_genPartFlav = NanoAODQuantity("Photon_genPartFlav")                                                                                                                                
"""dtype: UChar_t; description: Flavour of genParticle for MC matching to status==1 photons or electrons: 1 = prompt photon, 11 = prompt electron, 0 = unknown or unmatched"""
Photon_genPartIdx = NanoAODQuantity("Photon_genPartIdx")                                                                                                                                  
"""dtype: Int_t; description: Index into genParticle list for MC matching to status==1 photons or electrons"""
Photon_hoe = NanoAODQuantity("Photon_hoe")                                                                                                                                                
"""dtype: Float_t; description: H over E"""
Photon_isScEtaEB = NanoAODQuantity("Photon_isScEtaEB")                                                                                                                                    
"""dtype: Bool_t; description: is supercluster eta within barrel acceptance"""
Photon_isScEtaEE = NanoAODQuantity("Photon_isScEtaEE")                                                                                                                                    
"""dtype: Bool_t; description: is supercluster eta within endcap acceptance"""
Photon_jetIdx = NanoAODQuantity("Photon_jetIdx")                                                                                                                                          
"""dtype: Int_t; description: index of the associated jet (-1 if none)"""
Photon_mass = NanoAODQuantity("Photon_mass")                                                                                                                                              
"""dtype: Float_t; description: mass"""
Photon_mvaID = NanoAODQuantity("Photon_mvaID")                                                                                                                                            
"""dtype: Float_t; description: MVA ID score, Fall17V2"""
Photon_mvaID_Fall17V1p1 = NanoAODQuantity("Photon_mvaID_Fall17V1p1")                                                                                                                      
"""dtype: Float_t; description: MVA ID score, Fall17V1p1"""
Photon_mvaID_WP80 = NanoAODQuantity("Photon_mvaID_WP80")                                                                                                                                  
"""dtype: Bool_t; description: MVA ID WP80, Fall17V2"""
Photon_mvaID_WP90 = NanoAODQuantity("Photon_mvaID_WP90")                                                                                                                                  
"""dtype: Bool_t; description: MVA ID WP90, Fall17V2"""
Photon_pdgId = NanoAODQuantity("Photon_pdgId")                                                                                                                                            
"""dtype: Int_t; description: PDG code assigned by the event reconstruction (not by MC truth)"""
Photon_pfRelIso03_all = NanoAODQuantity("Photon_pfRelIso03_all")                                                                                                                          
"""dtype: Float_t; description: PF relative isolation dR=0.3, total (with rho*EA PU corrections)"""
Photon_pfRelIso03_chg = NanoAODQuantity("Photon_pfRelIso03_chg")                                                                                                                          
"""dtype: Float_t; description: PF relative isolation dR=0.3, charged component (with rho*EA PU corrections)"""
Photon_phi = NanoAODQuantity("Photon_phi")                                                                                                                                                
"""dtype: Float_t; description: phi"""
Photon_pixelSeed = NanoAODQuantity("Photon_pixelSeed")                                                                                                                                    
"""dtype: Bool_t; description: has pixel seed"""
Photon_pt = NanoAODQuantity("Photon_pt")                                                                                                                                                  
"""dtype: Float_t; description: p_{T}"""
Photon_r9 = NanoAODQuantity("Photon_r9")                                                                                                                                                  
"""dtype: Float_t; description: R9 of the supercluster, calculated with full 5x5 region"""
Photon_seedGain = NanoAODQuantity("Photon_seedGain")                                                                                                                                      
"""dtype: UChar_t; description: Gain of the seed crystal"""
Photon_sieie = NanoAODQuantity("Photon_sieie")                                                                                                                                            
"""dtype: Float_t; description: sigma_IetaIeta of the supercluster, calculated with full 5x5 region"""
Photon_vidNestedWPBitmap = NanoAODQuantity("Photon_vidNestedWPBitmap")                                                                                                                    
"""dtype: Int_t; description: Fall17V2 VID compressed bitmap (MinPtCut,PhoSCEtaMultiRangeCut,PhoSingleTowerHadOverEmCut,PhoFull5x5SigmaIEtaIEtaCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut,PhoGenericRhoPtScaledCut), 2 bits per cut"""

Pileup_gpudensity = NanoAODQuantity("Pileup_gpudensity")                                                                                                                                  
"""dtype: Float_t; description: Generator-level PU vertices / mm"""
Pileup_nPU = NanoAODQuantity("Pileup_nPU")                                                                                                                                                
"""dtype: Int_t; description: the number of pileup interactions that have been added to the event in the current bunch crossing"""
Pileup_nTrueInt = NanoAODQuantity("Pileup_nTrueInt")                                                                                                                                      
"""dtype: Float_t; description: the true mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled"""
Pileup_pudensity = NanoAODQuantity("Pileup_pudensity")                                                                                                                                    
"""dtype: Float_t; description: PU vertices / mm"""
Pileup_sumEOOT = NanoAODQuantity("Pileup_sumEOOT")                                                                                                                                        
"""dtype: Int_t; description: number of early out of time pileup"""
Pileup_sumLOOT = NanoAODQuantity("Pileup_sumLOOT")                                                                                                                                        
"""dtype: Int_t; description: number of late out of time pileup"""

PuppiMET_phi = NanoAODQuantity("PuppiMET_phi")                                                                                                                                            
"""dtype: Float_t; description: phi"""
PuppiMET_phiJERDown = NanoAODQuantity("PuppiMET_phiJERDown")                                                                                                                              
"""dtype: Float_t; description: JER down phi"""
PuppiMET_phiJERUp = NanoAODQuantity("PuppiMET_phiJERUp")                                                                                                                                  
"""dtype: Float_t; description: JER up phi"""
PuppiMET_phiJESDown = NanoAODQuantity("PuppiMET_phiJESDown")                                                                                                                              
"""dtype: Float_t; description: JES down phi"""
PuppiMET_phiJESUp = NanoAODQuantity("PuppiMET_phiJESUp")                                                                                                                                  
"""dtype: Float_t; description: JES up phi"""
PuppiMET_phiUnclusteredDown = NanoAODQuantity("PuppiMET_phiUnclusteredDown")                                                                                                              
"""dtype: Float_t; description: Unclustered down phi"""
PuppiMET_phiUnclusteredUp = NanoAODQuantity("PuppiMET_phiUnclusteredUp")                                                                                                                  
"""dtype: Float_t; description: Unclustered up phi"""
PuppiMET_pt = NanoAODQuantity("PuppiMET_pt")                                                                                                                                              
"""dtype: Float_t; description: pt"""
PuppiMET_ptJERDown = NanoAODQuantity("PuppiMET_ptJERDown")                                                                                                                                
"""dtype: Float_t; description: JER down pt"""
PuppiMET_ptJERUp = NanoAODQuantity("PuppiMET_ptJERUp")                                                                                                                                    
"""dtype: Float_t; description: JER up pt"""
PuppiMET_ptJESDown = NanoAODQuantity("PuppiMET_ptJESDown")                                                                                                                                
"""dtype: Float_t; description: JES down pt"""
PuppiMET_ptJESUp = NanoAODQuantity("PuppiMET_ptJESUp")                                                                                                                                    
"""dtype: Float_t; description: JES up pt"""
PuppiMET_ptUnclusteredDown = NanoAODQuantity("PuppiMET_ptUnclusteredDown")                                                                                                                
"""dtype: Float_t; description: Unclustered down pt"""
PuppiMET_ptUnclusteredUp = NanoAODQuantity("PuppiMET_ptUnclusteredUp")                                                                                                                    
"""dtype: Float_t; description: Unclustered up pt"""
PuppiMET_sumEt = NanoAODQuantity("PuppiMET_sumEt")                                                                                                                                        
"""dtype: Float_t; description: scalar sum of Et"""

RawMET_phi = NanoAODQuantity("RawMET_phi")                                                                                                                                                
"""dtype: Float_t; description: phi"""
RawMET_pt = NanoAODQuantity("RawMET_pt")                                                                                                                                                  
"""dtype: Float_t; description: pt"""
RawMET_sumEt = NanoAODQuantity("RawMET_sumEt")                                                                                                                                            
"""dtype: Float_t; description: scalar sum of Et"""

RawPuppiMET_phi = NanoAODQuantity("RawPuppiMET_phi")                                                                                                                                      
"""dtype: Float_t; description: phi"""
RawPuppiMET_pt = NanoAODQuantity("RawPuppiMET_pt")                                                                                                                                        
"""dtype: Float_t; description: pt"""
RawPuppiMET_sumEt = NanoAODQuantity("RawPuppiMET_sumEt")                                                                                                                                  
"""dtype: Float_t; description: scalar sum of Et"""

nSV = NanoAODQuantity("nSV")                                                                                                                                                              
"""dtype: UInt_t; description: """
SV_charge = NanoAODQuantity("SV_charge")                                                                                                                                                  
"""dtype: Int_t; description: sum of the charge of the SV tracks"""
SV_chi2 = NanoAODQuantity("SV_chi2")                                                                                                                                                      
"""dtype: Float_t; description: reduced chi2, i.e. chi/ndof"""
SV_dlen = NanoAODQuantity("SV_dlen")                                                                                                                                                      
"""dtype: Float_t; description: decay length in cm"""
SV_dlenSig = NanoAODQuantity("SV_dlenSig")                                                                                                                                                
"""dtype: Float_t; description: decay length significance"""
SV_dxy = NanoAODQuantity("SV_dxy")                                                                                                                                                        
"""dtype: Float_t; description: 2D decay length in cm"""
SV_dxySig = NanoAODQuantity("SV_dxySig")                                                                                                                                                  
"""dtype: Float_t; description: 2D decay length significance"""
SV_eta = NanoAODQuantity("SV_eta")                                                                                                                                                        
"""dtype: Float_t; description: eta"""
SV_mass = NanoAODQuantity("SV_mass")                                                                                                                                                      
"""dtype: Float_t; description: mass"""
SV_ndof = NanoAODQuantity("SV_ndof")                                                                                                                                                      
"""dtype: Float_t; description: number of degrees of freedom"""
SV_ntracks = NanoAODQuantity("SV_ntracks")                                                                                                                                                
"""dtype: UChar_t; description: number of tracks"""
SV_pAngle = NanoAODQuantity("SV_pAngle")                                                                                                                                                  
"""dtype: Float_t; description: pointing angle, i.e. acos(p_SV * (SV - PV)) """
SV_phi = NanoAODQuantity("SV_phi")                                                                                                                                                        
"""dtype: Float_t; description: phi"""
SV_pt = NanoAODQuantity("SV_pt")                                                                                                                                                          
"""dtype: Float_t; description: pt"""
SV_x = NanoAODQuantity("SV_x")                                                                                                                                                            
"""dtype: Float_t; description: secondary vertex X position, in cm"""
SV_y = NanoAODQuantity("SV_y")                                                                                                                                                            
"""dtype: Float_t; description: secondary vertex Y position, in cm"""
SV_z = NanoAODQuantity("SV_z")                                                                                                                                                            
"""dtype: Float_t; description: secondary vertex Z position, in cm"""

nSoftActivityJet = NanoAODQuantity("nSoftActivityJet")                                                                                                                                    
"""dtype: UInt_t; description: jets clustered from charged candidates compatible with primary vertex (charge()!=0 && pvAssociationQuality()>=5 && vertexRef().key()==0)"""
SoftActivityJet_eta = NanoAODQuantity("SoftActivityJet_eta")                                                                                                                              
"""dtype: Float_t; description: eta"""
SoftActivityJet_phi = NanoAODQuantity("SoftActivityJet_phi")                                                                                                                              
"""dtype: Float_t; description: phi"""
SoftActivityJet_pt = NanoAODQuantity("SoftActivityJet_pt")                                                                                                                                
"""dtype: Float_t; description: pt"""

nSubGenJetAK8 = NanoAODQuantity("nSubGenJetAK8")                                                                                                                                          
"""dtype: UInt_t; description: slimmedGenJetsAK8SoftDropSubJets, i.e. subjets of ak8 Jets made with visible genparticles"""
SubGenJetAK8_eta = NanoAODQuantity("SubGenJetAK8_eta")                                                                                                                                    
"""dtype: Float_t; description: eta"""
SubGenJetAK8_mass = NanoAODQuantity("SubGenJetAK8_mass")                                                                                                                                  
"""dtype: Float_t; description: mass"""
SubGenJetAK8_phi = NanoAODQuantity("SubGenJetAK8_phi")                                                                                                                                    
"""dtype: Float_t; description: phi"""
SubGenJetAK8_pt = NanoAODQuantity("SubGenJetAK8_pt")                                                                                                                                      
"""dtype: Float_t; description: pt"""

nSubJet = NanoAODQuantity("nSubJet")                                                                                                                                                      
"""dtype: UInt_t; description: slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"""
SubJet_btagCSVV2 = NanoAODQuantity("SubJet_btagCSVV2")                                                                                                                                    
"""dtype: Float_t; description:  pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)"""
SubJet_btagDeepB = NanoAODQuantity("SubJet_btagDeepB")                                                                                                                                    
"""dtype: Float_t; description: DeepCSV b+bb tag discriminator"""
SubJet_eta = NanoAODQuantity("SubJet_eta")                                                                                                                                                
"""dtype: Float_t; description: eta"""
SubJet_hadronFlavour = NanoAODQuantity("SubJet_hadronFlavour")                                                                                                                            
"""dtype: Int_t; description: flavour from hadron ghost clustering"""
SubJet_mass = NanoAODQuantity("SubJet_mass")                                                                                                                                              
"""dtype: Float_t; description: mass"""
SubJet_n2b1 = NanoAODQuantity("SubJet_n2b1")                                                                                                                                              
"""dtype: Float_t; description: N2 with beta=1"""
SubJet_n3b1 = NanoAODQuantity("SubJet_n3b1")                                                                                                                                              
"""dtype: Float_t; description: N3 with beta=1"""
SubJet_nBHadrons = NanoAODQuantity("SubJet_nBHadrons")                                                                                                                                    
"""dtype: UChar_t; description: number of b-hadrons"""
SubJet_nCHadrons = NanoAODQuantity("SubJet_nCHadrons")                                                                                                                                    
"""dtype: UChar_t; description: number of c-hadrons"""
SubJet_phi = NanoAODQuantity("SubJet_phi")                                                                                                                                                
"""dtype: Float_t; description: phi"""
SubJet_pt = NanoAODQuantity("SubJet_pt")                                                                                                                                                  
"""dtype: Float_t; description: pt"""
SubJet_rawFactor = NanoAODQuantity("SubJet_rawFactor")                                                                                                                                    
"""dtype: Float_t; description: 1 - Factor to get back to raw pT"""
SubJet_tau1 = NanoAODQuantity("SubJet_tau1")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (1 axis)"""
SubJet_tau2 = NanoAODQuantity("SubJet_tau2")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (2 axis)"""
SubJet_tau3 = NanoAODQuantity("SubJet_tau3")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (3 axis)"""
SubJet_tau4 = NanoAODQuantity("SubJet_tau4")                                                                                                                                              
"""dtype: Float_t; description: Nsubjettiness (4 axis)"""

nTau = NanoAODQuantity("nTau")                                                                                                                                                            
"""dtype: UInt_t; description: slimmedTaus after basic selection (pt > 18 && tauID('decayModeFindingNewDMs') && (tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') || (tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum'))<2.5) || tauID('byVVVLooseDeepTau2017v2p1VSjet')))"""
Tau_charge = NanoAODQuantity("Tau_charge")                                                                                                                                                
"""dtype: Int_t; description: electric charge"""
Tau_chargedIso = NanoAODQuantity("Tau_chargedIso")                                                                                                                                        
"""dtype: Float_t; description: charged isolation"""
Tau_cleanmask = NanoAODQuantity("Tau_cleanmask")                                                                                                                                          
"""dtype: UChar_t; description: simple cleaning mask with priority to leptons"""
Tau_decayMode = NanoAODQuantity("Tau_decayMode")                                                                                                                                          
"""dtype: Int_t; description: decayMode()"""
Tau_dxy = NanoAODQuantity("Tau_dxy")                                                                                                                                                      
"""dtype: Float_t; description: d_{xy} of lead track with respect to PV, in cm (with sign)"""
Tau_dz = NanoAODQuantity("Tau_dz")                                                                                                                                                        
"""dtype: Float_t; description: d_{z} of lead track with respect to PV, in cm (with sign)"""
Tau_eta = NanoAODQuantity("Tau_eta")                                                                                                                                                      
"""dtype: Float_t; description: eta"""
Tau_genPartFlav = NanoAODQuantity("Tau_genPartFlav")                                                                                                                                      
"""dtype: UChar_t; description: Flavour of genParticle for MC matching to status==2 taus: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched"""
Tau_genPartIdx = NanoAODQuantity("Tau_genPartIdx")                                                                                                                                        
"""dtype: Int_t; description: Index into genParticle list for MC matching to status==2 taus"""
Tau_idAntiEleDeadECal = NanoAODQuantity("Tau_idAntiEleDeadECal")                                                                                                                          
"""dtype: Bool_t; description: Anti-electron dead-ECal discriminator"""
Tau_idAntiMu = NanoAODQuantity("Tau_idAntiMu")                                                                                                                                            
"""dtype: UChar_t; description: Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight"""
Tau_idDecayModeOldDMs = NanoAODQuantity("Tau_idDecayModeOldDMs")                                                                                                                          
"""dtype: Bool_t; description: tauID('decayModeFinding')"""
Tau_idDeepTau2017v2p1VSe = NanoAODQuantity("Tau_idDeepTau2017v2p1VSe")                                                                                                                    
"""dtype: UChar_t; description: byDeepTau2017v2p1VSe ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight"""
Tau_idDeepTau2017v2p1VSjet = NanoAODQuantity("Tau_idDeepTau2017v2p1VSjet")                                                                                                                
"""dtype: UChar_t; description: byDeepTau2017v2p1VSjet ID working points (deepTau2017v2p1): bitmask 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight"""
Tau_idDeepTau2017v2p1VSmu = NanoAODQuantity("Tau_idDeepTau2017v2p1VSmu")                                                                                                                  
"""dtype: UChar_t; description: byDeepTau2017v2p1VSmu ID working points (deepTau2017v2p1): bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight"""
Tau_jetIdx = NanoAODQuantity("Tau_jetIdx")                                                                                                                                                
"""dtype: Int_t; description: index of the associated jet (-1 if none)"""
Tau_leadTkDeltaEta = NanoAODQuantity("Tau_leadTkDeltaEta")                                                                                                                                
"""dtype: Float_t; description: eta of the leading track, minus tau eta"""
Tau_leadTkDeltaPhi = NanoAODQuantity("Tau_leadTkDeltaPhi")                                                                                                                                
"""dtype: Float_t; description: phi of the leading track, minus tau phi"""
Tau_leadTkPtOverTauPt = NanoAODQuantity("Tau_leadTkPtOverTauPt")                                                                                                                          
"""dtype: Float_t; description: pt of the leading track divided by tau pt"""
Tau_mass = NanoAODQuantity("Tau_mass")                                                                                                                                                    
"""dtype: Float_t; description: mass"""
Tau_neutralIso = NanoAODQuantity("Tau_neutralIso")                                                                                                                                        
"""dtype: Float_t; description: neutral (photon) isolation"""
Tau_phi = NanoAODQuantity("Tau_phi")                                                                                                                                                      
"""dtype: Float_t; description: phi"""
Tau_photonsOutsideSignalCone = NanoAODQuantity("Tau_photonsOutsideSignalCone")                                                                                                            
"""dtype: Float_t; description: sum of photons outside signal cone"""
Tau_pt = NanoAODQuantity("Tau_pt")                                                                                                                                                        
"""dtype: Float_t; description: pt"""
Tau_puCorr = NanoAODQuantity("Tau_puCorr")                                                                                                                                                
"""dtype: Float_t; description: pileup correction"""
Tau_rawDeepTau2017v2p1VSe = NanoAODQuantity("Tau_rawDeepTau2017v2p1VSe")                                                                                                                  
"""dtype: Float_t; description: byDeepTau2017v2p1VSe raw output discriminator (deepTau2017v2p1)"""
Tau_rawDeepTau2017v2p1VSjet = NanoAODQuantity("Tau_rawDeepTau2017v2p1VSjet")                                                                                                              
"""dtype: Float_t; description: byDeepTau2017v2p1VSjet raw output discriminator (deepTau2017v2p1)"""
Tau_rawDeepTau2017v2p1VSmu = NanoAODQuantity("Tau_rawDeepTau2017v2p1VSmu")                                                                                                                
"""dtype: Float_t; description: byDeepTau2017v2p1VSmu raw output discriminator (deepTau2017v2p1)"""
Tau_rawIso = NanoAODQuantity("Tau_rawIso")                                                                                                                                                
"""dtype: Float_t; description: combined isolation (deltaBeta corrections)"""
Tau_rawIsodR03 = NanoAODQuantity("Tau_rawIsodR03")                                                                                                                                        
"""dtype: Float_t; description: combined isolation (deltaBeta corrections, dR=0.3)"""

TauEmbedding_SelectionNewMass = NanoAODQuantity("TauEmbedding_SelectionNewMass")                                                                                                          
"""dtype: Float_t; description: Mass of the Dimuon pair using the new selection algorithm (for internal studies only)"""
TauEmbedding_SelectionOldMass = NanoAODQuantity("TauEmbedding_SelectionOldMass")                                                                                                          
"""dtype: Float_t; description: Mass of the Dimuon pair using the old selection algorithm (for internal studies only)"""
TauEmbedding_initialMETEt = NanoAODQuantity("TauEmbedding_initialMETEt")                                                                                                                  
"""dtype: Float_t; description: MET Et of selected event"""
TauEmbedding_initialMETphi = NanoAODQuantity("TauEmbedding_initialMETphi")                                                                                                                
"""dtype: Float_t; description: MET phi of selected event"""
TauEmbedding_initialPuppiMETEt = NanoAODQuantity("TauEmbedding_initialPuppiMETEt")                                                                                                        
"""dtype: Float_t; description: PuppiMET Et of selected event"""
TauEmbedding_initialPuppiMETphi = NanoAODQuantity("TauEmbedding_initialPuppiMETphi")                                                                                                      
"""dtype: Float_t; description: PuppiMET phi of selected event"""
TauEmbedding_isMediumLeadingMuon = NanoAODQuantity("TauEmbedding_isMediumLeadingMuon")                                                                                                    
"""dtype: Bool_t; description: leading muon ID (medium)"""
TauEmbedding_isMediumTrailingMuon = NanoAODQuantity("TauEmbedding_isMediumTrailingMuon")                                                                                                  
"""dtype: Bool_t; description: trailing muon ID (medium)"""
TauEmbedding_isTightLeadingMuon = NanoAODQuantity("TauEmbedding_isTightLeadingMuon")                                                                                                      
"""dtype: Bool_t; description: leading muon ID (tight)"""
TauEmbedding_isTightTrailingMuon = NanoAODQuantity("TauEmbedding_isTightTrailingMuon")                                                                                                    
"""dtype: Bool_t; description: trailing muon ID (tight)"""
TauEmbedding_nInitialPairCandidates = NanoAODQuantity("TauEmbedding_nInitialPairCandidates")                                                                                              
"""dtype: Float_t; description: number of muons pairs suitable for selection (for internal studies only)"""

TkMET_phi = NanoAODQuantity("TkMET_phi")                                                                                                                                                  
"""dtype: Float_t; description: raw track MET phi"""
TkMET_pt = NanoAODQuantity("TkMET_pt")                                                                                                                                                    
"""dtype: Float_t; description: raw track MET pt"""
TkMET_sumEt = NanoAODQuantity("TkMET_sumEt")                                                                                                                                              
"""dtype: Float_t; description: raw track scalar sum of Et"""

nTrigObj = NanoAODQuantity("nTrigObj")                                                                                                                                                    
"""dtype: UInt_t; description: """
TrigObj_eta = NanoAODQuantity("TrigObj_eta")                                                                                                                                              
"""dtype: Float_t; description: eta"""
TrigObj_filterBits = NanoAODQuantity("TrigObj_filterBits")                                                                                                                                
"""dtype: Int_t; description: extra bits of associated information: 1 = CaloIdL_TrackIdL_IsoVL, 2 = 1e (WPTight), 4 = 1e (WPLoose), 8 = OverlapFilter PFTau, 16 = 2e, 32 = 1e-1mu, 64 = 1e-1tau, 128 = 3e, 256 = 2e-1mu, 512 = 1e-2mu, 1024 = 1e (32_L1DoubleEG_AND_L1SingleEGOr), 2048 = 1e (CaloIdVT_GsfTrkIdT), 4096 = 1e (PFJet), 8192 = 1e (Photon175_OR_Photon200) for Electron (PixelMatched e/gamma); 1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = 1mu, 16 = 2mu, 32 = 1mu-1e, 64 = 1mu-1tau, 128 = 3mu, 256 = 2mu-1e, 512 = 1mu-2e, 1024 = 1mu (Mu50), 2048 = 1mu (Mu100) for Muon; 1 = LooseChargedIso, 2 = MediumChargedIso, 4 = TightChargedIso, 8 = TightID OOSC photons, 16 = HPS, 32 = single-tau + tau+MET, 64 = di-tau, 128 = e-tau, 256 = mu-tau, 512 = VBF+di-tau for Tau; Jet bits: bit 0 for VBF cross-cleaned from loose iso PFTau, bit 1 for hltBTagCaloCSVp087Triple, bit 2 for hltDoubleCentralJet90, bit 3 for hltDoublePFCentralJetLooseID90, bit 4 for hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet, bit 5 for hltQuadCentralJet30, bit 6 for hltQuadPFCentralJetLooseID30, bit 7 for hltL1sQuadJetC50IorQuadJetC60IorHTT280IorHTT300IorHTT320IorTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBF or hltL1sQuadJetCIorTripleJetVBFIorHTT, bit 8 for hltQuadCentralJet45, bit 9 for hltQuadPFCentralJetLooseID45, bit 10  for hltL1sQuadJetC60IorHTT380IorHTT280QuadJetIorHTT300QuadJet or hltL1sQuadJetC50to60IorHTT280to500IorHTT250to340QuadJet bit 11 for hltBTagCaloCSVp05Double or hltBTagCaloDeepCSVp17Double, bit 12 for hltPFCentralJetLooseIDQuad30, bit 13 for hlt1PFCentralJetLooseID75, bit 14 for hlt2PFCentralJetLooseID60, bit 15 for hlt3PFCentralJetLooseID45, bit 16 for hlt4PFCentralJetLooseID40, bit 17 for hltBTagPFCSVp070Triple or hltBTagPFDeepCSVp24Triple or hltBTagPFDeepCSV4p5Triple  for Jet; HT bits: bit 0 for hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet, bit 1 for hltL1sQuadJetC50IorQuadJetC60IorHTT280IorHTT300IorHTT320IorTripleJet846848VBFIorTripleJet887256VBFIorTripleJet927664VBF or hltL1sQuadJetCIorTripleJetVBFIorHTT, bit 2 for hltL1sQuadJetC60IorHTT380IorHTT280QuadJetIorHTT300QuadJet or hltL1sQuadJetC50to60IorHTT280to500IorHTT250to340QuadJet, bit 3 for hltCaloQuadJet30HT300 or hltCaloQuadJet30HT320, bit 4 for hltPFCentralJetsLooseIDQuad30HT300 or hltPFCentralJetsLooseIDQuad30HT330 for HT; MHT bits: bit 0 for hltCaloQuadJet30HT300 or hltCaloQuadJet30HT320, bit 1 for hltPFCentralJetsLooseIDQuad30HT300 or hltPFCentralJetsLooseIDQuad30HT330 for MHT; """
TrigObj_id = NanoAODQuantity("TrigObj_id")                                                                                                                                                
"""dtype: Int_t; description: ID of the object: 11 = Electron (PixelMatched e/gamma), 22 = Photon (PixelMatch-vetoed e/gamma), 13 = Muon, 15 = Tau, 1 = Jet, 6 = FatJet, 2 = MET, 3 = HT, 4 = MHT"""
TrigObj_l1charge = NanoAODQuantity("TrigObj_l1charge")                                                                                                                                    
"""dtype: Int_t; description: charge of associated L1 seed"""
TrigObj_l1iso = NanoAODQuantity("TrigObj_l1iso")                                                                                                                                          
"""dtype: Int_t; description: iso of associated L1 seed"""
TrigObj_l1pt = NanoAODQuantity("TrigObj_l1pt")                                                                                                                                            
"""dtype: Float_t; description: pt of associated L1 seed"""
TrigObj_l1pt_2 = NanoAODQuantity("TrigObj_l1pt_2")                                                                                                                                        
"""dtype: Float_t; description: pt of associated secondary L1 seed"""
TrigObj_l2pt = NanoAODQuantity("TrigObj_l2pt")                                                                                                                                            
"""dtype: Float_t; description: pt of associated 'L2' seed (i.e. HLT before tracking/PF)"""
TrigObj_phi = NanoAODQuantity("TrigObj_phi")                                                                                                                                              
"""dtype: Float_t; description: phi"""
TrigObj_pt = NanoAODQuantity("TrigObj_pt")                                                                                                                                                
"""dtype: Float_t; description: pt"""

HLTriggerFinalPath = NanoAODQuantity("HLTriggerFinalPath")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
HLTriggerFirstPath = NanoAODQuantity("HLTriggerFirstPath")                                                                                                                                
"""dtype: Bool_t; description: Trigger/flag bit (process: HLT)"""
LHEPdfWeight = NanoAODQuantity("LHEPdfWeight")                                                                                                                                            
"""dtype: Float_t; description: LHE pdf variation weights (w_var / w_nominal) for LHA IDs 325300 - 325402"""
LHEReweightingWeight = NanoAODQuantity("LHEReweightingWeight")                                                                                                                            
"""dtype: Float_t; description: """
LHEScaleWeight = NanoAODQuantity("LHEScaleWeight")                                                                                                                                        
"""dtype: Float_t; description: LHE scale variation weights (w_var / w_nominal); [0] is MUF="0.5" MUR="0.5"; [1] is MUF="1.0" MUR="0.5"; [2] is MUF="2.0" MUR="0.5"; [3] is MUF="0.5" MUR="1.0"; [4] is MUF="1.0" MUR="1.0"; [5] is MUF="2.0" MUR="1.0"; [6] is MUF="0.5" MUR="2.0"; [7] is MUF="1.0" MUR="2.0"; [8] is MUF="2.0" MUR="2.0" """
PSWeight = NanoAODQuantity("PSWeight")                                                                                                                                                    
"""dtype: Float_t; description: PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;"""
SoftActivityJetHT = NanoAODQuantity("SoftActivityJetHT")                                                                                                                                  
"""dtype: Float_t; description: scalar sum of soft activity jet pt, pt>1"""
SoftActivityJetHT10 = NanoAODQuantity("SoftActivityJetHT10")                                                                                                                              
"""dtype: Float_t; description: scalar sum of soft activity jet pt , pt >10"""
SoftActivityJetHT2 = NanoAODQuantity("SoftActivityJetHT2")                                                                                                                                
"""dtype: Float_t; description: scalar sum of soft activity jet pt, pt >2"""
SoftActivityJetHT5 = NanoAODQuantity("SoftActivityJetHT5")                                                                                                                                
"""dtype: Float_t; description: scalar sum of soft activity jet pt, pt>5"""
SoftActivityJetNjets10 = NanoAODQuantity("SoftActivityJetNjets10")                                                                                                                        
"""dtype: Int_t; description: number of soft activity jet pt, pt >2"""
SoftActivityJetNjets2 = NanoAODQuantity("SoftActivityJetNjets2")                                                                                                                          
"""dtype: Int_t; description: number of soft activity jet pt, pt >10"""
SoftActivityJetNjets5 = NanoAODQuantity("SoftActivityJetNjets5")                                                                                                                          
"""dtype: Int_t; description: number of soft activity jet pt, pt >5"""
event = NanoAODQuantity("event")                                                                                                                                                          
"""dtype: ULong64_t; description: event/l"""
fixedGridRhoFastjetAll = NanoAODQuantity("fixedGridRhoFastjetAll")                                                                                                                        
"""dtype: Float_t; description: rho from all PF Candidates, used e.g. for JECs"""
fixedGridRhoFastjetCentral = NanoAODQuantity("fixedGridRhoFastjetCentral")                                                                                                                
"""dtype: Float_t; description: rho from all PF Candidates for central region, used e.g. for JECs"""
fixedGridRhoFastjetCentralCalo = NanoAODQuantity("fixedGridRhoFastjetCentralCalo")                                                                                                        
"""dtype: Float_t; description: rho from calo towers with |eta| < 2.5, used e.g. egamma PFCluster isolation"""
fixedGridRhoFastjetCentralChargedPileUp = NanoAODQuantity("fixedGridRhoFastjetCentralChargedPileUp")                                                                                      
"""dtype: Float_t; description: rho from charged PF Candidates for central region, used e.g. for JECs"""
fixedGridRhoFastjetCentralNeutral = NanoAODQuantity("fixedGridRhoFastjetCentralNeutral")                                                                                                  
"""dtype: Float_t; description: rho from neutral PF Candidates with |eta| < 2.5, used e.g. for rho corrections of some lepton isolations"""
genTtbarId = NanoAODQuantity("genTtbarId")                                                                                                                                                
"""dtype: Int_t; description: ttbar categorization"""
genWeight = NanoAODQuantity("genWeight")                                                                                                                                                  
"""dtype: Float_t; description: generator weight"""
luminosityBlock = NanoAODQuantity("luminosityBlock")                                                                                                                                      
"""dtype: UInt_t; description: luminosityBlock/i"""
run = NanoAODQuantity("run")                                                                                                                                                              
"""dtype: UInt_t; description: run/i"""

nboostedTau = NanoAODQuantity("nboostedTau")                                                                                                                                              
"""dtype: UInt_t; description: slimmedBoostedTaus after basic selection (pt > 40 && tauID('decayModeFindingNewDMs') && (tauID('byVVLooseIsolationMVArun2017v2DBoldDMwLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017') || tauID('byVVLooseIsolationMVArun2017v2DBnewDMwLT2017')))"""
boostedTau_charge = NanoAODQuantity("boostedTau_charge")                                                                                                                                  
"""dtype: Int_t; description: electric charge"""
boostedTau_chargedIso = NanoAODQuantity("boostedTau_chargedIso")                                                                                                                          
"""dtype: Float_t; description: charged isolation"""
boostedTau_decayMode = NanoAODQuantity("boostedTau_decayMode")                                                                                                                            
"""dtype: Int_t; description: decayMode()"""
boostedTau_eta = NanoAODQuantity("boostedTau_eta")                                                                                                                                        
"""dtype: Float_t; description: eta"""
boostedTau_genPartFlav = NanoAODQuantity("boostedTau_genPartFlav")                                                                                                                        
"""dtype: UChar_t; description: Flavour of genParticle for MC matching to status==2 taus: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched"""
boostedTau_genPartIdx = NanoAODQuantity("boostedTau_genPartIdx")                                                                                                                          
"""dtype: Int_t; description: Index into genParticle list for MC matching to status==2 taus"""
boostedTau_idAntiEle2018 = NanoAODQuantity("boostedTau_idAntiEle2018")                                                                                                                    
"""dtype: UChar_t; description: Anti-electron MVA discriminator V6 (2018): bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight"""
boostedTau_idAntiMu = NanoAODQuantity("boostedTau_idAntiMu")                                                                                                                              
"""dtype: UChar_t; description: Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight"""
boostedTau_idMVAnewDM2017v2 = NanoAODQuantity("boostedTau_idMVAnewDM2017v2")                                                                                                              
"""dtype: UChar_t; description: IsolationMVArun2017v2DBnewDMwLT ID working point (2017v2): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight"""
boostedTau_idMVAoldDM2017v2 = NanoAODQuantity("boostedTau_idMVAoldDM2017v2")                                                                                                              
"""dtype: UChar_t; description: IsolationMVArun2017v2DBoldDMwLT ID working point (2017v2): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight"""
boostedTau_idMVAoldDMdR032017v2 = NanoAODQuantity("boostedTau_idMVAoldDMdR032017v2")                                                                                                      
"""dtype: UChar_t; description: IsolationMVArun2017v2DBoldDMdR0p3wLT ID working point (2017v2): bitmask 1 = VVLoose, 2 = VLoose, 4 = Loose, 8 = Medium, 16 = Tight, 32 = VTight, 64 = VVTight"""
boostedTau_jetIdx = NanoAODQuantity("boostedTau_jetIdx")                                                                                                                                  
"""dtype: Int_t; description: index of the associated jet (-1 if none)"""
boostedTau_leadTkDeltaEta = NanoAODQuantity("boostedTau_leadTkDeltaEta")                                                                                                                  
"""dtype: Float_t; description: eta of the leading track, minus tau eta"""
boostedTau_leadTkDeltaPhi = NanoAODQuantity("boostedTau_leadTkDeltaPhi")                                                                                                                  
"""dtype: Float_t; description: phi of the leading track, minus tau phi"""
boostedTau_leadTkPtOverTauPt = NanoAODQuantity("boostedTau_leadTkPtOverTauPt")                                                                                                            
"""dtype: Float_t; description: pt of the leading track divided by tau pt"""
boostedTau_mass = NanoAODQuantity("boostedTau_mass")                                                                                                                                      
"""dtype: Float_t; description: mass"""
boostedTau_neutralIso = NanoAODQuantity("boostedTau_neutralIso")                                                                                                                          
"""dtype: Float_t; description: neutral (photon) isolation"""
boostedTau_phi = NanoAODQuantity("boostedTau_phi")                                                                                                                                        
"""dtype: Float_t; description: phi"""
boostedTau_photonsOutsideSignalCone = NanoAODQuantity("boostedTau_photonsOutsideSignalCone")                                                                                              
"""dtype: Float_t; description: sum of photons outside signal cone"""
boostedTau_pt = NanoAODQuantity("boostedTau_pt")                                                                                                                                          
"""dtype: Float_t; description: pt"""
boostedTau_puCorr = NanoAODQuantity("boostedTau_puCorr")                                                                                                                                  
"""dtype: Float_t; description: pileup correction"""
boostedTau_rawAntiEle2018 = NanoAODQuantity("boostedTau_rawAntiEle2018")                                                                                                                  
"""dtype: Float_t; description: Anti-electron MVA discriminator V6 raw output discriminator (2018)"""
boostedTau_rawAntiEleCat2018 = NanoAODQuantity("boostedTau_rawAntiEleCat2018")                                                                                                            
"""dtype: Int_t; description: Anti-electron MVA discriminator V6 category (2018)"""
boostedTau_rawIso = NanoAODQuantity("boostedTau_rawIso")                                                                                                                                  
"""dtype: Float_t; description: combined isolation (deltaBeta corrections)"""
boostedTau_rawIsodR03 = NanoAODQuantity("boostedTau_rawIsodR03")                                                                                                                          
"""dtype: Float_t; description: combined isolation (deltaBeta corrections, dR=0.3)"""
boostedTau_rawMVAnewDM2017v2 = NanoAODQuantity("boostedTau_rawMVAnewDM2017v2")                                                                                                            
"""dtype: Float_t; description: byIsolationMVArun2017v2DBnewDMwLT raw output discriminator (2017v2)"""
boostedTau_rawMVAoldDM2017v2 = NanoAODQuantity("boostedTau_rawMVAoldDM2017v2")                                                                                                            
"""dtype: Float_t; description: byIsolationMVArun2017v2DBoldDMwLT raw output discriminator (2017v2)"""
boostedTau_rawMVAoldDMdR032017v2 = NanoAODQuantity("boostedTau_rawMVAoldDMdR032017v2")                                                                                                    
"""dtype: Float_t; description: byIsolationMVArun2017v2DBoldDMdR0p3wLT raw output discriminator (2017v2)"""

btagWeight_CSVV2 = NanoAODQuantity("btagWeight_CSVV2")                                                                                                                                    
"""dtype: Float_t; description: b-tag event weight for CSVV2"""
btagWeight_DeepCSVB = NanoAODQuantity("btagWeight_DeepCSVB")                                                                                                                              
"""dtype: Float_t; description: b-tag event weight for DeepCSVB"""

nLHEPdfWeight = NanoAODQuantity("nLHEPdfWeight")                                                                                                                                          
"""dtype: UInt_t; description: """
nLHEReweightingWeight = NanoAODQuantity("nLHEReweightingWeight")                                                                                                                          
"""dtype: UInt_t; description: """
nLHEScaleWeight = NanoAODQuantity("nLHEScaleWeight")                                                                                                                                      
"""dtype: UInt_t; description: """
nPSWeight = NanoAODQuantity("nPSWeight")                                                                                                                                                  
"""dtype: UInt_t; description: """
