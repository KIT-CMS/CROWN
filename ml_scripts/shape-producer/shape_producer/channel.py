# -*- coding: utf-8 -*-

from .cutstring import Cuts, Cut
import logging
logger = logging.getLogger(__name__)
"""
"""


class Channel(object):
    @property
    def cuts(self):
        return self._cuts

    @property
    def name(self):
        return self._name


###########################################
# Common block
###########################################
class EE(Channel):
    def __init__(self):
        self._name = "ee"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.1 && iso_2<0.1", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_singleelectron==1 && pt_1>26 && pt_2>26)",
                "trg_singleelectron"))


# Common MM
class MM2016(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(pt_1 > 23 && trg_singlemuon==1)&&(0<1)", "trg_selection"))


class MM2017(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(trg_singlemuon_27==1 || trg_singlemuon_24==1)", "trg_selection"))


class MM2018(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(trg_singlemuon_27==1 || trg_singlemuon_24==1)", "trg_selection"))


# Common MT
class MT2017(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byTightDeepTau2017v2p1VSmu_2>0.5", "againstMuonDiscriminator"),
            Cut("byVVLooseDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((trg_singlemuon_27 == 1) || (trg_singlemuon_24 == 1) || (pt_1 < 25 && trg_crossmuon_mu20tau27 == 1))",
                "trg_selection"))


class MT2018(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byTightDeepTau2017v2p1VSmu_2>0.5", "againstMuonDiscriminator"),
            Cut("byVVLooseDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && (((trg_singlemuon_27 == 1) || (trg_singlemuon_24 == 1)) || (pt_1 < 25 && (trg_crossmuon_mu20tau27_hps == 1 || trg_crossmuon_mu20tau27 == 1)))",
                "trg_selection"))


class MT2016(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter==1", "met_filter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byTightDeepTau2017v2p1VSmu_2 > 0.5", "againstMuonDiscriminator"),
            Cut("byVVLooseDeepTau2017v2p1VSe_2>0.5","againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((pt_1 >= 23 && trg_singlemuon == 1) || (trg_mutaucross == 1 && pt_1 < 23 && abs(eta_2)<2.1))","trg_selection")
        )


# Common ET
class ET2017(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byVLooseDeepTau2017v2p1VSmu_2>0.5", "againstMuonDiscriminator"),
            Cut("byTightDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && pt_1 > 25 && ((((trg_singleelectron_35 == 1) || (trg_singleelectron_32 == 1) || ((trg_singleelectron_27 == 1))) || (abs(eta_1)>1.5 && pt_1 >= 28 && pt_1 < 40 && isEmbedded)) || (pt_1>25 && pt_1<28 && pt_2>35 && ((isEmbedded && (abs(eta_1)>1.5)) || (trg_crossele_ele24tau30 == 1))))",
                "trg_selection"))


class ET2018(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byVLooseDeepTau2017v2p1VSmu_2>0.5", "againstMuonDiscriminator"),
            Cut("byTightDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && (((trg_singleelectron_35 == 1) || (trg_singleelectron_32 == 1) || (pt_1>25 && pt_1<33 && pt_2>35 && (trg_crossele_ele24tau30_hps == 1 || trg_crossele_ele24tau30 == 1))))",
                "trg_selection"))


class ET2016(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter==1", "met_filter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byVLooseDeepTau2017v2p1VSmu_2>0.5", "againstMuonDiscriminator"),
            Cut("byTightDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "ele_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((pt_1>26 && (trg_singleelectron==1)) || (pt_1<26 && pt_1>25 && (trg_eletaucross==1)))", "trg_selection"))


# Common TT
class TT2016(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter==1", "met_filter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byVLooseDeepTau2017v2p1VSmu_1>0.5 && byVLooseDeepTau2017v2p1VSmu_2>0.5", "againstMuonDiscriminator"),
            Cut("byVVLooseDeepTau2017v2p1VSe_1>0.5 && byVVLooseDeepTau2017v2p1VSe_2>0.5", "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_1>0.5", "tau_1_iso"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5", "tau_2_iso"),
            Cut("q_1*q_2<0", "os"),
            # Cut("pt_tt>50", "pt_h"),
            Cut("trg_doubletau==1", "trg_doubletau"))


class TT2017(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byVLooseDeepTau2017v2p1VSmu_1>0.5 && byVLooseDeepTau2017v2p1VSmu_2>0.5",
                "againstMuonDiscriminator"),
            Cut("byVVLooseDeepTau2017v2p1VSe_1>0.5 && byVVLooseDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_1>0.5",
                "tau_1_iso"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5",
                "tau_2_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_doubletau_35_tightiso_tightid == 1) || (trg_doubletau_40_mediso_tightid == 1) || (trg_doubletau_40_tightiso == 1)",
                "trg_selection"))


class TT2018(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("byVLooseDeepTau2017v2p1VSmu_1>0.5 && byVLooseDeepTau2017v2p1VSmu_2>0.5",
                "againstMuonDiscriminator"),
            Cut("byVVLooseDeepTau2017v2p1VSe_1>0.5 && byVVLooseDeepTau2017v2p1VSe_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byMediumDeepTau2017v2p1VSjet_1>0.5",
                "tau_1_iso"),
            Cut("byMediumDeepTau2017v2p1VSjet_2>0.5",
                "tau_2_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(((!(isMC||isEmbedded) && run>=317509) || (isMC||isEmbedded)) && (trg_doubletau_35_mediso_hps == 1)) || (!(isMC||isEmbedded) && (run<317509) && ((trg_doubletau_35_tightiso_tightid == 1) || (trg_doubletau_40_mediso_tightid == 1) || (trg_doubletau_40_tightiso == 1)))",
                "trg_selection"))


# Common EM
class EM2016(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("abs(eta_1)<2.4", "electron_eta"),
            Cut("pt_1>15 && pt_2>15 && ((pt_1>15 && pt_2>24 && trg_muonelectron_mu23ele12 == 1) || (pt_1>24 && pt_2>15 && trg_muonelectron_mu8ele23 == 1))","trg_selection"))


class EM2017(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("abs(eta_1)<2.4","electron_eta"),
            Cut("pt_1>15 && pt_2>15 && ((trg_muonelectron_mu23ele12 == 1) || (trg_muonelectron_mu8ele23 == 1))",
                "trg_selection"))


class EM2018(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("(nbtag>0) *(kinfit_chi2>-1) *(kinfit_chi2<990)", "btag"),
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("abs(eta_1)<2.4", "electron_eta"),
            Cut("(trg_muonelectron_mu23ele12 == 1 && pt_1>15 && pt_2 > 24) || (trg_muonelectron_mu8ele23 == 1 && pt_1>24 && pt_2>15)",
                "trg_selection"))


###########################################
# SM block
###########################################
# SM 2016
class ETSM2016(ET2016):
    def __init__(self, **kvargs):
        super(ETSM2016, self).__init__(**kvargs)


class MTSM2016(MT2016):
    def __init__(self, **kvargs):
        super(MTSM2016, self).__init__(**kvargs)


class TTSM2016(TT2016):
    def __init__(self, **kvargs):
        super(TTSM2016, self).__init__(**kvargs)


class EMSM2016(EM2016):
    def __init__(self, **kvargs):
        super(EMSM2016, self).__init__(**kvargs)

class MMSM2016(MM2016):
    def __init__(self, **kvargs):
        super(MMSM2016, self).__init__(**kvargs)


class EESM2016(ET2016):
    def __init__(self, **kvargs):
        super(EESM2016, self).__init__(**kvargs)


# SM 2017
class ETSM2017(ET2017):
    def __init__(self, **kvargs):
        super(ETSM2017, self).__init__(**kvargs)


class MTSM2017(MT2017):
    def __init__(self, **kvargs):
        super(MTSM2017, self).__init__(**kvargs)


class TTSM2017(TT2017):
    def __init__(self, **kvargs):
        super(TTSM2017, self).__init__(**kvargs)


class EMSM2017(EM2017):
    def __init__(self, **kvargs):
        super(EMSM2017, self).__init__(**kvargs)

class MMSM2017(MM2017):
    def __init__(self, **kvargs):
        super(MMSM2017, self).__init__(**kvargs)


class EESM2017(ET2017):
    def __init__(self, **kvargs):
        super(EESM2017, self).__init__(**kvargs)


# SM 2018
class ETSM2018(ET2018):
    def __init__(self, **kvargs):
        super(ETSM2018, self).__init__(**kvargs)


class MTSM2018(MT2018):
    def __init__(self, **kvargs):
        super(MTSM2018, self).__init__(**kvargs)


class TTSM2018(TT2018):
    def __init__(self, **kvargs):
        super(TTSM2018, self).__init__(**kvargs)


class EMSM2018(EM2018):
    def __init__(self, **kvargs):
        super(EMSM2018, self).__init__(**kvargs)

class MMSM2018(MM2018):
    def __init__(self, **kvargs):
        super(MMSM2018, self).__init__(**kvargs)


class EESM2018(ET2018):
    def __init__(self, **kvargs):
        super(EESM2018, self).__init__(**kvargs)


###########################################
# MSSM block
###########################################
# MSSM 2016
class ETMSSM2016(ET2016):
    def __init__(self, **kvargs):
        super(ETMSSM2016, self).__init__(**kvargs)


class MTMSSM2016(MT2016):
    def __init__(self, **kvargs):
        super(MTMSSM2016, self).__init__(**kvargs)


class TTMSSM2016(TT2016):
    def __init__(self, **kvargs):
        super(TTMSSM2016, self).__init__(**kvargs)


class EMMSSM2016(EM2016):
    def __init__(self, **kvargs):
        super(EMMSSM2016, self).__init__(**kvargs)


# MSSM 2017
class ETMSSM2017(ET2017):
    def __init__(self, **kvargs):
        super(ETMSSM2017, self).__init__(**kvargs)


class MTMSSM2017(MT2017):
    def __init__(self, **kvargs):
        super(MTMSSM2017, self).__init__(**kvargs)


class TTMSSM2017(TT2017):
    def __init__(self, **kvargs):
        super(TTMSSM2017, self).__init__(**kvargs)


class EMMSSM2017(EM2017):
    def __init__(self, **kvargs):
        super(EMMSSM2017, self).__init__(**kvargs)


# MSSM 2018
class ETMSSM2018(ET2018):
    def __init__(self, **kvargs):
        super(ETMSSM2018, self).__init__(**kvargs)


class MTMSSM2018(MT2018):
    def __init__(self, **kvargs):
        super(MTMSSM2018, self).__init__(**kvargs)


class TTMSSM2018(TT2018):
    def __init__(self, **kvargs):
        super(TTMSSM2018, self).__init__(**kvargs)


class EMMSSM2018(EM2018):
    def __init__(self, **kvargs):
        super(EMMSSM2018, self).__init__(**kvargs)


# PU - no cuts
class PU(Channel):
    def __init__(self):
        self._name = "pu"
        self._cuts = Cuts()


# collection of channels an analysis can be ran on
class Channels(object):
    def __init__(self, name):
        self._name = name
        self._channels = []

    def add(self, channel):
        self._channels.append(channel)

    @property
    def name(self):
        return self._name
