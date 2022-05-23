# -*- coding: utf-8 -*-

import copy
import os

from .estimation_methods import EstimationMethod, SStoOSEstimationMethod, ABCDEstimationMethod, SumUpEstimationMethod, NewFakeEstimationMethodLT, NewFakeEstimationMethodTT
from .histogram import *
from .cutstring import *
from .process import *
from .systematics import *
from .systematic_variations import *
from .era import log_query

import logging
logger = logging.getLogger(__name__)

ggH_htxs = {
    "ggH125": "(htxs_stage1p1cat>=100)&&(htxs_stage1p1cat<=113)",
}

qqH_htxs = {
    "qqH125": "(htxs_stage1p1cat>=200)&&(htxs_stage1p1cat<=210)",
}


def get_triggerweight_for_channel(channel):
    weight = Weight("1.0", "triggerweight")

    singleMC = "singleTriggerMCEfficiencyWeightKIT_1"
    crossMCL = "crossTriggerMCEfficiencyWeightKIT_1"
    # MCTau_1 = "((abs(eta_2)<2.1)*((byMediumDeepTau2017v2p1VSjet_1<0.5 && byVVVLooseDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_vloose_DeepTau_1 + (byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_medium_DeepTau_1))"
    MCTau_1 = "((abs(eta_2)<2.1)*((byMediumDeepTau2017v2p1VSjet_1<0.5 && byVVVLooseDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_vloose_DeepTau_1 + (byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_medium_DeepTau_1))"  # Hotfix for old trigger weights.
    MCTau_2 = MCTau_1.replace("_1", "_2")

    if "mt" in channel:
        trig_sL = "(pt_1 >= 23 && trg_singlemuon)"
        trig_X = "(pt_1 < 23 && trg_mutaucross && abs(eta_2)<2.1)"

        MuTauMC = "*".join([trig_sL, singleMC]) + "+" + "*".join([trig_X, crossMCL, MCTau_2])
        MuTauData = MuTauMC.replace("MC", "Data")
        MuTau = "(" + MuTauData + ")/(" + MuTauMC + ")"
        weight = Weight(MuTau, "triggerweight")

    elif "et" in channel:
        trig_sL = "(trg_singleelectron)"
        trig_X = "(pt_1 > 25 && pt_1 < 26 && trg_eletaucross)"

        ElTauMC = "*".join([trig_sL, singleMC]) + "+" + "*".join([trig_X, crossMCL, MCTau_2])
        ElTauData = ElTauMC.replace("MC", "Data")
        ElTau = "(" + ElTauData + ")/(" + ElTauMC + ")"
        weight = Weight(ElTau, "triggerweight")

    elif "tt" in channel:
        DiTauMC = "*".join([MCTau_1, MCTau_2])
        DiTauData = DiTauMC.replace("MC", "Data")
        DiTau = "("+DiTauData+")/("+DiTauMC+")"
        weight = Weight(DiTau, "triggerweight")

    elif "em" in channel:
        ElMuData = "(trigger_23_data_Weight_2*trigger_12_data_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_data_Weight_1*trigger_8_data_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_data_Weight_2*trigger_23_data_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))"
        ElMuEmb = ElMuData.replace('data', 'mc')
        ElMu = "("+ElMuData+")/("+ElMuEmb+")"
        weight = Weight(ElMu, "triggerweight")
    return weight


def aiso_muon_correction(channel):
    if "em" in channel:
        return Weight("(iso_2 <= 0.15)*1.0+((iso_2 > 0.15 && iso_2 < 0.20)*(((abs(eta_2) > 0 && abs(eta_2) < 0.9)*(((pt_2 > 10.0 && pt_2 < 15.0)*0.9327831343)+((pt_2 > 15.0 && pt_2 < 20.0)*0.953716709495)+((pt_2 > 20.0 && pt_2 < 22.0)*0.970107931338)+((pt_2 > 22.0 && pt_2 < 24.0)*0.975194707069)+((pt_2 > 24.0 && pt_2 < 26.0)*0.987926954702)+((pt_2 > 26.0 && pt_2 < 28.0)*0.984899759186)+((pt_2 > 28.0 && pt_2 < 30.0)*0.98461236424)+((pt_2 > 30.0 && pt_2 < 32.0)*0.982678110647)+((pt_2 > 32.0 && pt_2 < 34.0)*0.969564600424)+((pt_2 > 34.0 && pt_2 < 36.0)*0.952545963996)+((pt_2 > 36.0 && pt_2 < 38.0)*0.93801027736)+((pt_2 > 38.0 && pt_2 < 40.0)*0.928180181431)+((pt_2 > 40.0 && pt_2 < 45.0)*0.936587036212)+((pt_2 > 45.0 && pt_2 < 50.0)*0.937996645301)+((pt_2 > 50.0 && pt_2 < 60.0)*0.906087339587)+((pt_2 > 60.0 && pt_2 < 80.0)*0.887611582681)+((pt_2 > 80.0 && pt_2 < 100.0)*0.8834356199)+((pt_2 > 100.0 && pt_2 < 200.0)*1.17717912783)+((pt_2 > 200.0)*0.782070939337)))+((abs(eta_2) > 0.9 && abs(eta_2) < 1.2)*(((pt_2 > 10.0 && pt_2 < 15.0)*0.886291162409)+((pt_2 > 15.0 && pt_2 < 20.0)*0.915805873893)+((pt_2 > 20.0 && pt_2 < 22.0)*0.928213984488)+((pt_2 > 22.0 && pt_2 < 24.0)*0.968808738285)+((pt_2 > 24.0 && pt_2 < 26.0)*1.00847685497)+((pt_2 > 26.0 && pt_2 < 28.0)*1.01823133239)+((pt_2 > 28.0 && pt_2 < 30.0)*0.992528525978)+((pt_2 > 30.0 && pt_2 < 32.0)*0.978795541905)+((pt_2 > 32.0 && pt_2 < 34.0)*0.942964386045)+((pt_2 > 34.0 && pt_2 < 36.0)*0.938710844744)+((pt_2 > 36.0 && pt_2 < 38.0)*0.922702562159)+((pt_2 > 38.0 && pt_2 < 40.0)*0.897758415445)+((pt_2 > 40.0 && pt_2 < 45.0)*0.909162491)+((pt_2 > 45.0 && pt_2 < 50.0)*0.90265167858)+((pt_2 > 50.0 && pt_2 < 60.0)*0.912325787246)+((pt_2 > 60.0 && pt_2 < 80.0)*0.897018870572)+((pt_2 > 80.0 && pt_2 < 100.0)*0.972647372742)+((pt_2 > 100.0 && pt_2 < 200.0)*1.38562213225)+((pt_2 > 200.0)*0.738304282781)))+((abs(eta_2) > 1.2 && abs(eta_2) < 2.1)*(((pt_2 > 10.0 && pt_2 < 15.0)*0.88678133381)+((pt_2 > 15.0 && pt_2 < 20.0)*0.855042730357)+((pt_2 > 20.0 && pt_2 < 22.0)*0.897842682768)+((pt_2 > 22.0 && pt_2 < 24.0)*0.905849165918)+((pt_2 > 24.0 && pt_2 < 26.0)*0.910626040493)+((pt_2 > 26.0 && pt_2 < 28.0)*0.952076550342)+((pt_2 > 28.0 && pt_2 < 30.0)*0.968869527514)+((pt_2 > 30.0 && pt_2 < 32.0)*0.942569376345)+((pt_2 > 32.0 && pt_2 < 34.0)*0.941205386066)+((pt_2 > 34.0 && pt_2 < 36.0)*0.925500794627)+((pt_2 > 36.0 && pt_2 < 38.0)*0.907300484346)+((pt_2 > 38.0 && pt_2 < 40.0)*0.87984390364)+((pt_2 > 40.0 && pt_2 < 45.0)*0.87339713294)+((pt_2 > 45.0 && pt_2 < 50.0)*0.87980130335)+((pt_2 > 50.0 && pt_2 < 60.0)*0.860066115116)+((pt_2 > 60.0 && pt_2 < 80.0)*0.857712710727)+((pt_2 > 80.0 && pt_2 < 100.0)*1.0645948221)+((pt_2 > 100.0 && pt_2 < 200.0)*1.18849162977)+((pt_2 > 200.0)*1.28784467602)))+((abs(eta_2) > 2.1 && abs(eta_2) < 2.4)*(((pt_2 > 10.0 && pt_2 < 15.0)*0.776167258269)+((pt_2 > 15.0 && pt_2 < 20.0)*0.770868349402)+((pt_2 > 20.0 && pt_2 < 22.0)*0.808779663589)+((pt_2 > 22.0 && pt_2 < 24.0)*0.812754474056)+((pt_2 > 24.0 && pt_2 < 26.0)*0.84667222665)+((pt_2 > 26.0 && pt_2 < 28.0)*0.837142139899)+((pt_2 > 28.0 && pt_2 < 30.0)*0.8356560823)+((pt_2 > 30.0 && pt_2 < 32.0)*0.888386540505)+((pt_2 > 32.0 && pt_2 < 34.0)*0.881083238091)+((pt_2 > 34.0 && pt_2 < 36.0)*0.872500048844)+((pt_2 > 36.0 && pt_2 < 38.0)*0.861737355714)+((pt_2 > 38.0 && pt_2 < 40.0)*0.872186406375)+((pt_2 > 40.0 && pt_2 < 45.0)*0.853060222605)+((pt_2 > 45.0 && pt_2 < 50.0)*0.927735148085)+((pt_2 > 50.0 && pt_2 < 60.0)*0.82749753618)+((pt_2 > 60.0 && pt_2 < 80.0)*0.924329437022)+((pt_2 > 80.0 && pt_2 < 100.0)*0.887073323216)+((pt_2 > 100.0 && pt_2 < 200.0)*1.15559449916)+((pt_2 > 200.0)*0.37887229649)))))", "m_aiso_correction")
    else:
        return  Weight("1.0", "m_aiso_correction")


def get_singlelepton_triggerweight_for_channel(channel):
    weight = Weight("1.0", "triggerweight_sl")

    # MCTau_1 = "((byMediumDeepTau2017v2p1VSjet_1<0.5 && byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_medium_DeepTau_1 + (byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_medium_DeepTau_1)"
    MCTau_1 = "((byMediumDeepTau2017v2p1VSjet_1<0.5 && byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_medium_DeepTau_1 + (byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerMCEfficiencyWeight_medium_DeepTau_1)"
    MCTau_2 = MCTau_1.replace("_1", "_2")

    # if "mt" in channel or "et" in channel:
    #     weight = Weight("singleTriggerDataEfficiencyWeightKIT_1/singleTriggerMCEfficiencyWeightKIT_1","triggerweight")
    # elif "tt" in channel:
    #     DiTauMC = "*".join([MCTau_1,MCTau_2])
    #     DiTauData = DiTauMC.replace("MC","Data")
    #     DiTau = "("+DiTauData+")/("+DiTauMC+")"
    #     weight = Weight(DiTau,"triggerweight")

    return weight


def get_tauByIsoIdWeight_for_channel(channel):
    # WPs: Tight 0.87, urrently used: SR mt,et Tight; SR tt Tight
    # Source: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_ID_efficiency
    weight = Weight("1.0", "taubyIsoIdWeight")
    if "ID" in channel.__class__.__name__:  # this is used for the TauID measurements
        return weight
    elif "mt" in channel.name or "et" in channel.name:
        weight = Weight("((gen_match_2 == 5)*tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2 + (gen_match_2 != 5))",
                        "taubyIsoIdWeight")
    elif "tt" in channel.name:
        # dm11_nom = 0.89484048
        # weight once dm11 is fixed:
        weight = Weight("((gen_match_1 == 5)*tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_1 + (gen_match_1 != 5))*((gen_match_2 == 5)*tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2 + (gen_match_2 != 5))", "taubyIsoIdWeight")
        # weight = Weight("(((gen_match_1 == 5)*(((decayMode_1!=11)*tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_1)+((decayMode_1==11)*{dm11_nom})) + (gen_match_1 != 5))*((gen_match_2 == 5)*(((decayMode_2!=11)*tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2)+((decayMode_2==11)*{dm11_nom})) + (gen_match_2 != 5)))".format(dm11_nom=dm11_nom), "taubyIsoIdWeight")
    return weight

def get_eleRecoWeight_for_channel(channel):
    weight = Weight("1.0", "eleRecoWeight")
    if "et" in channel:
        weight = Weight("(eleRecoWeight_1)", "eleRecoWeight")
    if "em" in channel:
        weight = Weight("(eleRecoWeight_1)",
                        "eleRecoWeight")
    return weight


class DataEstimation(EstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 friend_directory=None,
                 folder="nominal"):
        super(DataEstimation, self).__init__(name="data_obs",
                                             folder=folder,
                                             era=era,
                                             directory=directory,
                                             friend_directory=friend_directory,
                                             channel=channel,
                                             mc_campaign=None)
        self._channel = channel

    def get_files(self):
        return self.artus_file_names(self.era.data_files(self._channel))

    def get_cuts(self):
        return Cuts()


class FakeEstimationLT(DataEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(DataEstimation, self).__init__(
            name="fakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(Weight("ff2_nom", "fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5",
                "tau_aiso"))
        return super(FakeEstimationLT,
                     self).create_root_objects(aiso_systematic)


class NewFakeEstimationLT(NewFakeEstimationMethodLT):
    def __init__(
            self,
            era,
            directory,
            channel,
            nofake_processes,
            data_process,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(NewFakeEstimationLT, self).__init__(
            name="jetFakes",
            folder="nominal",
            era=era,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            nofake_processes=nofake_processes,
            data_process=data_process,
            aisoCut=Cut(
                "byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5",
                "tau_aiso"),
            fakeWeightstring="ff2_nom")


class NewFakeEstimationTT(NewFakeEstimationMethodTT):
    def __init__(
            self,
            era,
            directory,
            channel,
            nofake_processes,
            data_process,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(NewFakeEstimationTT, self).__init__(
            name="jetFakes",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            nofake_processes=nofake_processes,
            data_process=data_process,
            aisoCut=Cut(
                "(byMediumDeepTau2017v2p1VSjet_2>0.5&&byMediumDeepTau2017v2p1VSjet_1<0.5&&byVVVLooseDeepTau2017v2p1VSjet_1>0.5)||(byMediumDeepTau2017v2p1VSjet_1>0.5&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5)",
                "tau_aiso"),
            fakeWeightstring="(0.5*ff1_nom*(byMediumDeepTau2017v2p1VSjet_1<0.5)+0.5*ff2_nom*(byMediumDeepTau2017v2p1VSjet_2<0.5))"
        )

class FakeEstimationTT(DataEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(DataEstimation, self).__init__(
            name="fakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(
            Weight(
                "(0.5*ff1_nom*(byMediumDeepTau2017v2p1VSjet_1<0.5)+0.5*ff2_nom*(byMediumDeepTau2017v2p1VSjet_2<0.5))",
                "fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_1_iso")
        aiso_systematic.category.cuts.remove("tau_2_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "(byMediumDeepTau2017v2p1VSjet_2>0.5&&byMediumDeepTau2017v2p1VSjet_1<0.5&&byVVVLooseDeepTau2017v2p1VSjet_1>0.5)||(byMediumDeepTau2017v2p1VSjet_1>0.5&&byMediumDeepTau2017v2p1VSjet_2<0.5&&byVVVLooseDeepTau2017v2p1VSjet_2>0.5)",
                "tau_aiso"))
        return super(FakeEstimationTT,
                     self).create_root_objects(aiso_systematic)


class NMSSMEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal", heavy_mass = -1, light_mass= -1,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(NMSSMEstimation, self).__init__(
            name="NMSSM_{}_125_{}".format(heavy_mass, light_mass),
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            heavy_mass = heavy_mass,
            light_mass = light_mass,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("0.1*crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "NMSSMM{}h1M125tautauh2M{}$".format(self.heavy_mass, self.light_mass),
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        } 
        files = self.era.datasets_helper.get_nicks_with_query(query)
        files = [x for x in files if str(self.light_mass)+"_" in x]
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class HTTEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="HTT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(^GluGluHToTauTau.*125.*|^VBFHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,"generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHEstimation(HTTEstimation):
    htxs_dict = ggH_htxs

    def __init__(
            self,
            name,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(HTTEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        weights = super(ggHEstimation, self).get_weights()
        # weights.add(Weight("8.8384e-8/numberGeneratedEventsWeight", "ggh_stitching_weight")),
        weights.remove("numberGeneratedEventsWeight")
        weights.add(Weight("(numberGeneratedEventsWeight*(abs(crossSectionPerEventWeight - 3.0469376) > 1e-5)+1.0/(9673200 + 19939500 + 19977000)*(abs(crossSectionPerEventWeight - 3.0469376) < 1e-5))", "numberGeneratedEventsWeight"))  # 9673200 for inclusive sample and 19673200 for extension
        weights.add(Weight("ggh_NNLO_weight", "gghNNLO"))
        weights.add(Weight("1.01", "bbh_inclusion_weight"))
        return weights

    def get_cuts(self):
        return Cuts(
            Cut(self.htxs_dict[self.name],"htxs_match"))

    def get_files(self):
        query = {
            "process": "(^GluGluHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class qqHEstimation(HTTEstimation):
    htxs_dict = qqH_htxs

    def __init__(
            self,
            name,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(HTTEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut(self.htxs_dict[self.name],"htxs_match"))

    def get_weights(self):
        weights = super(qqHEstimation, self).get_weights()
        weights.remove("numberGeneratedEventsWeight")
        weights.add(Weight("(numberGeneratedEventsWeight*(abs(crossSectionPerEventWeight - 0.2370687)>1e-4)+1.0/(1499400 + 1999000 + 2997000)*(abs(crossSectionPerEventWeight - 0.2370687)<1e-4))", "numberGeneratedEventsWeight")) # 1499400 for inclusive sample and 1999000 for ext1 and 2997000 for ext2
        return weights

    def get_files(self):
        query = {
            "process":
            "(^VBFHToTauTau.*125.*|^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class VHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="VH",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=505)",
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="WH",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=305)",
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="ZH",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut("(htxs_stage1p1cat>=400)&&(htxs_stage1p1cat<=505)",
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^ZHToTauTau.*125.*|^ggZH.*ZToNuNu.*125.*|^ggZH.*ZToLL.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ttHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="ttH",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query = {
            "process": "(^ttHJetToTT.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class HWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,):
        super(HWWEstimation, self).__init__(
            name="HWW",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(VBF|GluGlu).*HToWWTo2L2Nu_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-JHUgenv628-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(ggHWWEstimation, self).__init__(
            name="ggHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "GluGluHToWWTo2L2Nu_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-JHUgenv628-pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class qqHWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(qqHWWEstimation, self).__init__(
            name="qqHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "VBFHToWWTo2L2Nu_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-JHUgenv628-pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WHWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(WHWWEstimation, self).__init__(
            name="WHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^HW(minus|plus)J_HToWW_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZHWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(ZHWWEstimation, self).__init__(
            name="ZHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(^HZJ_HToWW_M125|^GluGluZH_HToWW_M125)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class SUSYggHEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, mass, contribution, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,):
        super(SUSYggHEstimation, self).__init__(
            name="_".join(["gg"+contribution,str(mass)]),
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")
        self.mass = mass
        self.contribution = contribution

    def get_weights(self):
        contribution_weight = "1.0"
        if self.contribution in ["A_i", "A_t", "A_b", "H_i", "H_t", "H_b", "h_i", "h_t", "h_b"]:
            contribution_weight = "gg%s_weight"%self.contribution

        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight(contribution_weight, "contributionWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            #Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"), # not needed, since should be at 1.0

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^SUSYGluGluToHToTauTau_M{MASS}$".format(MASS=self.mass),
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class SUSYbbHEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, mass, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,):
        super(SUSYbbHEstimation, self).__init__(
            name="_".join(["bbH",str(mass)]),
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")
        self.mass = mass

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            #Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"), # not needed, since should be at 1.0

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^SUSYGluGluToBBHToTauTau_M{MASS}$".format(MASS=self.mass),
            "data": False,
            "campaign": self._mc_campaign,
            # "generator": "amcatnlo-pythia8", TODO: At the moment only LO samples available
            "generator": "^pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class bbH120Estimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="bbH120",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query = {
            "process": "(^SUSYGluGluToBBHToTauTau.*120$)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class bbH130Estimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="bbH130",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query = {
            "process": "(^SUSYGluGluToBBHToTauTau.*130$)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class DYJetsToLLEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        self.atNLO = atNLO
        name = "DYJetsToLLNLO" if self.atNLO else "DYJetsToLL"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        z_stitching_weight = Weight("(1.0)", "z_stitching_weight")
        if self.atNLO:
            z_stitching_weight = Weight(
                "((genbosonmass >= 50.0) * 5.0324e-05 + (genbosonmass < 50.0)*((abs(crossSectionPerEventWeight - 3.987) < 0.01)*4.6936e-06 + (abs(crossSectionPerEventWeight - 10.01) < 0.01)*3.7568e-06))",
                "z_stitching_weight"
            )  # xsec_NNLO [pb] = 2025.74*3, N_inclusive_NLO = 120762939, xsec_NNLO/N_inclusive_NLO = 5.0324e-05; fraction of negative events in 'generatorWeight'
        else:
            z_stitching_weight = Weight(
                "((genbosonmass >= 50.0) * 4.1545e-05*((npartons == 0 || npartons >= 5)*1.0+(npartons == 1)*0.32123574062076404+(npartons == 2)*0.3314444833963529+(npartons == 3)*0.3389929050626262+(npartons == 4)*0.2785338687268455) + (genbosonmass < 50.0)*(numberGeneratedEventsWeight * crossSectionPerEventWeight))",
                "z_stitching_weight")
        return Weights(
            Weight("generatorWeight", "generatorWeight"),
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
	        self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            z_stitching_weight,
            # lumi weight
            self.era.lumi_weight)

    def get_files(self):
        queryM10 = {
            "process": "DYJetsToLL_M10to50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v2"
        }
        queryM50_inclusive_2_3jet = {
            "process": "DY(|2|3)JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v2"
        }
        queryM50_1jet_v1 = {
            "process": "DY1JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "^$",
            "version": "v1"
        }
        queryM50_inc = {
            "process": "DYJetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "^$",
            "version": "v2"
        }
        queryM50_4jet = {
            "process": "DY4JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v2"
        }
        queryEWKZ = {
            "process": "^EWKZ",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "ext2"
        }
        queryM50NLO_inc = {
            "process": "DYJetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8",
        }
        if self.atNLO:
            files = self.era.datasets_helper.get_nicks_with_query(queryM50NLO_inc) + \
                    self.era.datasets_helper.get_nicks_with_query(queryEWKZ)
            log_query(self.name, queryM50NLO_inc, files)
        else:
            files = self.era.datasets_helper.get_nicks_with_query(queryM50_inc) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_inclusive_2_3jet) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_1jet_v1) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM10) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_4jet) + \
                    self.era.datasets_helper.get_nicks_with_query(queryEWKZ)
            log_query(self.name, queryM50_inc, files)
        return self.artus_file_names(files)


class EWKZEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        self.atNLO = atNLO
        super(DYJetsToLLEstimation, self).__init__(
            name="EWKZ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query_ewkz = {
            "process": "^EWKZ",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "ext2"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query_ewkz)
        log_query(self.name, query_ewkz, files)
        return self.artus_file_names(files)


class ZTTEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        self.atNLO = atNLO
        name = "ZTTNLO" if self.atNLO else "ZTT"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            ztt_cut = "gen_match_1==4 && gen_match_2==5"
        elif "et" in self.channel.name:
            ztt_cut = "gen_match_1==3 && gen_match_2==5"
        elif "tt" in self.channel.name:
            ztt_cut = "gen_match_1==5 && gen_match_2==5"
        elif "em" in self.channel.name:
            ztt_cut = "gen_match_1==3 && gen_match_2==4"
        elif "mm" in self.channel.name:
            ztt_cut = "gen_match_1==4 && gen_match_2==4"
        return Cuts(Cut(ztt_cut, "ztt_cut"))


class ZLEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        self.atNLO = atNLO
        name = "ZLNLO" if self.atNLO else "ZL"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    # def get_weights(self):
    #     # TODO remove this temp fix for the eletaufakerate as soon as the updated weights are in th ntuple
    #     weights = super(ZLEstimation, self).get_weights()
    #     weights.add(Weight("(gen_match_2==1 || gen_match_2==3)*(((abs(eta_1) < 1.46) * (1./0.6) * 1.22) + ((abs(eta_1) > 1.5588) * (1./0.88) * 1.47))+!(gen_match_2==1 || gen_match_2==3)", "eletauFakeRateWeightFix"))
    #     return weights

    def get_cuts(self):
        if "mt" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "et" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            emb_veto = "!(gen_match_1==5 && gen_match_2==5)"
            ff_veto = "!(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==4)"
            ff_veto = "(1.0)"
        elif "mm" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==4)"
            ff_veto = "(1.0)"
        return Cuts(Cut("%s && %s" % (emb_veto, ff_veto),
                        "dy_emb_and_ff_veto"))


class ZJEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        self.atNLO = atNLO
        name = "ZJNLO" if self.atNLO else "ZJ"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "gen_match_2 == 6"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0 == 1"
        return Cuts(Cut(ct, "dy_fakes"))


class ZTTEmbeddedEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(ZTTEmbeddedEstimation, self).__init__(
            name="EMB",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            friend_directory=friend_directory,
            directory=directory,
            channel=channel,
            mc_campaign=None)

    def emb_triggerweights(self):
        channel = self.channel.name
        weight = Weight("1.0", "triggerweight")

        singleEMB = "singleTriggerEmbeddedEfficiencyWeightKIT_1"
        crossEMBL = "crossTriggerEmbeddedEfficiencyWeightKIT_1"
        EMBTau_1 = "((byMediumDeepTau2017v2p1VSjet_1<0.5 && byVVVLooseDeepTau2017v2p1VSjet_1>0.5)*crossTriggerEMBEfficiencyWeight_vloose_DeepTau_1 + (byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerEMBEfficiencyWeight_medium_DeepTau_1)" # hotfix
        # EMBTau_1 = "((byMediumDeepTau2017v2p1VSjet_1<0.5 && byVVVLooseDeepTau2017v2p1VSjet_1>0.5)*crossTriggerEMBEfficiencyWeight_vloose_DeepTau_1 + (byMediumDeepTau2017v2p1VSjet_1>0.5)*crossTriggerEMBEfficiencyWeight_medium_DeepTau_1)"
        EMBTau_2 = EMBTau_1.replace("_1", "_2")

        if "mt" in channel:
            trig_sL = "(pt_1 >= 23 && trg_singlemuon)"
            trig_X = "(pt_1 < 23 && trg_mutaucross)"

            MuTauEMB = "{singletrigger} + {crosstrigger}".format(
                singletrigger="*".join([trig_sL, singleEMB]),
                crosstrigger="*".join([trig_X, crossEMBL, EMBTau_2]))
            MuTauData = MuTauEMB.replace("EMB", "Data").replace("Embedded", "Data")
            MuTau = "(" + MuTauData + ")/(" + MuTauEMB + ")"
            weight = Weight(MuTau, "triggerweight")

        elif "et" in channel:
            trig_sL = "(trg_singleelectron)"
            trig_X = "(pt_1 > 25 && pt_1 < 26 && trg_eletaucross)"

            ElTauEMB = "{singletrigger} + {crosstrigger}".format(
                singletrigger="*".join([trig_sL, singleEMB]),
                crosstrigger="*".join([trig_X, crossEMBL, EMBTau_2])
                )
            ElTauData = ElTauEMB.replace("EMB", "Data").replace("Embedded", "Data")
            ElTau = "(" + ElTauData + ")/(" + ElTauEMB + ")"
            weight = Weight(ElTau, "triggerweight")

        elif "tt" in channel:
            DiTauEMB = "*".join([EMBTau_1, EMBTau_2])
            DiTauData = DiTauEMB.replace("EMB", "Data").replace("Embedded", "Data")
            DiTau = "("+DiTauData+")/("+DiTauEMB+")"
            weight = Weight(DiTau, "triggerweight")

        elif "em" in channel:
            ElMuData = "(trigger_23_data_Weight_2*trigger_12_data_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_data_Weight_1*trigger_8_data_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_data_Weight_2*trigger_23_data_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))"
            ElMuEmb = ElMuData.replace('data', 'embed')
            ElMu = "("+ElMuData+")/("+ElMuEmb+")"
            weight = Weight(ElMu, "triggerweight")
        return weight

    def get_weights(self):
        emb_weights = Weights(
            Weight("generatorWeight*(generatorWeight<=1.0)", "simulation_sf"),
            Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2",
                   "scale_factor"),
            Weight("embeddedDecayModeWeight", "decayMode_SF"))
        if self.channel.name == "mt":
            emb_weights.add(Weight("idWeight_1*isoWeight_1", "lepton_sf"))
            emb_weights.add(self.get_tauByIsoIdWeight_for_channel(self.channel))
            emb_weights.add(
                Weight("gen_match_1==4 && gen_match_2==5", "emb_veto"))
            emb_weights.add(self.emb_triggerweights())

        elif self.channel.name == "et":
            emb_weights.add(Weight("idWeight_1*isoWeight_1", "lepton_sf"))
            emb_weights.add(self.get_tauByIsoIdWeight_for_channel(self.channel)),
            emb_weights.add(
                Weight("gen_match_1==3 && gen_match_2==5", "emb_veto"))
            emb_weights.add(self.emb_triggerweights())

        elif self.channel.name == "tt":
            emb_weights.add(self.emb_triggerweights())
            emb_weights.add(self.get_tauByIsoIdWeight_for_channel(self.channel))
            emb_weights.add(
                Weight("gen_match_1==5 && gen_match_2==5", "emb_veto"))
        elif self.channel.name == "em":
            emb_weights.add(
                Weight("(gen_match_1==3 && gen_match_2==4)", "emb_veto")
            )
            emb_weights.add(self.emb_triggerweights())
            emb_weights.add(
                Weight("idWeight_1*isoWeight_1*idWeight_2*looseIsoWeight_2",
                       "leptopn_sf"))
            emb_weights.remove(
                "decayMode_SF"
            )  # embeddedDecayModeWeight is only for tau decay modes
        return emb_weights

    def get_files(self):
        query = {"process": "Embedding2016(B|C|D|E|F|G|H)", "embedded": True}
        if self.channel.name == "mt":
            query["campaign"] = "MuTauFinalState"
            query["scenario"] = "inputDoubleMu_94X_Legacy_miniAODv"
        elif self.channel.name == "et":
            query["campaign"] = "ElTauFinalState"
            query["scenario"] = "inputDoubleMu_94X_Legacy_miniAODv"
        elif self.channel.name == "tt":
            query["campaign"] = "TauTauFinalState"
            query["scenario"] = "inputDoubleMu_94X_Legacy_miniAODv"
        elif self.channel.name == "em":
            query["campaign"] = "ElMuFinalState"
            query["scenario"] = "inputDoubleMu_94X_Legacy_miniAODv"
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        return Cuts(
            Cut(
                "((gen_match_1>2 && gen_match_1<6) &&  (gen_match_2>2 && gen_match_2<6))",
                "dy_genuine_tau"))


class WEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
    ):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC related weights
            #Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            # # xsec_NNLO [pb] = 61526.7, N_inclusive = 86916455, xsec_NNLO/N_inclusive = 0.00070788321 [pb] weights: [1.0, 0.2691615837248596, 0.1532341436287767, 0.03960756033932645, 0.03969970742404736]
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "((0.00070788321*((npartons <= 0 || npartons >= 5)*1.0 + (npartons == 1)*0.2691615837248596 + (npartons == 2)*0.1532341436287767 + (npartons == 3)*0.03960756033932645 + (npartons == 4)*0.03969970742404736)) * (genbosonmass>=0.0) + numberGeneratedEventsWeight * crossSectionPerEventWeight * (genbosonmass<0.0))",
                "wj_stitching_weight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.?JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "^EWKW",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "ext2"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKWpEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(EWKWpEstimation, self).__init__(
            name="EWKWp",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "(5.190747826298e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKWp_stitching_weight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKWPlus",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKWmEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(EWKWmEstimation, self).__init__(
            name="EWKWm",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "(4.200367267668e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKW_stitching_weight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKWMinus",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WEstimationRaw(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(WEstimationRaw, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "(((npartons == 0 || npartons >= 5)*7.09390278348407e-4) + ((npartons == 1)*1.90063898596475e-4) + ((npartons == 2)*5.8529964471165e-5) + ((npartons == 3)*1.9206444928444e-5) + ((npartons == 4)*1.923548021385e-5))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)",
                "wj_stitching_weight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.*JetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


# class WEstimation(SumUpEstimationMethod):
#     def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
#         super(WEstimation, self).__init__(
#             name="W",
#             folder=folder,
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
#             era=era,
#             directory=directory,
#             friend_directory=friend_directory,
#             channel=channel,
#             processes=[
#                 Process(
#                     "W",
#                     WEstimationRaw(
#                         era,
#                         directory,
#                         channel,
#                         friend_directory=friend_directory)),
#                 Process(
#                     "EWKWp",
#                     EWKWpEstimation(
#                         era,
#                         directory,
#                         channel,
#                         friend_directory=friend_directory)),
#                 Process(
#                     "EWKWm",
#                     EWKWmEstimation(
#                         era,
#                         directory,
#                         channel,
#                         friend_directory=friend_directory))
#             ])


class WTEstimation(WEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(Cut("gen_match_1==3||gen_match_1==4", "wt_genmatch"))


class WLEstimation(WEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(Cut("!(gen_match_1==3||gen_match_1==4)", "wl_genmatch"))


class TTEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("topPtReweightWeightTTH", "topPtReweightWeight"),
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            #self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight","leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)


    def get_files(self):
        query = {
            "process": "TTTo(2L2Nu|Hadronic|SemiLeptonic).*",
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTLEstimation(TTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "et" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            emb_veto = "!(gen_match_1==5 && gen_match_2==5)"
            ff_veto = "!(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==4)"
            ff_veto = "(1.0)"
        elif "mm" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==4)"
            ff_veto = "(1.0)"
        return Cuts(Cut("%s && %s" % (emb_veto, ff_veto),
                        "tt_emb_and_ff_veto"))


class TTTEstimation(TTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==5"
        elif "et" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==5"
        elif "tt" in self.channel.name:
            tt_cut = "gen_match_1==5 && gen_match_2==5"
        elif "em" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==4"
        elif "mm" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==4"
        return Cuts(Cut(tt_cut, "ttt_cut"))


class TTJEstimation(TTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TTJ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "(gen_match_2 == 6 && gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0.0 == 1.0"
        return Cuts(Cut(ct, "tt_fakes"))


class VVEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VV",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2","isoWeight"),
            Weight("idWeight_1*idWeight_2","idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            aiso_muon_correction(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("1.252790591041545e-07*(abs(crossSectionPerEventWeight - 118.7) < 0.01) + 5.029933132068942e-07*(abs(crossSectionPerEventWeight - 12.14) < 0.01) + 2.501519047441559e-07*(abs(crossSectionPerEventWeight - 22.82) < 0.01) + numberGeneratedEventsWeight*(abs(crossSectionPerEventWeight - 118.7) > 0.01 && abs(crossSectionPerEventWeight - 12.14) > 0.01 && abs(crossSectionPerEventWeight - 22.82) > 0.01)","numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(WZTo3LNu|WZTo2L2Q|ZZTo2L2Q|ZZTo4L)$",  # Query for Di-Boson samples
            "data": False,
            "generator": "amcatnlo-pythia8",
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "(VVTo2L2Nu)$",  # Query for Di-Boson samples
            "data": False,
            "extension": "ext1",
            "campaign": self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process":
            "(STt-channelantitop4finclusiveDecays|STt-channeltop4finclusiveDecays|STtWantitop5finclusiveDecays|STtWtop5finclusiveDecays)",
            "data":
            False,
            "campaign":
            self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, "<optimzed out>", files)
        return self.artus_file_names(files)


# class VVEstimation(
#         EstimationMethod
# ):  # old method of estimation not using inclusive diboson samples
#     def __init__(
#             self,
#             era,
#             directory,
#             channel,
#             friend_directory=None,
#             folder="nominal",
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
#         super(VVEstimation, self).__init__(
#             name="VV",
#             folder=folder,
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=
#             get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
#             era=era,
#             directory=directory,
#             friend_directory=friend_directory,
#             channel=channel,
#             mc_campaign="RunIISummer16MiniAODv3")

#     def get_weights(self):
#         # if self.channel.name=="em":
#         #     return Weights(
#         #         Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
#         #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
#         #         Weight("triggerWeight_1*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
#         #         Weight("puweight", "puweight"),
#         #         self.get_tauByIsoIdWeight_for_channel(self.channel),
#         #         Weight("prefiringweight", "prefireWeight"),
#         #         self.era.lumi_weight)
#         # else:
#         return Weights(
#             Weight("isoWeight_1*isoWeight_2", "isoWeight"),
#             Weight("idWeight_1*idWeight_2", "idWeight"),
#             self.get_tauByIsoIdWeight_for_channel(self.channel),
#             Weight("puweight", "puweight"),
#             Weight("trackWeight_1*trackWeight_2", "trackweight"),
#             self.get_triggerweight_for_channel(self.channel._name),
#             aiso_muon_correction(self.channel._name),
#             self.get_singlelepton_triggerweight_for_channel(self.channel.name),
#             Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
#                    "leptonTauFakeRateWeight"),
#             get_eleRecoWeight_for_channel(self.channel.name),
#             Weight("prefiringweight", "prefireWeight"),
#             # MC weights
#             Weight("numberGeneratedEventsWeight",
#                    "numberGeneratedEventsWeight"),
#             Weight(
#                 "(1.0+0.56*(abs(crossSectionPerEventWeight-75.769996)<0.00001))",
#                 "VV_NNLO_reweight"),
#             Weight("generatorWeight", "generatorWeight"),
#             self.era.lumi_weight)

#     def get_files(self):
#         query = {
#             "process":
#             "(WWTo1L1Nu2Q|" + "WZTo1L1Nu2Q|" + "WZTo1L3Nu|" + "WZTo2L2Q|" +
#             "ZZTo2L2Q" + ")",
#             "data":
#             False,
#             "campaign":
#             self._mc_campaign,
#             "generator":
#             "amcatnlo-pythia8"
#         }
#         files = self.era.datasets_helper.get_nicks_with_query(query)

#         query = {
#             "process": "(VVTo2L2Nu|ZZTo4L)",
#             "extension": "ext1",
#             "data": False,
#             "campaign": self._mc_campaign,
#             "generator": "amcatnlo-pythia8"
#         }
#         files += self.era.datasets_helper.get_nicks_with_query(query)

#         query = {
#             "process": "WZJToLLLNu",
#             "data": False,
#             "campaign": self._mc_campaign,
#             "generator": "pythia8"
#         }
#         files += self.era.datasets_helper.get_nicks_with_query(query)

#         query = {
#             "process":
#             "(STt-channelantitop4finclusiveDecays|STt-channeltop4finclusiveDecays|STtWantitop5finclusiveDecays|STtWtop5finclusiveDecays)",
#             "data": False,
#             "campaign": self._mc_campaign
#         }
#         files += self.era.datasets_helper.get_nicks_with_query(query)

#         log_query(self.name, "<optimzed out>", files)
#         return self.artus_file_names(files)


class VVLEstimation(VVEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VVL",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "et" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            emb_veto = "!(gen_match_1==5 && gen_match_2==5)"
            ff_veto = "!(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==4)"
            ff_veto = "(1.0)"
        elif "mm" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==4)"
            ff_veto = "(1.0)"
        return Cuts(Cut("%s && %s" % (emb_veto, ff_veto),
                        "vv_emb_and_ff_veto"))


class VVTEstimation(VVEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VVT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==5"
        elif "et" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==5"
        elif "tt" in self.channel.name:
            tt_cut = "gen_match_1==5 && gen_match_2==5"
        elif "em" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==4"
        elif "mm" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==4"
        return Cuts(Cut(tt_cut, "vvt_cut"))


class VVJEstimation(VVEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VVJ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "(gen_match_2 == 6 && gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0.0 == 1.0"
        return Cuts(Cut(ct, "vv_fakes"))


class QCDEstimation_SStoOS_MTETEM(SStoOSEstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            friend_directory=None,
            extrapolation_factor=1.0,
            folder="nominal",
            qcd_weight=Weight("1.0", "qcd_Weight"),
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(QCDEstimation_SStoOS_MTETEM, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            friend_directory=friend_directory,
            data_process=data_process,
            extrapolation_factor=extrapolation_factor,
            qcd_weight=qcd_weight)


class QCDEstimationTT(ABCDEstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(QCDEstimationTT, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=
            [  # cuts to be removed to include region for shape derivation
                "tau_2_iso"
            ],
            BD_cuts=
            [  # cuts to be applied to restrict to region for shape derivation
                Cut("byMediumDeepTau2017v2p1VSjet_2<0.5",
                    "tau_2_iso"),
                Cut("byLooseDeepTau2017v2p1VSjet_2>0.5",
                    "tau_2_iso_loose")
            ],
            AB_cut_names=
            [  # cuts to be removed to include region for the determination of the extrapolation derivation
                "os"
            ],
            CD_cuts=
            [  # cuts to be applied to restrict to region for the determination of the extrapolation derivation
                Cut("q_1*q_2>0", "ss")
            ])


class WEstimationWithQCD(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            w_process,
            qcd_ss_to_os_extrapolation_factor,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(WEstimationWithQCD, self).__init__(
            name="WJets",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for MC WJets shape estimation in the signal region
        signal_region = copy.deepcopy(systematic.category)
        signal_region.name = (signal_region.name +
                              "_for_wjets_mc").lstrip(self.channel.name + "_")

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("m_t")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70", "m_t"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(
            systematic.category)  # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("m_t")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70", "m_t"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("m_t")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("m_t")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [
                high_mt_os_control_region, low_mt_os_control_region,
                high_mt_ss_control_region, low_mt_ss_control_region,
                inclusive_os_control_region, inclusive_ss_control_region
        ]:
            s = Systematic(category=category,
                           process=self._w_process,
                           analysis=systematic.analysis,
                           era=self.era,
                           variation=systematic.variation,
                           mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [
                    high_mt_os_control_region, high_mt_ss_control_region
            ]:
                s = Systematic(category=category,
                               process=process,
                               analysis=systematic.analysis,
                               era=self.era,
                               variation=systematic.variation,
                               mass=125)
                if process == self._data_process:
                    direction = s.variation._direction
                    s.variation = Nominal()
                    s.variation._direction = direction
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for signal region shape
        s = Systematic(category=signal_region,
                       process=self._w_process,
                       analysis=systematic.analysis,
                       era=self.era,
                       variation=systematic.variation,
                       mass=125)
        systematic._WandQCD_systematics.append(s)
        s.create_root_objects()
        root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        wjets_mc_shape = None
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("for_wjets_mc"):
                wjets_mc_shape = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result / wjets_ss_cr_count.result

        wjets_integral_low_mt_os = wjets_low_mt_os_cr_count.result
        wjets_integral_high_mt_os = wjets_high_mt_os_cr_counts.pop(
            self._w_process.name).result
        logger.debug("Integral of WJets MC in low mt OS region: %s",
                     str(wjets_integral_low_mt_os))
        logger.debug("Integral of WJets MC in high mt OS region: %s",
                     str(wjets_integral_high_mt_os))

        R_high_to_low_mt_os = wjets_integral_low_mt_os / wjets_integral_high_mt_os
        R_high_to_low_mt_ss = wjets_low_mt_ss_cr_count.result / wjets_high_mt_ss_cr_counts.pop(
            self._w_process.name).result
        logger.debug("WJets SS to OS extrapolation factor: %s",
                     str(R_ss_to_os))
        logger.debug("WJets high to low mt os extrapolation factor: %s",
                     str(R_high_to_low_mt_os))
        logger.debug("WJets high to low mt ss extrapolation factor: %s",
                     str(R_high_to_low_mt_ss))

        # Determine yields in wjets CRs
        logger.debug(
            "Data yield in ss high mt region: %s",
            str(wjets_high_mt_ss_cr_counts[self._data_process.name].result))
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_ss_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_ss_cr_counts.values()])
        logger.debug("MC yield to be subtracted: %s", str(sum_mc))
        for name, s in wjets_high_mt_ss_cr_counts.items():
            logger.debug(name + " : " + str(s.result / sum_mc))

        logger.debug(
            "Data yield in os high mt region: %s",
            str(wjets_high_mt_os_cr_counts[self._data_process.name].result))
        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_os_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_os_cr_counts.values()])
        logger.debug("MC yield to be subtracted: %s", str(sum_mc))
        for name, s in wjets_high_mt_os_cr_counts.items():
            logger.debug(name + " : " + str(s.result / sum_mc))

        logger.debug("WJets + QCD yield in ss high mt region: %s",
                     str(high_mt_ss_yield))
        logger.debug("WJets + QCD yield in os high mt region: %s",
                     str(high_mt_os_yield))

        # Derive and normalize final shape
        logger.debug("WJets MC yield in signal region: %s",
                     str(wjets_integral_low_mt_os))
        sf = R_ss_to_os * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor
            ) / wjets_integral_high_mt_os
        estimated_yield = R_high_to_low_mt_os * R_ss_to_os * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor)
        logger.debug("WJets Estimated yield in signal region: %s",
                     str(estimated_yield))
        logger.debug("Scale WJets by %s", str(sf))
        wjets_shape = copy.deepcopy(wjets_mc_shape)
        wjets_shape.result.Scale(sf)

        # Rename root object accordingly
        wjets_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        wjets_shape.replace_negative_entries_and_renormalize(tolerance=100.5)

        return wjets_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class QCDEstimationWithW(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            w_process,
            qcd_ss_to_os_extrapolation_factor,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel):
        super(QCDEstimationWithW, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for WJets and QCD shape estimation in the qcd control region
        qcd_control_region = copy.deepcopy(systematic.category)
        qcd_control_region.name = (qcd_control_region.name +
                                   "_ss_for_qcd").lstrip(self.channel.name +
                                                         "_")

        qcd_control_region.cuts.remove("os")

        qcd_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("m_t")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70", "m_t"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(
            systematic.category)  # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("m_t")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70", "m_t"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("m_t")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("m_t")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [
                high_mt_os_control_region, low_mt_os_control_region,
                high_mt_ss_control_region, low_mt_ss_control_region,
                inclusive_os_control_region, inclusive_ss_control_region
        ]:
            s = Systematic(category=category,
                           process=self._w_process,
                           analysis=systematic.analysis,
                           era=self.era,
                           variation=systematic.variation,
                           mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [
                    high_mt_os_control_region, high_mt_ss_control_region
            ]:
                s = Systematic(category=category,
                               process=process,
                               analysis=systematic.analysis,
                               era=self.era,
                               variation=systematic.variation,
                               mass=125)
                if process == self._data_process:
                    direction = s.variation._direction
                    s.variation = Nominal()
                    s.variation._direction = direction
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for Wjets and QCD shape
        for process in [self._data_process, self._w_process
                        ] + self._bg_processes:
            s = Systematic(category=qcd_control_region,
                           process=process,
                           analysis=systematic.analysis,
                           era=self.era,
                           variation=systematic.variation,
                           mass=125)
            if process == self._data_process:
                direction = s.variation._direction
                s.variation = Nominal()
                s.variation._direction = direction
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()

            root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        qcd_control_region_shapes = {}
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("ss_for_qcd"):
                qcd_control_region_shapes[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result / wjets_ss_cr_count.result

        wjets_integral_low_mt_ss = wjets_low_mt_ss_cr_count.result
        wjets_integral_high_mt_ss = wjets_high_mt_ss_cr_counts.pop(
            self._w_process.name).result

        R_high_to_low_mt_os = wjets_low_mt_os_cr_count.result / wjets_high_mt_os_cr_counts.pop(
            self._w_process.name).result
        R_high_to_low_mt_ss = wjets_integral_low_mt_ss / wjets_integral_high_mt_ss

        # Determine yields in wjets CRs
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_ss_cr_counts.values()])

        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_os_cr_counts.values()])

        # Derive and normalize final shape for QCD
        wjets_shape = qcd_control_region_shapes.pop(self._w_process.name)
        logger.debug("WJets MC yield in qcd control region: %s",
                     str(wjets_integral_low_mt_ss))
        sf = (high_mt_os_yield -
              self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                  R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor
              ) / wjets_integral_high_mt_ss
        estimated_yield = R_high_to_low_mt_ss * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor)
        logger.debug("WJets Estimated yield in qcd control region: %s",
                     str(estimated_yield))
        logger.debug("Scale WJets by %s", str(sf))
        wjets_shape.result.Scale(sf)
        wjets_shape._result.Write()

        qcd_shape = copy.deepcopy(
            qcd_control_region_shapes.pop(self._data_process.name))
        qcd_shape.result.Add(wjets_shape.result, -1.0)
        for sh in qcd_control_region_shapes.values():
            qcd_shape.result.Add(sh.result, -1.0)
        # Saving QCD shape in ss control region
        qcd_ss_shape = copy.deepcopy(qcd_shape)
        ss_category_name = ""
        for s in systematic._WandQCD_systematics:
            if s.category.name.endswith("ss_for_qcd"):
                ss_category_name = s.category._name
        qcd_ss_shape.name = systematic.name.replace(systematic.category._name,
                                                    ss_category_name)
        qcd_ss_shape._result.Write()

        # Rescale QCD shape for signal region
        qcd_shape.result.Scale(self._qcd_ss_to_os_extrapolation_factor)

        # Rename root object accordingly
        qcd_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        qcd_shape.replace_negative_entries_and_renormalize(tolerance=100.5)

        return qcd_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError
