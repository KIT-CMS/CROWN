import logging

log = logging.getLogger(__name__)


class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = {}
        self.children = {}
        log.debug("Setting up new Quantity {}".format(self.name))

    def get_leaf(self, shift, scope):
        if shift in self.get_shifts(scope):
            return self.name + shift
        return self.name

    def get_leafs_of_scope(self, scope):
        return [self.name] + [self.name + shift for shift in self.get_shifts(scope)]

    def shift(self, name, scope):
        log.debug("Adding shift {} to quantity {}".format(name, self.name))
        if not scope in self.shifts.keys():
            self.shifts[scope] = set()
        if not name in self.shifts[scope]:
            self.shifts[scope].add(name)
            if scope == "global":  # shift children in all scopes if scope is global
                for any_scope in self.children.values():
                    for c in any_scope:
                        c.shift(name, scope)
            else:
                if scope in self.children.keys():
                    for c in self.children[scope]:
                        c.shift(name, scope)

    def copy(self, name):
        copy = Quantity(name)
        copy.shifts = self.shifts
        copy.children = self.children
        return copy

    def adopt(self, child, scope):
        log.debug(
            "Adopting child quantity {} to quantity {}".format(child.name, self.name)
        )
        if not scope in self.children.keys():
            self.children[scope] = []
        self.children[scope].append(child)

    def get_shifts(self, scope):
        if scope in self.shifts.keys():
            if "global" in self.shifts.keys():
                return list(self.shifts[scope].union(self.shifts["global"]))
            else:
                return list(self.shifts[scope])
        else:
            if "global" in self.shifts.keys():
                return list(self.shifts["global"])
            else:
                return []


m_vis = Quantity("m_vis")
mt_1 = Quantity("mt_1")

good_taus_mask = Quantity("good_taus_mask")
Tau_pt = Quantity("Tau_pt")
Tau_eta = Quantity("Tau_eta")
Tau_phi = Quantity("Tau_phi")
Tau_mass = Quantity("Tau_mass")
Tau_dz = Quantity("Tau_dz")
Tau_IDraw = Quantity("Tau_rawDeepTau2017v2p1VSjet")

good_muons_mask = Quantity("good_muons_mask")
Muon_pt = Quantity("Muon_pt")
Muon_eta = Quantity("Muon_eta")
Muon_phi = Quantity("Muon_phi")
Muon_mass = Quantity("Muon_mass")
Muon_iso = Quantity("Muon_pfRelIso04_all")

electron_veto_mask = Quantity("electron_veto_mask")
electron_veto_flag = Quantity("extraelec_veto")
Electron_pt = Quantity("Electron_pt")
Electron_eta = Quantity("Electron_eta")
Electron_phi = Quantity("Electron_phi")
Electron_mass = Quantity("Electron_mass")
Electron_iso = Quantity("Electron_pfRelIso03_all")

jet_id_mask = Quantity("jet_id_mask")
jet_overlap_veto_mask = Quantity("jet_overlap_veto_mask")
good_jets_mask = Quantity("good_jets_mask")
good_bjets_mask = Quantity("good_bjets_mask")
Jet_eta = Quantity("Jet_eta")
Jet_phi = Quantity("Jet_phi")
Jet_pt = Quantity("Jet_pt")
Jet_mass = Quantity("Jet_mass")
ditaupair = Quantity("ditaupair")
good_jet_collection = Quantity("good_jet_collection")
good_bjet_collection = Quantity("good_bjet_collection")

p4_1 = Quantity("p4_1")
pt_1 = Quantity("pt_1")
eta_1 = Quantity("eta_1")
phi_1 = Quantity("phi_1")
p4_2 = Quantity("p4_2")
pt_2 = Quantity("pt_2")
eta_2 = Quantity("eta_2")
phi_2 = Quantity("phi_2")

njets = Quantity("njets")
nbtag = Quantity("nbtag")
jet_p4_1 = Quantity("jet_p4_1")
jpt_1 = Quantity("jpt_1")
jeta_1 = Quantity("jeta_1")
jphi_1 = Quantity("jphi_1")
jet_p4_2 = Quantity("jet_p4_2")
jpt_2 = Quantity("jpt_2")
jeta_2 = Quantity("jeta_2")
jphi_2 = Quantity("jphi_2")
mjj = Quantity("mjj")
bjet_p4_1 = Quantity("bjet_p4_1")
bpt_1 = Quantity("bpt_1")
beta_1 = Quantity("beta_1")
bphi_1 = Quantity("bphi_1")
bjet_p4_2 = Quantity("bjet_p4_2")
bpt_2 = Quantity("bpt_2")
beta_2 = Quantity("beta_2")
bphi_2 = Quantity("bphi_2")
