import logging

log = logging.getLogger(__name__)


class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = {}
        self.children = {}
        self.output_scopes = []
        log.debug("Setting up new Quantity {}".format(self.name))

    # the scopes, in which a quantity is used as an output is tracked in the output_scopes list.
    # This check is triggered for every producer.
    # If a quantity is already used within a given scope as output, this will result in an exception.
    def check_scope(self, scope):
        log.debug("Checking {} / scope {}".format(self.name, scope))
        if scope not in self.output_scopes:
            self.output_scopes.append(scope)
        else:
            log.error(
                "Quantity {} is already used as output in {} scope !".format(
                    self.name, scope
                )
            )
            raise Exception

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


class NanoAODQuantity(Quantity):
    def __init__(self, name):
        self.name = name
        super().__init__(name)

    # Quantities from the NanoAOD are not designed to be directly usable as output
    def check_scope(self, scope):
        log.error(
            "Quantity {} is a NanoAOD quantity and cant be used as output !".format(
                self.name
            )
        )
        raise Exception


m_vis = Quantity("m_vis")
mt_1 = Quantity("mt_1")

good_taus_mask = Quantity("good_taus_mask")
good_muons_mask = Quantity("good_muons_mask")
electron_veto_mask = Quantity("electron_veto_mask")
electron_veto_flag = Quantity("extraelec_veto")
jet_id_mask = Quantity("jet_id_mask")
jet_overlap_veto_mask = Quantity("jet_overlap_veto_mask")
good_jets_mask = Quantity("good_jets_mask")
good_bjets_mask = Quantity("good_bjets_mask")
Jet_pt_corrected = Quantity("Jet_pt_corrected")
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
mass_1 = Quantity("mass_1")
mass_2 = Quantity("mass_2")
dxy_1 = Quantity("dxy_1")
dxy_2 = Quantity("dxy_2")
dz_1 = Quantity("dz_1")
dz_2 = Quantity("dz_2")
q_1 = Quantity("q_1")
q_2 = Quantity("q_2")
iso_1 = Quantity("iso_1")
iso_2 = Quantity("iso_2")
decaymode_1 = Quantity("decaymode_1")
decaymode_2 = Quantity("decaymode_2")
gen_match_1 = Quantity("gen_match_1")
gen_match_2 = Quantity("gen_match_2")

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
