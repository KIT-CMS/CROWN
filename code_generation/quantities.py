class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = []
        self.children = []

    def get_leaf(self, shift):
        # add sth here that creates name of shifted quantity
        if shift in self.shifts:
            return self.name + shift
        return self.name

    def shift(self, name):
        if not name in self.shifts:
            self.shifts.append(name)
            for c in self.children:
                c.shift(name)

    def copy(self, name):
        copy = Quantity(name)
        copy.shifts = self.shifts
        copy.children = self.children
        return copy


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

ditaupair = Quantity("ditaupair")

p4_1 = Quantity("p4_1")
pt_1 = Quantity("pt_1")
eta_1 = Quantity("eta_1")
phi_1 = Quantity("phi_1")
p4_2 = Quantity("p4_2")
pt_2 = Quantity("pt_2")
eta_2 = Quantity("eta_2")
phi_2 = Quantity("phi_2")
