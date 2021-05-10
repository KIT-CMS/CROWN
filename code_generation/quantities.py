class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = {}
        self.children = {}

    def get_leaf(self, shift, scope):
        if shift in self.get_shifts(scope):
            return self.name + shift
        return self.name

    def get_leafs_of_scope(self, scope):
        return [self.name] + [self.name + shift for shift in self.get_shifts(scope)]

    def shift(self, name, scope):
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

ditaupair = Quantity("ditaupair")

p4_1 = Quantity("p4_1")
pt_1 = Quantity("pt_1")
eta_1 = Quantity("eta_1")
phi_1 = Quantity("phi_1")
p4_2 = Quantity("p4_2")
pt_2 = Quantity("pt_2")
eta_2 = Quantity("eta_2")
phi_2 = Quantity("phi_2")
