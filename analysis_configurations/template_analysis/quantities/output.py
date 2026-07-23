from code_generation.quantity import Quantity

lumi = Quantity()

base_muons_mask = Quantity()
good_muons_mask = Quantity()
nmuons = Quantity()

## MM pair quantities
dileptonpair = Quantity()

p4_1 = Quantity()
pt_1 = Quantity()
eta_1 = Quantity()
phi_1 = Quantity()
mass_1 = Quantity()
q_1 = Quantity()

p4_2 = Quantity()
pt_2 = Quantity()
eta_2 = Quantity()
phi_2 = Quantity()
mass_2 = Quantity()
q_2 = Quantity()

p4_mm_pair = Quantity()
mm_pair_mass = Quantity()
mm_pair_pt = Quantity()


## Gen Quantities
gen_dileptonpair = Quantity()

gen_p4_1 = Quantity()
gen_pt_1 = Quantity()
gen_eta_1 = Quantity()
gen_phi_1 = Quantity()
gen_mass_1 = Quantity()
gen_pdgid_1 = Quantity()

gen_p4_2 = Quantity()
gen_pt_2 = Quantity()
gen_eta_2 = Quantity()
gen_phi_2 = Quantity()
gen_mass_2 = Quantity()
gen_pdgid_2 = Quantity()

gen_p4_mm_pair = Quantity()
gen_mm_pair_mass = Quantity()

# sample flags
is_data = Quantity()
is_embedding = Quantity()
is_ttbar = Quantity()
is_dyjets = Quantity()
is_wjets = Quantity()
is_diboson = Quantity()

# Muon weights
reco_wgt_mu_1 = Quantity()
reco_wgt_mu_2 = Quantity()
id_wgt_mu_1 = Quantity()
id_wgt_mu_2 = Quantity()
iso_wgt_mu_1 = Quantity()
iso_wgt_mu_2 = Quantity()
