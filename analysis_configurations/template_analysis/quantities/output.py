from code_generation.quantity import Quantity

lumi = Quantity("lumi")

base_muons_mask = Quantity("base_muons_mask")
good_muons_mask = Quantity("good_muons_mask")
nmuons = Quantity("nmuons")

## MM pair quantities
dileptonpair = Quantity("dileptonpair")

p4_1 = Quantity("p4_1")
pt_1 = Quantity("pt_1")
eta_1 = Quantity("eta_1")
phi_1 = Quantity("phi_1")
mass_1 = Quantity("mass_1")
q_1 = Quantity("q_1")

p4_2 = Quantity("p4_2")
pt_2 = Quantity("pt_2")
eta_2 = Quantity("eta_2")
phi_2 = Quantity("phi_2")
mass_2 = Quantity("mass_2")
q_2 = Quantity("q_2")

mm_pair_mass = Quantity("mm_pair_mass")
mm_pair_pt = Quantity("mm_pair_pt")

## Gen Quantities
gen_dileptonpair = Quantity("gen_dileptonpair")

gen_p4_1 = Quantity("gen_p4_1")
gen_pt_1 = Quantity("gen_pt_1")
gen_eta_1 = Quantity("gen_eta_1")
gen_phi_1 = Quantity("gen_phi_1")
gen_mass_1 = Quantity("gen_mass_1")
gen_pdgid_1 = Quantity("gen_pdgid_1")

gen_p4_2 = Quantity("gen_p4_2")
gen_pt_2 = Quantity("gen_pt_2")
gen_eta_2 = Quantity("gen_eta_2")
gen_phi_2 = Quantity("gen_phi_2")
gen_mass_2 = Quantity("gen_mass_2")
gen_pdgid_2 = Quantity("gen_pdgid_2")

gen_mm_pair_mass = Quantity("gen_mm_pair_mass")

# sample flags
is_data = Quantity("is_data")
is_embedding = Quantity("is_embedding")
is_ttbar = Quantity("is_ttbar")
is_dyjets = Quantity("is_dyjets")
is_wjets = Quantity("is_wjets")
is_diboson = Quantity("is_diboson")

# Muon weights
reco_wgt_mu_1 = Quantity("reco_wgt_mu_1")
reco_wgt_mu_2 = Quantity("reco_wgt_mu_2")
id_wgt_mu_1 = Quantity("id_wgt_mu_1")
id_wgt_mu_2 = Quantity("id_wgt_mu_2")
iso_wgt_mu_1 = Quantity("iso_wgt_mu_1")
iso_wgt_mu_2 = Quantity("iso_wgt_mu_2")
