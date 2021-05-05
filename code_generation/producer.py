import code_generation.quantities as q


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


class Producer:
    def __init__(self, name, call, inputs, output, scopes):
        self.name = name
        self.call = call
        self.inputs = inputs
        self.output = output
        self.scopes = scopes

        # keep track of variable dependencies
        if self.output != None:
            for quantity in self.inputs:
                for scope in self.scopes:
                    quantity.adopt(self.output, scope)

    def shift(self, name, scope="global"):
        if not scope in self.scopes:
            print(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if isinstance(self.output, list):
            for entry in self.output:
                entry.shift(name, scope)
        else:
            self.output.shift(
                name, scope
            )  # crashes on purpose if output is None. This method should not be called for a producer without output

    def writecall(self, config, scope, shift=""):
        config[shift]["output"] = (
            ""
            if self.output == None
            else '{"'
            + '","'.join([x.get_leaf(shift, scope) for x in self.output])
            + '"}'
            if isinstance(self.output, list)
            else self.output.get_leaf(shift, scope)
        )
        config[shift]["input"] = ",".join(
            [x.get_leaf(shift, scope) for x in self.inputs]
        )
        config[shift]["input_coll"] = (
            '{"' + '","'.join([x.get_leaf(shift, scope) for x in self.inputs]) + '"}'
        )
        config[shift]["df"] = "{df}"
        return self.call.format(
            **config[shift]
        )  # use format (not format_map here) such that missing config entries cause an error

    def writecalls(self, config, scope):
        calls = [self.writecall(config, scope)]
        if self.output != None:
            list_of_shifts = (
                self.output[0].get_shifts(scope)
                if isinstance(self.output, list)
                else self.output.get_shifts(scope)
            )  # all entries must have same shifts
            for shift in list_of_shifts:
                calls.append(self.writecall(config, scope, shift))
        return calls


class VectorProducer(Producer):
    def __init__(self, name, call, inputs, output, scopes, vec_configs):
        self.name = name
        super().__init__(name, call, inputs, output, scopes)
        self.vec_configs = vec_configs

    def writecalls(self, config, scope):
        basecall = self.call
        calls = []
        shifts = [""]
        if self.output != None:
            shifts.extend(self.output[0].get_shifts(scope))
        for shift in shifts:
            # check that all config lists (and output if applicable) have same length
            n_versions = len(config[shift][self.vec_configs[0]])
            for key in self.vec_configs:
                if n_versions != len(config[shift][key]):
                    print(
                        "Following lists in config must have same length: %s, %s"
                        % (self.vec_configs[0], key)
                    )
                    raise Exception
            if self.output != None and len(self.output) != n_versions:
                print(
                    "VectorProducer expects either no output or same amount as entries in config lists (e.g. %s)!"
                    % self.vec_configs[0]
                )
            for i in range(n_versions):
                helper_dict = {}
                for key in self.vec_configs:
                    helper_dict[key] = config[shift][key][i]
                if self.output != None:
                    helper_dict["output"] = self.output[i].get_leaf(shift, scope)
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, scope, shift))
        self.call = basecall
        return calls


class ProducerGroup:
    PG_count = 1  # counter for internal quantities used by ProducerGroups

    def __init__(self, name, call, inputs, output, scopes, subproducers):
        self.name = name
        self.call = call
        self.inputs = inputs
        self.output = output
        self.producers = subproducers
        self.scopes = scopes
        # If call is provided, this is supposed to consume output of subproducers. Creating these internal products below:
        if self.call != None:
            for producer in self.producers:
                # check that output quantities of subproducers are not yet filled
                if producer.output != None:
                    print("Output of subproducers must be None!")
                    raise Exception
                # skip producers without output
                if not "output" in producer.call:
                    continue
                # create quantities that are produced by subproducers and then collected by the final call of the producer group
                producer.output = q.Quantity(
                    "PG_internal_quantity_%i" % self.__class__.PG_count
                )  # quantities of vector producers will be duplicated later on when config is known
                for quantity in producer.inputs:
                    for scope in producer.scopes:
                        quantity.adopt(producer.output, scope)
                self.__class__.PG_count += 1
                self.inputs.append(producer.output)
            # treat own collection function as subproducer
            self.producers.append(
                Producer(self.name, self.call, self.inputs, self.output, self.scopes)
            )

    def shift(self, name, scope="global"):
        for producer in self.producers:
            producer.shift(name, scope)

    def writecalls(self, config, scope):
        calls = []
        for producer in self.producers:
            # duplicate outputs of vector subproducers if they were generated automatically
            if self.call != None and hasattr(producer, "vec_configs"):
                producer.output = [producer.output]
                for i in range(len(config[""][producer.vec_configs[0]]) - 1):
                    producer.output.append(
                        producer.output[0].copy(
                            "PG_internal_quantity_%i" % self.__class__.PG_count
                        )
                    )
                    self.__class__.PG_count += 1
                    self.inputs.append(producer.output[-1])
            # retrieve calls of subproducers
            calls.extend(producer.writecalls(config, scope))
        return calls


MetFilter = VectorProducer(
    name="MetFilter",
    call='metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    inputs=[],
    output=None,
    scopes=["global"],
    vec_configs=["met_filters"],
)

TauPtCut = Producer(
    name="TauPtCut",
    call='physicsobject::CutPt({df}, "{input}", "{output}", {min_tau_pt})',
    inputs=[q.Tau_pt],
    output=None,
    scopes=["global"],
)
TauEtaCut = Producer(
    name="TauEtaCut",
    call='physicsobject::CutEta({df}, "{input}", "{output}", {max_tau_eta})',
    inputs=[q.Tau_eta],
    output=None,
    scopes=["global"],
)
TauDzCut = Producer(
    name="TauDzCut",
    call='physicsobject::CutDz({df}, "{input}", "{output}", {max_tau_dz})',
    inputs=[q.Tau_dz],
    output=None,
    scopes=["global"],
)
TauIDFilters = VectorProducer(
    name="TauIDFilters",
    call='physicsobject::tau::FilterTauID({df}, "{output}", "{tau_id}", {tau_id_idx})',
    inputs=[],
    output=None,
    scopes=["global"],
    vec_configs=["tau_id", "tau_id_idx"],
)
GoodTaus = ProducerGroup(
    name="GoodTaus",
    call='physicsobject::CombineMasks({df}, "{output}", {input_coll})',
    inputs=[],
    output=q.good_taus_mask,
    scopes=["global"],
    subproducers=[TauPtCut, TauEtaCut, TauDzCut, TauIDFilters],
)

MuonPtCut = Producer(
    name="MuonPtCut",
    call='physicsobject::CutPt({df}, "{input}", "{output}", {min_muon_pt})',
    inputs=[q.Muon_pt],
    output=None,
    scopes=["global"],
)
MuonEtaCut = Producer(
    name="MuonEtaCut",
    call='physicsobject::CutEta({df}, "{input}", "{output}", {max_muon_eta})',
    inputs=[q.Muon_eta],
    output=None,
    scopes=["global"],
)
MuonIDFilter = Producer(
    name="MuonIDFilter",
    call='physicsobject::muon::FilterID({df}, "{output}", "{muon_id}")',
    inputs=[],
    output=None,
    scopes=["global"],
)
MuonIsoFilter = Producer(
    name="MuonIsoFilter",
    call='physicsobject::muon::FilterIsolation({df}, "{output}", "{input}", {muon_iso_cut})',
    inputs=[q.Muon_iso],
    output=None,
    scopes=["global"],
)
GoodMuons = ProducerGroup(
    name="GoodMuons",
    call='physicsobject::CombineMasks({df}, "{output}", {input_coll})',
    inputs=[],
    output=q.good_muons_mask,
    scopes=["global"],
    subproducers=[MuonPtCut, MuonEtaCut, MuonIDFilter, MuonIsoFilter],
)

RequireObjects = VectorProducer(
    name="RequireObjects",
    call='physicsobject::FilterObjects({df}, "{require_candidate}", {require_candidate_number}, "{require_candidate}")',
    inputs=[],
    output=None,
    scopes=["global"],
    vec_configs=["require_candidate", "require_candidate_number"],
)

MTPairSelection = Producer(
    name="MTPairSelection",
    call='pairselection::mutau::PairSelection({df}, {input_coll}, "{output}")',
    inputs=[
        q.Tau_pt,
        q.Tau_IDraw,
        q.Muon_pt,
        q.Muon_iso,
        q.good_taus_mask,
        q.good_muons_mask,
    ],
    output=q.ditaupair,
    scopes=["mt"],
)

GoodMTPairFilter = Producer(
    name="GoodMTPairFilter",
    call='pairselection::filterGoodPairs({df}, "{input}", "GoodMuTauPairs")',
    inputs=[q.ditaupair],
    output=None,
    scopes=["mt"],
)

LVMu1 = Producer(
    name="LVMu1",
    call='lorentzvectors::build({df}, {input_coll}, 0, "{output}")',
    inputs=[q.ditaupair, q.Muon_pt, q.Muon_eta, q.Muon_phi, q.Muon_mass],
    output=q.p4_1,
    scopes=["mt"],
)
LVMu2 = Producer(
    name="LVMu2",
    call='lorentzvectors::build({df}, {input_coll}, 1, "{output}")',
    inputs=[q.ditaupair, q.Muon_pt, q.Muon_eta, q.Muon_phi, q.Muon_mass],
    output=q.p4_2,
    scopes=["mt"],
)
LVTau1 = Producer(
    name="LVTau1",
    call='lorentzvectors::build({df}, {input_coll}, 0, "{output}")',
    inputs=[q.ditaupair, q.Tau_pt, q.Tau_eta, q.Tau_phi, q.Tau_mass],
    output=q.p4_1,
    scopes=["mt"],
)
LVTau2 = Producer(
    name="LVTau2",
    call='lorentzvectors::build({df}, {input_coll}, 1, "{output}")',
    inputs=[q.ditaupair, q.Tau_pt, q.Tau_eta, q.Tau_phi, q.Tau_mass],
    output=q.p4_2,
    scopes=["mt"],
)

pt_1 = Producer(
    name="pt_1",
    call='quantities::pt({df}, varSet, "{output}", "{input}")',
    inputs=[q.p4_1],
    output=q.pt_1,
    scopes=["mt"],
)
pt_2 = Producer(
    name="pt_2",
    call='quantities::pt({df}, varSet, "{output}", "{input}")',
    inputs=[q.p4_2],
    output=q.pt_2,
    scopes=["mt"],
)
eta_1 = Producer(
    name="eta_1",
    call='quantities::eta({df}, varSet, "{output}", "{input}")',
    inputs=[q.p4_1],
    output=q.eta_1,
    scopes=["mt"],
)
eta_2 = Producer(
    name="eta_2",
    call='quantities::eta({df}, varSet, "{output}", "{input}")',
    inputs=[q.p4_2],
    output=q.eta_2,
    scopes=["mt"],
)
phi_1 = Producer(
    name="phi_1",
    call='quantities::phi({df}, varSet, "{output}", "{input}")',
    inputs=[q.p4_1],
    output=q.phi_1,
    scopes=["mt"],
)
phi_2 = Producer(
    name="phi_2",
    call='quantities::phi({df}, varSet, "{output}", "{input}")',
    inputs=[q.p4_2],
    output=q.phi_2,
    scopes=["mt"],
)
UnrollLV1 = ProducerGroup(
    name="UnrollLV1",
    call=None,
    inputs=None,
    output=None,
    scopes=["mt"],
    subproducers=[pt_1, eta_1, phi_1],
)
UnrollLV2 = ProducerGroup(
    name="UnrollLV2",
    call=None,
    inputs=None,
    output=None,
    scopes=["mt"],
    subproducers=[pt_2, eta_2, phi_2],
)

m_vis = Producer(
    name="m_vis",
    call='quantities::m_vis({df}, varSet, "{output}", {input_coll})',
    inputs=[q.p4_1, q.p4_2],
    output=q.m_vis,
    scopes=["mt"],
)
DiTauPairQuantities = ProducerGroup(
    name="DiTauPairQuantities",
    call=None,
    inputs=None,
    output=None,
    scopes=["mt"],
    subproducers=[UnrollLV1, UnrollLV2, m_vis],
)
