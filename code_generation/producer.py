import code_generation.quantities as q


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


class Producer:
    def __init__(self, call, inputs, output):
        self.call = call
        self.inputs = inputs
        self.output = output

        # keep track of variable dependencies
        if self.output != None:
            for qy in self.inputs:
                qy.children.append(self.output)

    def shift(self, name):
        if isinstance(self.output, list):
            for entry in self.output:
                entry.shift(name)
        else:
            self.output.shift(
                name
            )  # crashes on purpose if output is None. This method should not be called for a producer without output

    def writecall(self, config, shift=""):
        config[shift]["output"] = (
            ""
            if self.output == None
            else '{"' + '","'.join([x.get_leaf(shift) for x in self.output]) + '"}'
            if isinstance(self.output, list)
            else self.output.get_leaf(shift)
        )
        config[shift]["input"] = ",".join([x.get_leaf(shift) for x in self.inputs])
        config[shift]["input_coll"] = (
            '{"' + '","'.join([x.get_leaf(shift) for x in self.inputs]) + '"}'
        )
        config[shift]["df"] = "{df}"
        return self.call.format(
            **config[shift]
        )  # use format (not format_map here) such that missing config entries cause an error

    def writecalls(self, config):
        calls = [self.writecall(config)]
        if self.output != None:
            list_of_shifts = (
                self.output[0].shifts
                if isinstance(self.output, list)
                else self.output.shifts
            )  # all entries must have same shifts
            for shift in list_of_shifts:
                calls.append(self.writecall(config, shift))
        return calls


class VectorProducer(Producer):
    def __init__(self, call, inputs, output, vec_configs):
        super().__init__(call, inputs, output)
        self.vec_configs = vec_configs

    def writecalls(self, config):
        basecall = self.call
        calls = []
        shifts = [""]
        if self.output != None:
            shifts.extend(self.output[0].shifts)
        for shift in shifts:
            # check that all config lists (and output if applicable) have same length
            l = len(config[shift][self.vec_configs[0]])
            for key in self.vec_configs:
                if l != len(config[shift][key]):
                    print(
                        "Following lists in config must have same length: %s, %s"
                        % (self.vec_configs[0], key)
                    )
                    raise Exception
            if self.output != None and len(self.output) != l:
                print(
                    "VectorProducer expects either no output or same amount as entries in config lists (e.g. %s)!"
                    % self.vec_configs[0]
                )
            for i in range(l):
                helper_dict = {}
                for key in self.vec_configs:
                    helper_dict[key] = config[shift][key][i]
                if self.output != None:
                    helper_dict["output"] = self.output[i].get_leaf(shift)
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, shift))
        self.call = basecall
        return calls


class ProducerGroup:
    PG_count = 1  # counter for internal quantities used by ProducerGroups

    def __init__(self, call, inputs, output, subproducers):
        self.call = call
        self.inputs = inputs
        self.output = output
        self.producers = subproducers
        # If call is provided, this is supposed to consume output of subproducers. Creating these internal products below:
        if self.call != None:
            for p in self.producers:
                # check that output quantities of subproducers are not yet filled
                if p.output != None:
                    print("Output of subproducers must be None!")
                    raise Exception
                # skip producers without output
                if not "output" in p.call:
                    continue
                # create quantities that are produced by subproducers and then collected by the final call of the producer group
                p.output = q.Quantity(
                    "PG_internal_quantity_%i" % self.__class__.PG_count
                )  # quantities of vector producers will be duplicated later on when config is known
                for qy in p.inputs:
                    qy.children.append(p.output)
                self.__class__.PG_count += 1
                self.inputs.append(p.output)
            # treat own collection function as subproducer
            self.producers.append(Producer(self.call, self.inputs, self.output))

    def shift(self, name):
        for producer in self.producers:
            producer.shift(name)

    def writecalls(self, config):
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
            calls.extend(producer.writecalls(config))
        return calls


MetFilter = VectorProducer(
    'metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    [],
    None,
    ["met_filters"],
)

TauPtCut = Producer(
    'physicsobject::CutPt({df}, "{input}", "{output}", {min_tau_pt})', [q.Tau_pt], None
)
TauEtaCut = Producer(
    'physicsobject::CutEta({df}, "{input}", "{output}", {max_tau_eta})',
    [q.Tau_eta],
    None,
)
TauDzCut = Producer(
    'physicsobject::CutDz({df}, "{input}", "{output}", {max_tau_dz})', [q.Tau_dz], None
)
TauIDFilters = VectorProducer(
    'physicsobject::tau::FilterTauID({df}, "{output}", "{tau_id}", {tau_id_idx})',
    [],
    None,
    ["tau_id", "tau_id_idx"],
)
GoodTaus = ProducerGroup(
    'physicsobject::CombineMasks({df}, "{output}", {input_coll})',
    [],
    q.good_taus_mask,
    [TauPtCut, TauEtaCut, TauDzCut, TauIDFilters],
)

MuonPtCut = Producer(
    'physicsobject::CutPt({df}, "{input}", "{output}", {min_muon_pt})',
    [q.Muon_pt],
    None,
)
MuonEtaCut = Producer(
    'physicsobject::CutEta({df}, "{input}", "{output}", {max_muon_eta})',
    [q.Muon_eta],
    None,
)
MuonIDFilter = Producer(
    'physicsobject::muon::FilterID({df}, "{output}", "{muon_id}")', [], None
)
MuonIsoFilter = Producer(
    'physicsobject::muon::FilterIsolation({df}, "{output}", "{input}", {muon_iso_cut})',
    [q.Muon_iso],
    None,
)
GoodMuons = ProducerGroup(
    'physicsobject::CombineMasks({df}, "{output}", {input_coll})',
    [],
    q.good_muons_mask,
    [MuonPtCut, MuonEtaCut, MuonIDFilter, MuonIsoFilter],
)

RequireObjects = VectorProducer(
    'physicsobject::FilterObjects({df}, "{require_candidate}", {require_candidate_number}, "{require_candidate}")',
    [],
    None,
    ["require_candidate", "require_candidate_number"],
)

MTPairSelection = Producer(
    'pairselection::mutau::PairSelection({df}, {input_coll}, "{output}")',
    [q.Tau_pt, q.Tau_IDraw, q.Muon_pt, q.Muon_iso, q.good_taus_mask, q.good_muons_mask],
    q.ditaupair,
)

GoodMTPairFilter = Producer(
    'pairselection::filterGoodPairs({df}, "{input}", "GoodMuTauPairs")',
    [q.ditaupair],
    None,
)

LVMu1 = Producer(
    'lorentzvectors::build({df}, {input_coll}, 0, "{output}")',
    [q.ditaupair, q.Muon_pt, q.Muon_eta, q.Muon_phi, q.Muon_mass],
    q.p4_1,
)
LVMu2 = Producer(
    'lorentzvectors::build({df}, {input_coll}, 1, "{output}")',
    [q.ditaupair, q.Muon_pt, q.Muon_eta, q.Muon_phi, q.Muon_mass],
    q.p4_2,
)
LVTau1 = Producer(
    'lorentzvectors::build({df}, {input_coll}, 0, "{output}")',
    [q.ditaupair, q.Tau_pt, q.Tau_eta, q.Tau_phi, q.Tau_mass],
    q.p4_1,
)
LVTau2 = Producer(
    'lorentzvectors::build({df}, {input_coll}, 1, "{output}")',
    [q.ditaupair, q.Tau_pt, q.Tau_eta, q.Tau_phi, q.Tau_mass],
    q.p4_2,
)

pt_1 = Producer('quantities::pt({df}, varSet, "{output}", "{input}")', [q.p4_1], q.pt_1)
pt_2 = Producer('quantities::pt({df}, varSet, "{output}", "{input}")', [q.p4_2], q.pt_2)
eta_1 = Producer(
    'quantities::eta({df}, varSet, "{output}", "{input}")', [q.p4_1], q.eta_1
)
eta_2 = Producer(
    'quantities::eta({df}, varSet, "{output}", "{input}")', [q.p4_2], q.eta_2
)
phi_1 = Producer(
    'quantities::phi({df}, varSet, "{output}", "{input}")', [q.p4_1], q.phi_1
)
phi_2 = Producer(
    'quantities::phi({df}, varSet, "{output}", "{input}")', [q.p4_2], q.phi_2
)
UnrollLV1 = ProducerGroup(None, None, None, [pt_1, eta_1, phi_1])
UnrollLV2 = ProducerGroup(None, None, None, [pt_2, eta_2, phi_2])

m_vis = Producer(
    'quantities::m_vis({df}, varSet, "{output}", {input_coll})', [q.p4_1, q.p4_2], q.m_vis
)
DiTauPairQuantities = ProducerGroup(None, None, None, [UnrollLV1, UnrollLV2, m_vis])
