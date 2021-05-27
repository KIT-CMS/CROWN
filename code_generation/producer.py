import code_generation.quantities as q
import logging

log = logging.getLogger(__name__)


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


class Producer:
    def __init__(self, name, call, input, output, scopes):
        log.debug("Setting up a new producer {}".format(name))

        # sanity checks
        if not isinstance(input, list):
            print("Exception (%s): Argument 'input' must be a list!" % name)
            raise Exception
        if not isinstance(output, list) and output != None:
            print("Exception (%s): Argument 'output' must be a list or None!" % name)
            raise Exception
        self.name = name
        self.call = call
        self.input = input
        self.output = output
        self.scopes = scopes
        self.output_checked = False
        # keep track of variable dependencies
        if self.output != None:
            for input_quantity in self.input:
                for scope in self.scopes:
                    for output_quantity in self.output:
                        input_quantity.adopt(output_quantity, scope)
        log.debug("-----------------------------------------")
        log.debug("| Producer: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        if self.input == None:
            log.debug("| Inputs: None")
        else:
            log.debug("| Inputs: {}".format([input.name for input in self.input]))
        if self.output == None:
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug("| scopes: {}".format(self.scopes))
        log.debug("-----------------------------------------")

    def check_output(self):
        if self.output != None and not self.output_checked:
            for scope in self.scopes:
                for output_quantity in self.output:
                    output_quantity.check_scope(scope)
        self.output_checked = True

    def shift(self, name, scope="global"):
        if not scope in self.scopes:
            print(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if self.output is None:
            print(
                "Exception (%s): output None cannot be shifted ! How did you end up here ?"
                % name
            )
            raise Exception
        for entry in self.output:
            entry.shift(
                name, scope
            )  # crashes on purpose if output is None. This method should not be called for a producer without output

    def writecall(self, config, scope, shift=""):
        if self.output == None:
            config[shift]["output"] = ""
            config[shift]["output_vec"] = ""
        else:
            config[shift]["output"] = (
                '"' + '","'.join([x.get_leaf(shift, scope) for x in self.output]) + '"'
            )
            config[shift]["output_vec"] = (
                '{"'
                + '","'.join([x.get_leaf(shift, scope) for x in self.output])
                + '"}'
            )
        config[shift]["input"] = (
            '"' + '", "'.join([x.get_leaf(shift, scope) for x in self.input]) + '"'
        )
        config[shift]["input_vec"] = (
            '{"' + '","'.join([x.get_leaf(shift, scope) for x in self.input]) + '"}'
        )
        config[shift]["df"] = "{df}"
        return self.call.format(
            **config[shift]
        )  # use format (not format_map here) such that missing config entries cause an error

    def writecalls(self, config, scope):
        calls = [self.writecall(config, scope)]
        if self.output != None:
            list_of_shifts = self.output[0].get_shifts(
                scope
            )  # all entries must have same shifts
            for shift in list_of_shifts:
                calls.append(self.writecall(config, scope, shift))
        return calls


class VectorProducer(Producer):
    def __init__(self, name, call, input, output, scopes, vec_configs):
        self.name = name
        super().__init__(name, call, input, output, scopes)
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
                    helper_dict["output"] = (
                        '"' + self.output[i].get_leaf(shift, scope) + '"'
                    )
                    helper_dict["output_vec"] = (
                        '{"' + self.output[i].get_leaf(shift, scope) + '"}'
                    )
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, scope, shift))
        self.call = basecall
        return calls


class Filter(Producer):
    def __init__(self, name, call, input, scopes):
        super().__init__(name, call, input, None, scopes)

    def writecall(self, config, scope, shift=""):
        log.critical("{}: Filters do not support method writecall!".format(self.name))
        raise Exception

    def writecalls(self, config, scope):
        config[shift]["output"] = ""
        config[shift]["output_vec"] = ""
        inputs = []
        for quantity in self.input:
            inputs.extend(quantity.get_leafs_of_scope(scope))
        config[shift]["input"] = '"' + '", "'.join(inputs) + '"'
        config[shift]["input_vec"] = '{"' + '","'.join(inputs) + '"}'
        config[shift]["df"] = "{df}"
        return [
            self.call.format(**config[shift])
        ]  # use format (not format_map here) such that missing config entries cause an error


class ProducerGroup:
    PG_count = 1  # counter for internal quantities used by ProducerGroups

    def __init__(self, name, call, input, output, scopes, subproducers):
        self.name = name
        self.call = call
        self.input = input
        self.output = output
        self.producers = subproducers
        self.scopes = scopes
        # If call is provided, this is supposed to consume output of subproducers. Creating these internal products below:
        if self.call != None:
            for subproducer in self.producers:
                # skip producers without output
                if subproducer.output == None:
                    continue
                # if the subproducer does not have an output quantity, we assign an internally tracked quantity
                if subproducer.output == []:
                    # create quantities that are produced by subproducers and then collected by the final call of the producer group
                    subproducer.output = [
                        q.Quantity("PG_internal_quantity_%i" % self.__class__.PG_count)
                    ]  # quantities of vector producers will be duplicated later on when config is known
                    self.__class__.PG_count += 1
                    if isinstance(subproducer, ProducerGroup):
                        subproducer.producers[-1].output = subproducer.output
                    for quantity in subproducer.input:
                        for scope in subproducer.scopes:
                            quantity.adopt(subproducer.output[0], scope)
                for output_quantity in subproducer.output:
                    self.input.append(output_quantity)
            # treat own collection function as subproducer
            self.producers.append(
                Producer(self.name, self.call, self.input, self.output, self.scopes)
            )
        log.debug("-----------------------------------------")
        log.debug("| ProducerGroup: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        if self.input == None:
            log.debug("| Inputs: None")
        else:
            log.debug("| Inputs: {}".format([input.name for input in self.input]))
        if self.output == None:
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug(
            "| Producers: {}".format([producer.name for producer in self.producers])
        )
        log.debug("| scopes: {}".format(self.scopes))
        log.debug("-----------------------------------------")

    def check_output(self):
        for subproducer in self.producers:
            subproducer.check_output()

    def shift(self, name, scope="global"):
        for producer in self.producers:
            producer.shift(name, scope)

    def writecalls(self, config, scope):
        calls = []
        for producer in self.producers:
            # duplicate outputs of vector subproducers if they were generated automatically
            if self.call != None and hasattr(producer, "vec_configs"):
                for i in range(len(config[""][producer.vec_configs[0]]) - 1):
                    producer.output.append(
                        producer.output[0].copy(
                            "PG_internal_quantity_%i" % self.__class__.PG_count
                        )
                    )
                    self.__class__.PG_count += 1
                    self.input.append(producer.output[-1])
            # retrieve calls of subproducers
            calls.extend(producer.writecalls(config, scope))
        return calls


MetFilter = VectorProducer(
    name="MetFilter",
    call='metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    input=[],
    output=None,
    scopes=["global"],
    vec_configs=["met_filters"],
)

TauPtCut = Producer(
    name="TauPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_tau_pt})",
    input=[q.Tau_pt],
    output=[],
    scopes=["global"],
)
TauEtaCut = Producer(
    name="TauEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_tau_eta})",
    input=[q.Tau_eta],
    output=[],
    scopes=["global"],
)
TauDzCut = Producer(
    name="TauDzCut",
    call="physicsobject::CutDz({df}, {input}, {output}, {max_tau_dz})",
    input=[q.Tau_dz],
    output=[],
    scopes=["global"],
)
TauIDFilters = VectorProducer(
    name="TauIDFilters",
    call='physicsobject::tau::FilterTauID({df}, {output}, "{tau_id}", {tau_id_idx})',
    input=[],
    output=[],
    scopes=["global"],
    vec_configs=["tau_id", "tau_id_idx"],
)
GoodTaus = ProducerGroup(
    name="GoodTaus",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_taus_mask],
    scopes=["global"],
    subproducers=[TauPtCut, TauEtaCut, TauDzCut, TauIDFilters],
)

MuonPtCut = Producer(
    name="MuonPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_muon_pt})",
    input=[q.Muon_pt],
    output=[],
    scopes=["global"],
)
MuonEtaCut = Producer(
    name="MuonEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_muon_eta})",
    input=[q.Muon_eta],
    output=[],
    scopes=["global"],
)
MuonIDFilter = Producer(
    name="MuonIDFilter",
    call='physicsobject::muon::FilterID({df}, {output}, "{muon_id}")',
    input=[],
    output=[],
    scopes=["global"],
)
MuonIsoFilter = Producer(
    name="MuonIsoFilter",
    call="physicsobject::muon::FilterIsolation({df}, {output}, {input}, {muon_iso_cut})",
    input=[q.Muon_iso],
    output=[],
    scopes=["global"],
)
GoodMuons = ProducerGroup(
    name="GoodMuons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_muons_mask],
    scopes=["global"],
    subproducers=[MuonPtCut, MuonEtaCut, MuonIDFilter, MuonIsoFilter],
)

VetoElectronPtCut = Producer(
    name="ElectronPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_VetoElectron_pt})",
    input=[q.Electron_pt],
    output=[],
    scopes=["global"],
)
VetoElectronEtaCut = Producer(
    name="ElectronEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_VetoElectron_eta})",
    input=[q.Electron_eta],
    output=[],
    scopes=["global"],
)
VetoElectronIDFilter = Producer(
    name="ElectronIDFilter",
    call='physicsobject::electron::FilterID({df}, {output}, "{VetoElectron_id}")',
    input=[],
    output=[],
    scopes=["global"],
)
VetoElectronIsoFilter = Producer(
    name="ElectronIsoFilter",
    call="physicsobject::electron::FilterIsolation({df}, {output}, {input}, {max_VetoElectron_iso})",
    input=[q.Electron_iso],
    output=[],
    scopes=["global"],
)
VetoElectrons = ProducerGroup(
    name="VetoElectrons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[],
    scopes=["global"],
    subproducers=[
        VetoElectronPtCut,
        VetoElectronEtaCut,
        VetoElectronIDFilter,
        VetoElectronIsoFilter,
    ],
)

GoodElectronsVeto = ProducerGroup(
    name="GoodElectronsVeto",
    call="physicsobject::LeptonVetoFlag({df}, {output}, {input})",
    input=[],
    output=[q.electron_veto_flag],
    scopes=["global"],
    subproducers=[VetoElectrons],
)


RequireObjects = VectorProducer(
    name="RequireObjects",
    call='physicsobject::FilterObjects({df}, "{require_candidate}", {require_candidate_number}, "{require_candidate}")',
    input=[],
    output=[],
    scopes=["global"],
    vec_configs=["require_candidate", "require_candidate_number"],
)

MTPairSelection = Producer(
    name="MTPairSelection",
    call="pairselection::mutau::PairSelection({df}, {input_vec}, {output})",
    input=[
        q.Tau_pt,
        q.Tau_IDraw,
        q.Muon_pt,
        q.Muon_iso,
        q.good_taus_mask,
        q.good_muons_mask,
    ],
    output=[q.ditaupair],
    scopes=["mt"],
)

GoodMTPairFilter = Producer(
    name="GoodMTPairFilter",
    call='pairselection::filterGoodPairs({df}, {input}, "GoodMuTauPairs")',
    input=[q.ditaupair],
    output=None,
    scopes=["mt"],
)

LVMu1 = Producer(
    name="LVMu1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[q.ditaupair, q.Muon_pt, q.Muon_eta, q.Muon_phi, q.Muon_mass],
    output=[q.p4_1],
    scopes=["mt"],
)
LVMu2 = Producer(
    name="LVMu2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[q.ditaupair, q.Muon_pt, q.Muon_eta, q.Muon_phi, q.Muon_mass],
    output=[q.p4_2],
    scopes=["mt"],
)
LVTau1 = Producer(
    name="LVTau1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[q.ditaupair, q.Tau_pt, q.Tau_eta, q.Tau_phi, q.Tau_mass],
    output=[q.p4_1],
    scopes=["mt"],
)
LVTau2 = Producer(
    name="LVTau2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[q.ditaupair, q.Tau_pt, q.Tau_eta, q.Tau_phi, q.Tau_mass],
    output=[q.p4_2],
    scopes=["mt"],
)
