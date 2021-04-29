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
        for q in self.inputs:
            q.children.append(self.output)

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
            "" if self.output == None else '{"' + '","'.join([x.get_leaf(shift) for x in self.output]) + '"}' if isinstance(self.output, list) else self.output.get_leaf(shift)
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
            list_of_shifts = self.output[0].shifts if isinstance(self.output, list) else self.output.shifts #all entries must have same shifts
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
            #check that all config lists (and output if applicable) have same length
            l = len(config[shift][self.vec_configs[0]])
            for key in self.vec_configs:
                if l!=len(config[shift][key]):
                    print("Following lists in config must have same length: %s, %s"%(self.vec_configs[0], key))
                    raise Exception
            if self.output!=None and len(self.output)!=l:
                    print("VectorProducer expects either no output or same amount as entries in config lists (e.g. %s)!"%self.vec_configs[0])
            for i in range(l):
                helper_dict = {}
                for key in self.vec_configs:
                    helper_dict[key] = config[shift][key][i]
                if self.output!=None:
                    helper_dict["output"] = self.output[i].get_leaf(shift)
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, shift))
        self.call = basecall
        return calls

class ProducerGroup():
    PG_count = 1 #counter for internal quantities used by ProducerGroups

    def __init__(self, call, inputs, output, subproducers):
        self.call = call
        self.inputs = inputs
        self.output = output
        self.producers = subproducers
        #link producers
        for p in self.producers:
            #skip producers without output
            if not "output" in p.call:
                continue
            #check that output quantities of subproducers are not yet filled
            if p.output!=None:
                print("Output of subproducers must be None!")
                raise Exception
            #create quantities that are produced by subproducers and then collected by the final call of the producer group
            p.output = q.Quantity("PG_internal_quantity%i"%self.__class__.PG_count) #quantities of vector producers will be duplicated later on when config is known
            self.__class__.PG_count += 1
            self.inputs.append(p.output)
        #treat own collection function as subproducer
        self.producers.append(Producer(self.call, self.inputs, self.output))

    def shift(self, name):
        for producer in self.producers:
            producer.shift(name)

    def writecalls(self, config):
        calls = []
        for producer in self.producers:
            #duplicate outputs of vector subproducers
            if hasattr(producer, "vec_configs"):
                producer.output = [producer.output]
                for i in range(len(config[""][producer.vec_configs[0]]) - 1):
                    producer.output.append(producer.output[0].copy("PG_internal_quantity%i"%self.__class__.PG_count))
                    self.__class__.PG_count += 1
                    self.inputs.append(producer.output[-1])
            #retrieve calls of subproducers
            calls.extend(producer.writecalls(config))
        return calls


prod1 = Producer(
    "phyticsd::FilterID(auto {df}, const std::string {output}, std::string isolationName)",
    [],
    q.pt_1,
)
prod2 = Producer(
    "FilterID(auto {df}, const std::string {output}, std::string isolationName)",
    [],
    q.pt_2,
)
prod3 = Producer(
    'FilterID({df}, "{output}", {input_coll}, {ptcut})', [q.pt_1, q.pt_2], q.m_vis
)

MetFilter = VectorProducer(
    'metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    [],
    None,
    ["met_filters"],
)

TauPtCut = Producer('physicsobject::CutPt({df}, "{input}", "{output}", {ptcut})', [q.Tau_pt], None)
TauEtaCut = Producer('physicsobject::CutEta({df}, "{input}", "{output}", {etacut})', [q.Tau_eta], None)
TauDzCut = Producer('physicsobject::CutDz({df}, "{input}", "{output}", {ptcut})', [q.Tau_dz], None)
TauIDFilters = VectorProducer('physicsobject::tau::FilterTauID({df}, "{output}", "{tau_id}", {tau_id_idx})', [], None, ["tau_id", "tau_id_idx"])
GoodTaus = ProducerGroup('physicsobject::CombineMasks({df}, "{output}", {input_coll})', [], q.good_taus_mask, [TauPtCut, TauEtaCut, TauDzCut, TauIDFilters])
