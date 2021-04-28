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
            shifts.extend(self.output.shifts)
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
                    helper_dict["output"] = self.output[i]
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, shift))
        self.call = basecall
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
