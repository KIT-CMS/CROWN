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
        self.output.shift(
            name
        )  # crashes on purpose if output is None. This method should not be called for a producer without output

    def writecall(self, config, shift=""):
        config[shift]["output"] = (
            "" if self.output == None else self.output.get_leaf(shift)
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
            for shift in self.output.shifts:
                calls.append(self.writecall(config, shift))
        return calls


class producer_vec(producer):
    def __init__(self, call, inputs, output, vec_config):
        super().__init__(call, inputs, output)
        self.vec_config = vec_config

    def writecalls(self, config):
        basecall = self.call
        calls = []
        shifts = [""]
        if self.output != None:
            shifts.extend(self.output.shifts)
        for shift in shifts:
            for val in config[shift][self.vec_config]:
                self.call = basecall.format_map(SafeDict({self.vec_config: val}))
                calls.append(self.writecall(config, shift))
        self.call = basecall
        return calls


prod1 = producer(
    "phyticsd::FilterID(auto {df}, const std::string {output}, std::string isolationName)",
    [],
    q.pt_1,
)
prod2 = producer(
    "FilterID(auto {df}, const std::string {output}, std::string isolationName)",
    [],
    q.pt_2,
)
prod3 = producer(
    'FilterID({df}, "{output}", {input_coll}, {ptcut})', [q.pt_1, q.pt_2], q.m_vis
)

MetFilter = producer_vec(
    'metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    [],
    None,
    "met_filters",
)
