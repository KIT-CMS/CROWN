import code_generation.quantity as q
import logging

log = logging.getLogger(__name__)


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


class Producer:
    def __init__(self, name, call, input, output, scopes):
        log.debug("Setting up a new producer {}".format(name))

        # sanity checks
        if not isinstance(input, list) and not isinstance(input, dict):
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
        # if input not given as dict and therfore not scope specific transform into dict with all scopes
        if not isinstance(self.input, dict):
            inputdict = {}
            for scope in self.scopes:
                inputdict[scope] = self.input
            self.input = inputdict
        # keep track of variable dependencies
        if self.output != None:
            for scope in self.scopes:
                for input_quantity in self.input[scope]:
                    for output_quantity in self.output:
                        input_quantity.adopt(output_quantity, scope)
        log.debug("-----------------------------------------")
        log.debug("| Producer: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        for scope in self.scopes:
            if self.input[scope] == None:
                log.debug("| Inputs ({}): None".format(scope))
            else:
                log.debug(
                    "| Inputs ({}): {}".format(
                        scope, [input.name for input in self.input[scope]]
                    )
                )
        if self.output == None:
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug("| scopes: {}".format(self.scopes))
        log.debug("-----------------------------------------")

    # Check if a output_quantity is already used as an output by
    # another producer within the same scope.
    # If this occurs, a Exception is thrown, since this is not possible with dataframes

    def reserve_output(self, scope):
        if self.output != None:
            for output_quantity in self.output:
                output_quantity.reserve_scope(scope)

    def shift(self, name, scope="global"):
        if not scope in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if self.output is None:
            log.error(
                "Exception (%s): output None cannot be shifted ! How did you end up here ?"
                % name
            )
            raise Exception
        for entry in self.output:
            entry.shift(name, scope)

    def ignore_shift(self, name, scope="global"):
        if not scope in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if self.output is None:
            log.error(
                "Exception (%s): output None cannot be shifted ! How did you end up here ?"
                % name
            )
            raise Exception
        for entry in self.output:
            entry.ignore_shift(name, scope)

    def writecall(self, config, scope, shift=""):
        if self.output == None:
            config[shift][scope]["output"] = ""
            config[shift][scope]["output_vec"] = ""
        else:
            config[shift][scope]["output"] = (
                '"' + '","'.join([x.get_leaf(shift, scope) for x in self.output]) + '"'
            )
            config[shift][scope]["output_vec"] = (
                '{"'
                + '","'.join([x.get_leaf(shift, scope) for x in self.output])
                + '"}'
            )
        config[shift][scope]["input"] = (
            '"'
            + '", "'.join([x.get_leaf(shift, scope) for x in self.input[scope]])
            + '"'
        )
        config[shift][scope]["input_vec"] = (
            '{"'
            + '","'.join([x.get_leaf(shift, scope) for x in self.input[scope]])
            + '"}'
        )
        config[shift][scope]["df"] = "{df}"
        config[shift][scope]["vec_open"] = "{vec_open}"
        config[shift][scope]["vec_close"] = "{vec_close}"
        return self.call.format(
            **config[shift][scope]
        )  # use format (not format_map here) such that missing config entries cause an error

    def writecalls(self, config, scope):
        if scope not in self.scopes:
            log.error(
                "Exception (%s): Tried to use producer in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
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
            n_versions = len(config[shift][scope][self.vec_configs[0]])
            for key in self.vec_configs:
                if n_versions != len(config[shift][scope][key]):
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
                    helper_dict[key] = config[shift][scope][key][i]
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


class BaseFilter(Producer):
    def __init__(self, name, call, input, scopes):
        super().__init__(name, call, input, None, scopes)

    def writecall(self, config, scope, shift=""):
        log.critical("{}: Filters do not support method writecall!".format(self.name))
        raise Exception

    def writecalls(self, config, scope):
        inputs = []
        for quantity in self.input[scope]:
            inputs.extend(quantity.get_leafs_of_scope(scope))
        formatdict = {}
        formatdict["input"] = '"' + '", "'.join(inputs) + '"'
        formatdict["input_vec"] = '{"' + '","'.join(inputs) + '"}'
        formatdict["df"] = "{df}"
        return [
            self.call.format(**formatdict)
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
        # if input not given as dict and therfore not scope specific transform into dict with all scopes
        if not isinstance(self.input, dict):
            inputdict = {}
            for scope in self.scopes:
                inputdict[scope] = self.input
            self.input = inputdict
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
                    for scope in subproducer.scopes:
                        for quantity in subproducer.input[scope]:
                            quantity.adopt(subproducer.output[0], scope)
                for scope in self.scopes:
                    for output_quantity in subproducer.output:
                        self.input[scope].append(output_quantity)
            # treat own collection function as subproducer
            self.setup_own_producer()
        log.debug("-----------------------------------------")
        log.debug("| ProducerGroup: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        for scope in self.scopes:
            if self.input[scope] == None:
                log.debug("| Inputs ({}): None".format(scope))
            else:
                log.debug(
                    "| Inputs ({}): {}".format(
                        scope, [input.name for input in self.input[scope]]
                    )
                )
        if self.output == None:
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug(
            "| Producers: {}".format([producer.name for producer in self.producers])
        )
        log.debug("| scopes: {}".format(self.scopes))
        log.debug("-----------------------------------------")

    def setup_own_producer(self):
        self.producers.append(
            Producer(self.name, self.call, self.input, self.output, self.scopes)
        )

    # for a producer group, step iteratively
    # through the subproducers and reserve the output there

    def reserve_output(self, scope):
        for subproducer in self.producers:
            subproducer.reserve_output(scope)

    def shift(self, name, scope="global"):
        for producer in self.producers:
            producer.shift(name, scope)

    def ignore_shift(self, name, scope="global"):
        for producer in self.producers:
            producer.ignore_shift(name, scope)

    def writecalls(self, config, scope):
        calls = []
        for producer in self.producers:
            # duplicate outputs of vector subproducers if they were generated automatically
            if self.call != None and hasattr(producer, "vec_configs"):
                for i in range(len(config[""][scope][producer.vec_configs[0]]) - 1):
                    producer.output.append(
                        producer.output[0].copy(
                            "PG_internal_quantity_%i" % self.__class__.PG_count
                        )
                    )
                    self.__class__.PG_count += 1
                    self.input[scope].append(producer.output[-1])
            # retrieve calls of subproducers
            calls.extend(producer.writecalls(config, scope))
        return calls


class Filter(ProducerGroup):
    def __init__(self, name, call, input, scopes, subproducers):
        self.__class__.PG_count = ProducerGroup.PG_count
        super().__init__(name, call, input, None, scopes, subproducers)
        ProducerGroup.PG_count = self.__class__.PG_count

    def setup_own_producer(self):
        self.producers.append(BaseFilter(self.name, self.call, self.input, self.scopes))
