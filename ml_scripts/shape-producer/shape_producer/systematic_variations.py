# -*- coding: utf-8 -*-
"""
"""
import logging
logger = logging.getLogger(__name__)
import pprint

# this helper function can be used in case the systematic variation's name ends with "Down" and "Up"
def create_systematic_variations(name, property_name, systematic_variation):
    results = [] 
    results.append(systematic_variation(name, property_name, "Down"))
    results.append(systematic_variation(name, property_name, "Up"))
    return results


# class performing the systematic variation
class SystematicVariation(object):
    def __init__(self, name, direction):
        self._name = name
        self._direction = direction

    def __str__(self):
        return "SystematicVariation(name=%s, direction=%s)" % (self._name, self._direction)

    # TODO: init name vs getter name is different, confusing!
    @property
    def name(self):
        return self._name + self._direction

    def change_histogram_name(self, h_settings, direction):
        if isinstance(h_settings["name"], list):
            h_settings["name"].append(direction)
        else:
            h_settings["name"] = [h_settings["name"], direction]
        return h_settings

    def shifted_root_objects(self, h_settings):
        return h_settings

    def is_nominal(self):
        return False


class Nominal(SystematicVariation):
    # TODO: Do this with super?
    def __init__(self, direction=None):
        self._name = "Nominal"
        self._direction = direction

    @property
    def name(self):
        name = self._name
        if self._direction:
            name += "_" + self._direction
        return name

    def change_histogram_name(self, h_settings, direction):
        return h_settings

    def shifted_root_objects(self, h_settings):
        return h_settings

    def is_nominal(self):
        return True


class DifferentPipeline(SystematicVariation):
    def __init__(self, name, pipeline, direction):
        super(DifferentPipeline, self).__init__(name, direction)
        self._pipeline = pipeline

    def __str__(self):
        return super(DifferentPipeline, self).__str__()[:-1] + ', pipeline=' + self._pipeline + ')'

    def shifted_root_objects(self, h_settings):
        logger.debug('DifferentPipeline::shifted_root_objects:')
        for index in range(len(h_settings)):
            h_settings[index]["folder"][2] = self._pipeline + self._direction
        logger.debug('    ' + pprint.pformat(h_settings, indent=4))  # ; exit(1)
        return h_settings


from shape_producer.cutstring import Cut


class ReplaceExpressions(SystematicVariation):
    def __init__(self, name, direction, replace_dict):
        super(ReplaceExpressions, self).__init__(name, direction)
        self.replace_dict = replace_dict

    def __str__(self):
        return super(ReplaceExpressions, self).__str__()[:-1] + ', replace_dict=' + self.replace_dict + ')'

    def shifted_root_objects(self, h_settings):
        logger.debug('ReplaceExpressions::shifted_root_objects:')

        # what is the case for having a list of settings?
        returned_h_settings = []
        for index in range(len(h_settings)):
            affected = False
            for old, new in self.replace_dict.iteritems():

                # Replace if found in variable expression
                affected = affected or old in h_settings[index]['variable']._expression
                h_settings[index]['variable']._expression = h_settings[index]['variable']._expression.replace(old, new)

                # Replace if found in weights
                for w in h_settings[index]['weights']._weightstrings:
                    affected = affected or old in w._weightstring
                    w._weightstring = w._weightstring.replace(old, new)

                # Replace if found in cuts
                cuts = h_settings[index]['cuts']
                for c in cuts._cutstrings:
                    affected = affected or old in c._weightstring
                    cut_name = c.name
                    cut_value = c._weightstring.replace(old, new)
                    if old in c.weightstring:
                        try:
                            cuts.remove(cut_name)
                        except:
                            logger.error('Couldn\'t remove cut:', cut_name, ' from list of cuts:', cuts)
                            raise
                        try:
                            cuts.add(Cut(cut_value, cut_name))
                        except:
                            print ('Couldn\'t add the updated cut', cut_name, ' : ', cut_value, 'to the cuts of category')
                            raise
            if affected:
                returned_h_settings.append(h_settings[index])

        return returned_h_settings  # h_settings


class SquareAndRemoveWeight(SystematicVariation):
    def __init__(self, name, weight_name, direction):
        super(SquareAndRemoveWeight, self).__init__(name, direction)
        self._weight_name = weight_name
    def shifted_root_objects(self, h_settings):
        for index in range(len(h_settings)):
            if self._direction == "Up":
                h_settings[index]["weights"] = h_settings[index][
                    "weights"].square(self._weight_name)
            elif self._direction == "Down":
                h_settings[index]["weights"] = h_settings[index][
                    "weights"].remove(self._weight_name)
        return h_settings


class Relabel(SystematicVariation):
    # TODO: Do this with super?
    def __init__(self, name, direction=None):
        self._name = name
        self._direction = direction

    @property
    def name(self):
        name = self._name
        if self._direction:
            name += self._direction
        return name

    def change_histogram_name(self, h_settings, direction):
        return h_settings

    def shifted_root_objects(self, h_settings):
        return h_settings


class ReplaceWeight(SystematicVariation):
    def __init__(self, name, weight_name, new_weight, direction):
        super(ReplaceWeight, self).__init__(name, direction)
        self._weight_name = weight_name
        self._new_weight = new_weight

    def shifted_root_objects(self, h_settings):
        for index in range(len(h_settings)):
            h_settings[index]["weights"] = h_settings[index][
                "weights"].remove(self._weight_name)  #.add(self._new_weight)
            h_settings[index]["weights"].add(self._new_weight)
        return h_settings


class AddWeight(SystematicVariation):
    def __init__(self, name, weight_name, new_weight, direction):
        super(AddWeight, self).__init__(name, direction)
        self._weight_name = weight_name
        self._new_weight = new_weight

    def shifted_root_objects(self, h_settings):
        for index in range(len(h_settings)):
            h_settings[index]["weights"] = h_settings[index]["weights"]
            h_settings[index]["weights"].add(self._new_weight)
        return h_settings
