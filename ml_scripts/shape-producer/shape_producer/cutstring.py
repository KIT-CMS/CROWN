# -*- coding: utf-8 -*-

import sys
import logging
import copy
logger = logging.getLogger(__name__)
"""
"""

# Cut -> base class for sorting them
# Cuts -> holder for selection steps
# Weight -> Weight, constant
# Weights -> holder for weight expression

supported_operators = ['<', '>', '&&', '||', '==', '!=', '<=', '>=']
inverted_operators = ['>=', '<=', '||', '&&', '!=', '==', '>', '<']


# Base class all others inherit from.
# Usage: For any weight string that does not follow one of the special definitions listed below
class Weight(object):
    def __init__(self, weightstring, name=False):
        self._weightstring = weightstring
        if name == False:
            logger.fatal(
                "No appropriate name has been assigned to weight with value %s. Please use an explicit name for this weight string.",
                str(weightstring))
            raise ValueError
        self._name = name

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'Weight("%s" : "%s")' % (self._name, self._weightstring)

    @property
    def name(self):
        return self._name

    def extract(self):
        return self.embrace(self._weightstring)

    def embrace(self, c):
        return "(" + c + ")"


# Class for a constant number, e.g. '0.9' or a constant variable e.g. "generatorWeight"
class Constant(Weight):
    def __init__(self, weightstring, name=False):
        self.is_float = False
        self._weightstring = weightstring
        try:
            float(weightstring)
            is_float = True
            self._name = name
        except:
            self._name = weightstring if name == False else name

        if name == False:
            logger.fatal(
                "No appropriate name has been assigned to weight with value %s. Please use an explicit name for this weight string.",
                str(weightstring))
            raise ValueError

    # TODO: Remove the replace("1.0/") feature if not used, too hacky!
    def invert(self):
        if is_float:
            self._weightstring = str(1.0 / float(self._weightstring))
        else:
            if self._weightstring.startswith("1.0/"):
                self._weightstring = self._weightstring.replace("1.0/", "")
            else:
                self._weightstring = "1.0/" + self._weightstring


# holder class for weight/cutstring objects
class Weights(object):
    def __init__(self, *args):
        self._weightstrings = []
        if args != False:
            for w in args:
                self.add(w)

    def __str__(self):
        print_str = super(type(self), self).__repr__().split(' object at')[0][1:] + '\n'  # ': [' + str(hex(id(self))) + ']\n'

        if len(self._weightstrings) == 0:
            print_str = 'Weights: empty'
        else:
            for i in self._weightstrings:
                print_str += '\t' + str(i).replace('\n', '\n\t') + '\n'
        return print_str

    def __repr__(self):
        print_str = super(type(self), self).__repr__().split(' object at')[0][1:] + ': [' + str(hex(id(self))) + ']\n'

        if len(self._weightstrings) == 0:
            print_str = super(type(self), self).__repr__().split(' object at')[0][1:] + str(hex(id(self))) + ': empty'
        else:
            for i in self._weightstrings:
                print_str += '\t' + str(i).replace('\n', '\n\t') + '\n'
        return print_str

    def add(self, weightstring):
        if (issubclass(type(weightstring), Weight)):
            if weightstring.name in self.names:
                logger.fatal(
                    "Not possible to add the weightstring %s since its name is not unique.",
                    weightstring)
                raise LookupError
            else:
                self._weightstrings.append(weightstring)

    def __add__(self, other):
        new = copy.deepcopy(self)
        for c in other._weightstrings:
            new.add(c)
        return new

    def extract(self):
        if len(self._weightstrings) > 0:
            full_weightstring = "*".join(
                [c.extract() for c in self._weightstrings])
            return full_weightstring
        else:
            logger.fatal(
                "Can't extract full weightstring because no weights are given."
            )
            raise Exception

    @property
    def names(self):
        return [w.name for w in self._weightstrings]

    def get(self, name):
        for w in self._weightstrings:
            if w.name == name:
                return w
        logger.fatal("The name %s is not part of this Weights object.", name)
        raise KeyError

    # TODO: Implement this using pop or delete
    def remove(self, name):
        if name in self.names:
            self._weightstrings = [
                w for w in self._weightstrings if not w.name == name
            ]
        else:
            logger.fatal(
                "Error while trying to remove weightstring with key %s", name)
            raise KeyError
        return self

    def square(self, name):
        new = copy.deepcopy(self.get(name))
        new._name = new._name + "_squared"
        self._weightstrings.append(new)
        return self


# Class for a simple cut expression e.g. 'pt_1>22'
class Cut(object):
    def __init__(self, cutstring, name=False):
        self._varleft = None
        self._varright = None
        self._operator = None
        self._weightstring = cutstring
        if name != False:
            self._name = name
        else:
            self._name = filter(str.isalnum,
                                cutstring.replace(">", "Gt").replace(
                                    "<", "St"))

        # test if simple, parseable cutstring
        # TODO: Corner cases? They are all checked?
        # TODO: error handling?
        try:
            logger.debug('\t cutstring:' + cutstring)
            operators = [s for s in supported_operators if s in cutstring]
            self._operator = operators[0]
            tmpcutstring = cutstring.split(self._operator)
            if len(tmpcutstring) == 2 and len(operators) == 1:
                self._varleft = tmpcutstring[0]
                try:
                    self._varright = int(tmpcutstring[1])
                except ValueError:
                    self._varright = float(tmpcutstring[1])
            self.update_weightstring()
        except:
            logger.fatal(
                "Failed to compose cut from string \'{}\'.".format(cutstring))
            raise Exception

    def __str__(self):
        return '{%s : %s}' % (self._name, self._weightstring)

    @property
    def weightstring(self):
        return self._weightstring

    def invert(self):
        # TODO: error handling?
        self._operator = inverted_operators[supported_operators.index(
            self._operator)]
        self.update_weightstring()
        return self  # TODO: Remove this

    def update_weightstring(self):
        if (self._varleft != None) and (self._operator !=
                                        None) and (self._varright != None):
            self._weightstring = "".join(
                [self._varleft, self._operator,
                 str(self._varright)])
        return self  # TODO: remove this

    @property
    def value(self):
        return self._varright

    @value.setter
    def value(self, value):
        self._varright = float(value)
        self.update_weightstring()

    @property
    def variable(self):
        return self._varleft

    @variable.setter
    def variable(self, variable):
        self._varleft = variable
        self.update_weightstring()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not value.isalnum():
            logger.fatal("Given name is not alphanumeric.")
            raise Exception
        self._name = value

    def extract(self):
        return self.embrace(self._weightstring)

    def embrace(self, c):
        return "(" + c + ")"


# holder class for cutstring objects
class Cuts(object):
    def __init__(self, *args):
        self._cutstrings = []
        if args != False:
            for w in args:
                self.add(w)

    def __str__(self, indent=1):
        s = '\n'
        for cut in self._cutstrings:
            s += '\t' * indent + cut.name + ' : ' + cut._weightstring + '\n'
        return s

    # TODO: Remove this magic due to consistence
    def __add__(self, other):
        new_cuts = Cuts()
        for c in self._cutstrings:
            new_cuts.add(c)
        for c in other._cutstrings:
            new_cuts.add(c)
        return new_cuts

    def add(self, cutstring):
        if (issubclass(type(cutstring), Cut)) or isinstance(
                cutstring,
                Cut):  # TODO: Make this easier, remove one if not needed
            if cutstring.name in self.names:
                if (self.get(cutstring.name).weightstring).replace(" ", "") != (cutstring.weightstring).replace(" ", ""):
                    logger.fatal(
                        "Not possible to add the cutstring %s since its name is not unique.",
                        cutstring)
                    raise LookupError
            else:
                self._cutstrings.append(cutstring)
        else:
            logger.fatal("Cut %s is not of type `Cut`.", cutstring)
            raise TypeError

    def extract(self):
        return self._cutstrings

    def expand(self):
        if len(self._cutstrings) > 0:
            return "*".join([c.extract() for c in self._cutstrings])
        else:
            return "(1.0)"

    @property
    def names(self):
        return [w.name for w in self._cutstrings]

    def get(self, name):
        for w in self._cutstrings:
            if w.name == name:
                return w
        logger.fatal("The name %s is not part of this Cuts object.", name)
        raise KeyError

    # TODO: Use pop or remove
    def remove(self, name):
        if name in self.names:
            self._cutstrings = [
                w for w in self._cutstrings if not w.name == name
            ]
        return self  # TODO: Remove this
