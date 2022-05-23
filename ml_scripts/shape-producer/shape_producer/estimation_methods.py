# -*- coding: utf-8 -*-

import os

from .histogram import *
from .cutstring import *
from .systematics import *
from .systematic_variations import *

import logging
logger = logging.getLogger(__name__)


class EstimationMethod(object):
    def __init__(self,
                 name,
                 folder,
                 era,
                 directory,
                 channel,
                 mc_campaign,
                 friend_directory=None,
                 get_triggerweight_for_channel=None,
                 get_singlelepton_triggerweight_for_channel=None,
                 get_tauByIsoIdWeight_for_channel=None,
                 get_eleHLTZvtxWeight_for_channel=None,
                 heavy_mass=None,
                 light_mass=None):
        self._directory = directory
        self._folder = folder
        self._name = name
        self._mc_campaign = mc_campaign
        self._channel = channel
        self._era = era
        import six
        self._friend_directories = [friend_directory] if isinstance(
            friend_directory, six.string_types) else friend_directory

        self.get_triggerweight_for_channel = get_triggerweight_for_channel
        self.get_singlelepton_triggerweight_for_channel = get_singlelepton_triggerweight_for_channel
        self.get_tauByIsoIdWeight_for_channel = get_tauByIsoIdWeight_for_channel
        self.get_eleHLTZvtxWeight_for_channel = get_eleHLTZvtxWeight_for_channel
        self.heavy_mass = heavy_mass
        self.light_mass = light_mass
        for i in ['get_triggerweight_for_channel', 'get_singlelepton_triggerweight_for_channel', 'get_tauByIsoIdWeight_for_channel', 'get_eleHLTZvtxWeight_for_channel']:
            if getattr(self, i) is None:
                setattr(self, i, lambda x=None: (_ for _ in ()).throw(
                    Exception('Call of undefined method "%s" with arguments: %s' % (i, str(x))))
                )

    def get_path(self, systematic, folder):
        logger.debug('-------->EstimationMethod::get_path: ' + systematic.category.channel.name + "_" + folder + "/ntuple")
        return systematic.category.channel.name + "_" + folder + "/ntuple"

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def era(self):
        return self._era

    @property
    def channel(self):
        return self._channel

    def get_weights(self):
        return Weights(Weight("1.0", "constant"))

    def get_cuts(self):
        return Cuts()

    # function parsing the datasets helper to return the files
    # overwrite this function
    # TODO: Error handling for missing files!
    def get_files(self):
        raise NotImplementedError

    def get_friend_files(self):
        return [[
            filename.replace(self._directory, friend_directory)
            for filename in self.get_files()
        ] for friend_directory in self._friend_directories]

    def artus_file_names(self, files):
        return [os.path.join(self._directory, f, "%s.root" % f) for f in files]

    # wrapper function for the Histogram creation performing the systematic shifts
    def apply_systematic_variations(self, systematic, settings):
        return systematic.variation.shifted_root_objects(settings)

    def define_root_objects(self, systematic):
        logger.debug('------>EstimationMethod::define_root_objects: folder : ' + self._folder)
        histogram_settings = []
        histogram_settings.append({
            "name":
            systematic.name,
            "inputfiles":
            self.get_files,
            "folder": [self.get_path, systematic, self._folder],
            "cuts":
            systematic.category.cuts +
            self.get_cuts(),  # TODO: get_cuts() is correct here?
            "weights":
            systematic._additionalWeights + self.get_weights()
        })
        if systematic.category.variable != None:
            histogram_settings[-1]["variable"] = systematic.category.variable
        if self._friend_directories != None:
            histogram_settings[-1][
                "friend_inputfiles_collection"] = self.get_friend_files
        return histogram_settings

    # TODO: Make this less magic
    def create_root_objects(self, systematic):
        logger.debug('---->EstimationMethod::create_root_objects: systematic.name : ' + systematic.name)
        root_object_settings = self.define_root_objects(systematic)

        # execute the underlying functions
        root_object_settings = self.apply_systematic_variations(
            systematic, root_object_settings)
        for setting in root_object_settings:
            for key, value in setting.iteritems():
                if logger.getEffectiveLevel() == 10: print ('k,v:', key, value)
                if isinstance(value, list):
                    setting[key] = value[0](*value[1:])  # returned: <pipeline>/ntuple
                elif callable(value):
                    setting[key] = value()
                    if logger.getEffectiveLevel() == 10: print (setting[key])

        root_objects = []
        for setting in root_object_settings:
            root_objects.append(create_root_object(**setting))
        return root_objects

    # doing nothing, shape is exactly the histogram as default
    def do_estimation(self, systematic):
        if len(systematic.root_objects) > 1:
            logger.fatal(
                "There are %d histograms associated to the systematic with name %s, but not exactly 1.",
                len(systematic.root_objects), systematic.name)
            raise Exception
        elif len(systematic.root_objects) == 0:
            logger.warning(
                "There are 0 histograms associated to the systematic with name %s -> skipped",
                systematic.name)
            return None
        return systematic.root_objects[0]


class SStoOSEstimationMethod(EstimationMethod):
    def __init__(self,
                 name,
                 folder,
                 era,
                 directory,
                 channel,
                 bg_processes,
                 data_process,
                 friend_directory=None,
                 extrapolation_factor=1.0,
                 qcd_weight = Weight("1.0","qcd_weight"),
                 get_triggerweight_for_channel=None,
                 get_singlelepton_triggerweight_for_channel=None,
                 get_tauByIsoIdWeight_for_channel=None,
                 get_eleHLTZvtxWeight_for_channel=None,):
        super(SStoOSEstimationMethod, self).__init__(
            name=name,
            folder=folder,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,)
        self._bg_processes = [copy.deepcopy(p) for p in bg_processes]
        self._data_process = copy.deepcopy(data_process)
        self._extrapolation_factor = extrapolation_factor
        self._qcd_weight = qcd_weight

    def create_root_objects(self, systematic):
        ss_category = copy.deepcopy(systematic.category)
        ss_category.cuts.get("os").name = "ss"
        ss_category.cuts.get("ss").invert()
        ss_category.name = ss_category._name + "_ss"
        root_objects = []
        systematic._qcd_systematics = []
        for process in [self._data_process] + self._bg_processes:
            s = Systematic(
                category=ss_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125,
                additionalWeights=Weights(self._qcd_weight))
            systematic._qcd_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_qcd_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _qcd_systematics needed for QCD estimation.",
                systematic.name)
            raise Exception

        # Create shapes
        for s in systematic._qcd_systematics:
            s.do_estimation()

        # Data shape
        shape = copy.deepcopy(systematic._qcd_systematics[0].shape)

        # Subtract MC shapes from data shape
        for s in systematic._qcd_systematics[1:]:
            if not isinstance(shape.result,float):
                shape.result.Add(s.shape.result, -1.0)
            else:
                shape._result -= s.shape.result

        final_shape = copy.deepcopy(shape)
        # Saving shape in ss region
        if not isinstance(final_shape.result,float):
            shape._result = final_shape.result.Clone()
        else:
            shape._result = final_shape.result

        shape.name = systematic.name.replace(
            systematic.category._name,
            systematic._qcd_systematics[0].category._name)

        if not isinstance(final_shape.result,float):
            shape._result.Write()
            if shape._result.Integral() > 0.0:
                # Scale shape with extrapolation factor
                final_shape.result.Scale(self._extrapolation_factor)
            else:
                logger.warning("No data in same-sign region for systematic %s. Setting extrapolation factor to 0.0", systematic.name)
                final_shape.result.Scale(0.0)
        else:
            if shape._result > 0.0:
                # Scale shape with extrapolation factor
                final_shape._result *= self._extrapolation_factor
            else:
                logger.warning("No data in same-sign region for systematic %s. Setting extrapolation factor to 0.0", systematic.name)
                final_shape._result *= 0.0

        # Rename root object accordingly
        final_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        if not isinstance(final_shape.result,float):
            final_shape.replace_negative_entries_and_renormalize(tolerance=100.05)

        return final_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class ABCDEstimationMethod(EstimationMethod):
    def __init__(
            self,
            name,
            folder,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            AC_cut_names,
            BD_cuts,
            AB_cut_names,
            CD_cuts,
            friend_directory=None,
            get_triggerweight_for_channel=None,
            get_singlelepton_triggerweight_for_channel=None,
            get_tauByIsoIdWeight_for_channel=None,
            get_eleHLTZvtxWeight_for_channel=None,
    ):  #last four arguments correspond to 1. list of names of cuts to be removed in order to include the sideband for the shape derivation 2. list of cuts to be applied to restrict on that sideband and 3.,4. accordingly for the extrapolation factor sideband
        super(ABCDEstimationMethod, self).__init__(
            name=name,
            folder=folder,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,)
        self._bg_processes = [copy.deepcopy(p) for p in bg_processes]
        self._data_process = copy.deepcopy(data_process)
        self._AC_cut_names = AC_cut_names
        self._AB_cut_names = AB_cut_names
        self._BD_cuts = BD_cuts
        self._CD_cuts = CD_cuts

    def create_root_objects(self, systematic):
        # prepare sideband categories
        B_category = copy.deepcopy(systematic.category)
        B_category._name += "_B"
        C_category = copy.deepcopy(systematic.category)
        C_category._name += "_C"
        D_category = copy.deepcopy(systematic.category)
        D_category._name += "_D"
        # C and D are to be created as counts in order to determine the extrapolation factor
        C_category._variable = None
        D_category._variable = None
        # check whether cuts limiting signal region A are applied and remove them from B,C,D
        for cut_name in self._AB_cut_names:
            if cut_name in systematic.category.cuts.names:
                C_category.cuts.remove(cut_name)
                D_category.cuts.remove(cut_name)
            else:
                logger.fatal(
                    "The name %s is not part of the selection cuts of the signal region.",
                    cut_name)
                raise KeyError
        for cut_name in self._AC_cut_names:
            if cut_name in systematic.category.cuts.names:
                B_category.cuts.remove(cut_name)
                D_category.cuts.remove(cut_name)
            else:
                logger.fatal(
                    "The name %s is not part of the selection cuts of the signal region.",
                    cut_name)
                raise KeyError
        # apply other cuts instead
        for cut in self._BD_cuts:
            B_category.cuts.add(cut)
            D_category.cuts.add(cut)
        for cut in self._CD_cuts:
            C_category.cuts.add(cut)
            D_category.cuts.add(cut)

        root_objects = []
        systematic._ABCD_systematics = []
        for process in [self._data_process] + self._bg_processes:
            for category in [B_category, C_category, D_category]:
                s = Systematic(
                    category=category,
                    process=process,
                    analysis=systematic.analysis,
                    era=self.era,
                    variation=systematic.variation,
                    mass=125)
                if process == self._data_process:
                    direction = s.variation._direction
                    s.variation = Nominal()
                    s.variation._direction = direction
                systematic._ABCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_ABCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _ABCD_systematics needed for ABCD estimation.",
                systematic.name)
            raise Exception

        # Create shapes
        B_shapes = {}
        C_shapes = {}
        D_shapes = {}
        for s in systematic._ABCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("_B"):
                B_shapes[s.process.name] = s.shape
            elif s.category.name.endswith("_C"):
                C_shapes[s.process.name] = s.shape
            elif s.category.name.endswith("_D"):
                D_shapes[s.process.name] = s.shape

        # Determine extrapolation factor
        data_C = C_shapes[self._data_process.name].result
        C_yield = C_shapes.pop(self._data_process.name).result - sum(
            [s.result for s in C_shapes.values()])
        data_D = D_shapes[self._data_process.name].result
        D_yield = D_shapes.pop(self._data_process.name).result - sum(
            [s.result for s in D_shapes.values()])

        if data_C == 0.0 or data_D == 0.0:
            logger.warning("No data in C or D region in ABCD method for systematic %s. Setting extrapolation factor to 0.0", systematic.name)
            extrapolation_factor = 0.0
        elif not D_yield > 0.0:
            logger.warning("D_yield in ABCD method for systematic %s is %f.",
                         systematic.name, D_yield)
            extrapolation_factor = 0.0
        else:
            extrapolation_factor = C_yield / D_yield
        logger.debug("D to C extrapolation factor: %s",
                     str(extrapolation_factor))

        # Derive final shape
        derived_shape = copy.deepcopy(B_shapes.pop(self._data_process.name))
        for s in B_shapes.values():
            if not isinstance(derived_shape.result,float):
                derived_shape.result.Add(s.result, -1.0)
            else:
                derived_shape._result -= s.result
        if not isinstance(derived_shape.result,float):
            derived_shape.result.Scale(extrapolation_factor)
        else:
            derived_shape.result *= extrapolation_factor

        # Rename root object accordingly
        derived_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        derived_shape.replace_negative_entries_and_renormalize(
            tolerance=100.05)

        return derived_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class AddHistogramEstimationMethod(EstimationMethod):
    def __init__(self, name, folder, era, directory, channel, add_processes,
                 add_weights,
             get_triggerweight_for_channel=None,
             get_singlelepton_triggerweight_for_channel=None,
             get_tauByIsoIdWeight_for_channel=None,
             get_eleHLTZvtxWeight_for_channel=None,):
        super(AddHistogramEstimationMethod, self).__init__(
            name=name,
            folder=folder,
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,)
        self._add_processes = [copy.deepcopy(p) for p in add_processes]
        self._add_weights = copy.deepcopy(add_weights)

    def create_root_objects(self, systematic):

        root_objects = []
        systematic._add_systematics = []
        for process in self._add_processes:
            s = Systematic(
                category=systematic.category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._add_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):

        # Create shapes
        for s in systematic._add_systematics:
            s.do_estimation()

        # First shape
        final_shape = copy.deepcopy(systematic._add_systematics[0].shape)

        # Add/subtract additional shapes from first shape
        for s in systematic._add_systematics[1:]:
            if not isinstance(final_shape.result,float):
                final_shape.result.Add(s.shape.result, self._add_weights[-1])
            else:
                final_shape.result += s.shape.result * self._add_weights[-1]

        final_shape.name = systematic.name

        return final_shape

    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class SumUpEstimationMethod(EstimationMethod):
    def __init__(self,
                 name,
                 folder,
                 era,
                 directory,
                 channel,
                 processes,
                 factors=None,
                 friend_directory=None,
                 get_triggerweight_for_channel=None,
                 get_singlelepton_triggerweight_for_channel=None,
                 get_tauByIsoIdWeight_for_channel=None,
                 get_eleHLTZvtxWeight_for_channel=None,):
        super(SumUpEstimationMethod, self).__init__(
            name=name,
            folder=folder,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            channel=channel,
            mc_campaign=None)
        self._processes = [copy.deepcopy(p) for p in processes]
        if factors != None:
            if len(processes) != len(factors):
                logger.fatal(
                    "In SumUpEstimationMethod, number of factors must match number of processes!"
                )
                raise Exception
            self._factors = factors
        else:
            self._factors = [1.0] * len(processes)

    def create_root_objects(self, systematic):
        root_objects = []
        systematic._sumUp_systematics = []
        sum_category = copy.deepcopy(systematic.category)
        sum_category._name += "_sum" + self._name
        for process in self._processes:
            s = Systematic(
                category=sum_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._sumUp_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_sumUp_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _sumUp_systematics needed for summation.",
                systematic.name)
            raise Exception

        # Create shapes
        shapes = []
        for s in systematic._sumUp_systematics:
            s.do_estimation()
            shapes.append(s.shape)
        derived_shape = None
        for shape, factor in zip(shapes, self._factors):
            if isinstance(shape, Histogram):
                if derived_shape == None:
                    derived_shape = shape
                    derived_shape.result.Scale(factor)
                else:
                    derived_shape.result.Add(shape.result, factor)
            elif isinstance(shape, Count):
                if derived_shape == None:
                    derived_shape = shape
                    derived_shape._result *= factor
                else:
                    derived_shape._result += shape.result * factor
            else:
                logger.fatal(
                    "SumUpEstimationMethod expects Histogram or Count")
                raise Exception

        # Rename root object accordingly
        derived_shape.name = systematic.name

        return derived_shape


class NewFakeEstimationMethodLT(EstimationMethod):
    def __init__(self,
                 name,
                 folder,
                 era,
                 directory,
                 channel,
                 nofake_processes,
                 data_process,
                 aisoCut,
                 fakeWeightstring,
                 friend_directory=None,
                 get_triggerweight_for_channel=None,
                 get_singlelepton_triggerweight_for_channel=None,
                 get_tauByIsoIdWeight_for_channel=None,
                 get_eleHLTZvtxWeight_for_channel=None,):
        super(NewFakeEstimationMethodLT, self).__init__(
            name=name,
            folder=folder,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            channel=channel,
            mc_campaign=None)
        self._nofake_processes = [copy.deepcopy(p) for p in nofake_processes]
        self._data_process = copy.deepcopy(data_process)
        self._cutsToRemove = ["tau_iso"]
        self._aisoCut = aisoCut
        self._FFweight = Weight(fakeWeightstring, "fake_factor")

    def create_root_objects(self, systematic):
        FF_category = copy.deepcopy(systematic.category)
        for cut in self._cutsToRemove:
            FF_category.cuts.remove(cut)
        FF_category.cuts.add(self._aisoCut)
        FF_category.name = FF_category._name + "_FF"

        '''sdata = Systematic(
            category=FF_category,
            process=self._data_process,
            analysis=systematic.analysis,
            era=self.era,
            variation=systematic.variation if "_ff_" in systematic.variation.name else Nominal(), # only apply ff systematics to the data shape
            mass=125,
            additionalWeights=Weights(self._FFweight))
        systematic._ff_systematics = [sdata]
        sdata.create_root_objects()
        root_objects = sdata.root_objects'''
        systematic._ff_systematics = []
        root_objects = []
        for process in [self._data_process] + self._nofake_processes:
            process.estimation_method._friend_directories = self._friend_directories
            s = Systematic(
                category=FF_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation if ("_ff_" in systematic.variation.name or not "data" in process.name) else Nominal(),
                mass=125,
                additionalWeights=Weights(self._FFweight))
            systematic._ff_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_ff_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _ff_systematics needed for FF estimation.",
                systematic.name)
            raise Exception

        # Create shapes
        for s in systematic._ff_systematics:
            s.do_estimation()

        # Data shape
        shape = copy.deepcopy(systematic._ff_systematics[0].shape)

        # Subtract MC shapes from data shape
        for s in systematic._ff_systematics[1:]:
            if not isinstance(shape.result,float):
                shape.result.Add(s.shape.result, -1.0)
            else:
                shape._result -= s.shape.result

        # Rename root object accordingly
        shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        #shape.replace_negative_entries_and_renormalize(tolerance=100.05)

        return shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class NewFakeEstimationMethodTT(NewFakeEstimationMethodLT):
    def __init__(self,
                 name,
                 folder,
                 era,
                 directory,
                 channel,
                 nofake_processes,
                 data_process,
                 aisoCut,
                 fakeWeightstring,
                 friend_directory=None,
                 get_triggerweight_for_channel=None,
                 get_singlelepton_triggerweight_for_channel=None,
                 get_tauByIsoIdWeight_for_channel=None,
                 get_eleHLTZvtxWeight_for_channel=None,):
        super(NewFakeEstimationMethodLT, self).__init__(
            name=name,
            folder=folder,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            channel=channel,
            mc_campaign=None)
        self._nofake_processes = [copy.deepcopy(p) for p in nofake_processes]
        self._data_process = copy.deepcopy(data_process)
        self._cutsToRemove = ["tau_1_iso", "tau_2_iso"]
        self._aisoCut = aisoCut
        self._FFweight = Weight(fakeWeightstring, "fake_factor")
