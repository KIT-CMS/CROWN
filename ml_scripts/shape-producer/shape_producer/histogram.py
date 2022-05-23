# -*- coding: utf-8 -*-

import ROOT
from array import array
import hashlib
import logging
from . import binning
import os
from XRootD import client
logger = logging.getLogger(__name__)
"""
"""


# Base class for Histogram and Count
class TTreeContent(object):
    def __init__(self,
                 name,
                 inputfiles,
                 folder,
                 cuts,
                 weights,
                 friend_inputfiles_collection=None):  # empty histogram
        self._name = name
        self._inputfiles = [inputfiles] if isinstance(inputfiles,
                                                      str) else inputfiles
        self._friend_inputfiles_collection = [
            ([friend_inputfiles]
             if isinstance(friend_inputfiles, str) else friend_inputfiles)
            for friend_inputfiles in friend_inputfiles_collection
        ]

        self._cuts = cuts
        self._weights = weights
        self._weight_name = "weight_" + self._name  # internal name needed for TDFs
        self._result = None
        self._folder = folder

    def is_present(self):
        return self._result != None

    def files_folders(self):
        return (self._inputfiles, self._folder)

    def apply_cuts_on_dataframe(self, dataframe):
        for cutstring in self._cuts.extract():
            dataframe = dataframe.Filter(cutstring.extract(), cutstring.name)
        return dataframe

    def produce_eventweight(self, dataframe):
        new_dataframe = dataframe.Define(self._weight_name,
                                         self._weights.extract())
        return new_dataframe

    def update(self):
        raise NotImplementedError

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        self.update()

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __hash__(self):
        m = hashlib.md5()
        m.update(self._name)
        for i_file in self._inputfiles:
            m.update(i_file)
        m.update(self._cuts.expand())
        m.update(self._weights.extract())
        m.update(self._folder)
        if isinstance(self, Histogram):
            m.update(self._variable.name)
        return int(m.hexdigest(), 16)

    @property
    def result(self):
        if not self.is_present():
            logger.fatal(
                "The result object of %s called %s has not yet been produced.",
                self, self.name)
            raise Exception
        return self._result


class Histogram(TTreeContent):
    def __init__(self,
                 name,
                 inputfiles,
                 folder,
                 cuts,
                 weights,
                 variable,
                 friend_inputfiles_collection=None):  # empty histogram
        self._variable = variable
        super(Histogram, self).__init__(
            name,
            inputfiles,
            folder,
            cuts,
            weights,
            friend_inputfiles_collection=friend_inputfiles_collection)

    def create_result(self, dataframe=False):
        logger.debug("----->Histogram::create_result::" + self.name)
        if dataframe:
            if not isinstance(self._variable.binning, binning.ConstantBinning):
                logger.fatal("TDataFrames work only with a constant binning.")
                raise Exception
            self._result = dataframe.Histo1D(
                ("", self._cuts.expand() + "*" + self._weights.extract(),
                 self._variable.binning.nbinsx, self._variable.binning.xlow,
                 self._variable.binning.xhigh), self._variable.name,
                self._weight_name)
        else:  # classic way
            # combine files to a single tree using TChain
            logger.debug("------>combine files to a single tree using TChain")
            tree = ROOT.TChain()
            for inputfile in self._inputfiles:
                folder = self._folder
                logger.debug("------>inputfile: " + inputfile + " / " + self._folder)
                check_filepath(inputfile)
                tree.Add(inputfile + "/" + folder)
            # repeat this for friends if applicable
            friend_trees = []
            if self._friend_inputfiles_collection is not None:
                logger.debug(self._friend_inputfiles_collection)
                for i, friend_inputfiles in enumerate(self._friend_inputfiles_collection):
                    logger.debug("------>friend_inputfiles: [" + ', '.join(friend_inputfiles) + ']')
                    friend_tree = ROOT.TChain()
                    for friend_inputfile in friend_inputfiles:
                        folder = self._folder
                        logger.debug("-------->Add friend_inputfile:" + friend_inputfile + " / " + self._folder)
                        check_filepath(friend_inputfile)
                        friend_tree.Add(friend_inputfile + "/" + folder)
                    tree.AddFriend(friend_tree, 'fr%d' % i)
                    friend_trees.append(friend_tree)

            # create unfilled template histogram
            logger.debug("------>create unfilled template histogram")
            hist = ROOT.TH1D(self._name, self._name,
                             self._variable.binning.nbinsx,
                             self._variable.binning.bin_borders)
            # draw histogram and pipe result in the template histogram
            logger.debug("------>draw histogram and pipe result in the template histogram: ")
            logger.debug(self._variable.expression + ">>" + self._name + ' ; \n' + self._cuts.expand() + "*" + self._weights.extract())
            tree.Draw(self._variable.expression + ">>" + self._name,
                      self._cuts.expand() + "*" + self._weights.extract(),
                      "goff")
            # write out result
            logger.debug("------>write out result")
            self._result = ROOT.gDirectory.Get(self._name)
            # reset the chain to close the open files explicitely
            tree.Reset()
            for ft in friend_trees:
                ft.Reset()

        logger.debug("----->Histogram::create_result:: done")
        return self

    def update(self):
        if self.is_present():
            self._result.SetName(self._name)
        else:
            # TODO: What else?
            pass

    def save(self, output_tree, counts):
        self._result.Write("",ROOT.TObject.kOverwrite)

    def has_negative_entries(self):
        if not self.is_present():
            logger.fatal("Histogram %s is not produced.", self.name)
            raise Exception
        for i_bin in range(1, self._result.GetNbinsX() + 1):
            if self._result.GetBinContent(i_bin) < 0.0:
                logger.debug(
                    "Negative value found in histogram %s for bin %d.",
                    self._result, i_bin)
                return True
        return False

    def replace_negative_entries_and_renormalize(self, tolerance):
        if not self.is_present():
            logger.fatal("Histogram %s is not produced.", self.name)
            raise Exception

        # Find negative entries and calculate norm
        norm_all = 0.0
        norm_positive = 0.0
        for i_bin in range(1, self._result.GetNbinsX() + 1):
            this_bin = self._result.GetBinContent(i_bin)
            if this_bin < 0.0:
                self._result.SetBinContent(i_bin, 0.0)
            else:
                norm_positive += this_bin
            norm_all += this_bin

        if norm_all == 0.0 and norm_positive != 0.0:
            logger.fatal(
                "Aborted renormalization because initial normalization is zero, but positive normalization not. . Check histogram %s",
                self.name )
            raise Exception

        if norm_all < 0.0:
            logger.fatal(
                "Aborted renormalization because initial normalization is negative: %f. Check histogram %s ",
                norm_all, self.name)
            raise Exception

        if abs(norm_all - norm_positive) > tolerance * norm_all:
            logger.warning(
                "Renormalization failed because the normalization changed by %f, which is above the tolerance %f. Check histogram %s",
                abs(norm_all - norm_positive), tolerance * norm_all, self.name)
            #raise Exception

        # Renormalize histogram if negative entries are found
        if norm_all != norm_positive:
            if norm_positive == 0.0:
                logger.fatal(
                    "Renormalization failed because all bins have negative entries."
                )
                raise Exception
            for i_bin in range(1, self._result.GetNbinsX() + 1):
                this_bin = self._result.GetBinContent(i_bin)
                self._result.SetBinContent(i_bin,
                                           this_bin * norm_all / norm_positive)


# class to count the (weighted) number of events in a selection
class Count(TTreeContent):
    def __init__(self,
                 name,
                 inputfiles,
                 folder,
                 cuts,
                 weights,
                 friend_inputfiles_collection=None):
        super(Count, self).__init__(
            name,
            inputfiles,
            folder,
            cuts,
            weights,
            friend_inputfiles_collection=friend_inputfiles_collection)
        self._inputfiles = [inputfiles] if isinstance(inputfiles,
                                                      str) else inputfiles

        self._weights = weights
        self._result = False

    def create_result(self, dataframe=False):
        logger.debug("----->Count::create_result::" +
            '\n name: ' + self._name +
            '\n _inputfiles: [' + ', '.join(self._inputfiles) + ']' +
            '\n folder: ' + self._folder
        )
        if dataframe:
            self._result = dataframe.Define("flat", "1").Histo1D(
                "flat", self._weight_name)
        else:  # classic way
            tree = ROOT.TChain()
            for inputfile in self._inputfiles:
                check_filepath(inputfile)
                tree.Add(inputfile + "/" + self._folder)
            # repeat this for friends if applicable
            friend_trees = []
            if self._friend_inputfiles_collection != None:
                for friend_inputfiles in self._friend_inputfiles_collection:
                    friend_tree = ROOT.TChain()
                    for friend_inputfile in friend_inputfiles:
                        check_filepath(friend_inputfile)
                        friend_tree.Add(friend_inputfile + "/" + self._folder)
                    tree.AddFriend(friend_tree)
                    friend_trees.append(friend_tree)

            counthist = ROOT.TH1D(self._name, self._name, 1, 0.0,2.0)
            tree.Draw("1>>" + self._name,
                      self._cuts.expand() + "*" + self._weights.extract(),
                      "goff")
            self._result = counthist.GetBinContent(1)

            # reset the chain and friend trees to close the open files explicitely
            tree.Reset()
            for ft in friend_trees:
                ft.Reset()

        return self

    def save(self, output_tree, counts):
        counts.append(array("f", [self._result]))
        name = self._name
        output_tree.Branch(name, counts[-1], name + "/F")

    # TODO: FIXME: What is this?
    def update(self):
        if not isinstance(self._result, float):
            self._result = self._result.GetBinContent(59)


def check_filepath(inputfile):
    logger.debug("Checking if file is available:{}".format(inputfile))
    if "root://" in inputfile:
        # check if file exists via xrootd
        serverurl = inputfile.split("/")[0] + "//"+inputfile.split("/")[2]
        filepath = "//"+inputfile.strip(serverurl) + ".root"
        myclient = client.FileSystem(serverurl)
        status, info = myclient.stat(filepath)
        if info is None:
            logger.fatal("File not found: {}".format(inputfile))
            raise Exception
    else:
        # check is file is available locally
        if os.path.isfile(inputfile) is False:
            logger.fatal("File not found: {}".format(inputfile))
            raise Exception


# automatic determination of the type
# TODO: Remove me: Too much magic involved
def create_root_object(**kwargs):
    if "variable" in kwargs.keys():
        return Histogram(**kwargs)
    else:
        return Count(**kwargs)


# Helper function to use multiprocessing.Pool with class methods
def root_object_create_result(root_object):
    logger.debug("root_object_create_result")
    return root_object.create_result()


class RootObjects(object):
    def __init__(self, output_filename):
        self._root_objects = []
        self._counts = []
        self._produced = False
        self._output_filename = output_filename
        self._duplicate_root_objects = {}

    def add(self, root_object):
        if self._produced:
            logger.fatal(
                "A produce function has already been called. No more histograms can be added."
            )
            raise Exception
        else:
            if isinstance(root_object, list):
                self._root_objects += root_object
            else:
                self._root_objects.append(root_object)

    def add_unique(self, root_object):
        if self._produced:
            logger.fatal(
                "A produce function has already been called. No more histograms can be added."
            )
            raise Exception
        else:
            if isinstance(root_object, list):
                for r in root_object:
                    add_object = True
                    for ur in self._root_objects:
                        if r == ur:
                            self._duplicate_root_objects.setdefault(
                                ur, []).append(r)
                            add_object = False
                            break
                    if add_object:
                        self._root_objects.append(r)

            else:
                add_object = True
                for ur in self._root_objects:
                    if r == ur:
                        self._duplicate_root_objects.setdefault(ur,
                                                                []).append(r)
                        add_object = False
                        break
                if add_object:
                    self._root_objects.append(r)

    def new_histogram(self, **kwargs):
        self.add(Histogram(**kwargs))

    def new_count(self, **kwargs):
        self.add(Count(**kwargs))

    # get all possible files/folders combinations to determine how many data frames are needed
    def get_combinations(self, *args):
        files_folders = []
        for obj in args:
            for o in obj:
                if not o.files_folders() in files_folders:
                    files_folders.append(o.files_folders())
        return files_folders

    #getter function depending on the histogram name
    def get(self, name):
        for index in range(len(self._root_objects)):
            if self._root_objects[index].name == name:
                return self._root_objects[index]

    def create_output_file(self):
        logger.info("Create output file %s.", self._output_filename)
        self._output_file = ROOT.TFile(self._output_filename, "recreate")
        self._output_tree = ROOT.TTree("output_tree", "output_tree")

    def produce_tdf(self, num_threads):
        self.create_output_file()
        self._produced = True
        # determine how many data frames have to be created; sort by inputfiles and trees
        files_folders = self.get_combinations(self._root_objects, self._counts)

        for files_folder in files_folders:
            # create the dataframe
            common_dataframe = ROOT.Experimental.TDataFrame(
                str(files_folder[1]), str(files_folder[0][0]))
            # loop over the corresponding histograms and create an own dataframe for each histogram -> TODO
            for h in [
                    h for h in self._root_objects
                    if h.files_folders() == files_folder
            ]:
                # find overlapping cut selections -> dummy atm
                special_dataframe = h.apply_cuts_on_dataframe(common_dataframe)
                special_dataframe = h.produce_eventweight(special_dataframe)
                h.create_result(dataframe=special_dataframe)
            # create the histograms
            for h in [
                    h for h in self._root_objects
                    if h.files_folders() == files_folder
            ]:
                h.update()
                h.save(self._output_tree, self._counts)

    def produce_classic(self, num_threads):
        self.create_output_file()
        self._produced = True

        if len(self._root_objects) == 0:
            logger.info(
                "produce_classic : %u shapes in classic mode and %u processes. -> nothing to produce, returning",
                len(self._root_objects), num_threads)
            return self

        logger.info(
            "produce_classic : Start to produce %u shapes in classic mode and %u processes.",
            len(self._root_objects), num_threads)
        if num_threads == 1:
            for ro in self._root_objects:
                ro.create_result()
        else:
            from multiprocessing import Pool
            pool = Pool(processes=num_threads)
            root_objects_new = pool.map(root_object_create_result,
                                        [ro for ro in self._root_objects])
            pool.close()
            pool.join()

            # Because the new objects have different addresses in memory,
            # the result objects have to be copied.
            # Otherwise, the systematic's root_object has no result associated.
            for i_ro in range(len(root_objects_new)):
                self._root_objects[i_ro]._result = root_objects_new[
                    i_ro]._result

        for h in self._root_objects:  # write sequentially to prevent race conditions
            h.save(self._output_tree, self._counts)
        return self

    def set_duplicates(self):
        logger.debug(
            "set_duplicates : Setting duplicates to the corresponding produced ROOT objects",
        )
        for ro in self._root_objects:
            if ro in self._duplicate_root_objects:
                logger.debug(
                    "Setting following duplicates for produced object %s with address %s",
                    ro.name, ro)
                for dup_ro in self._duplicate_root_objects[ro]:
                    logger.debug("duplicate %s with address %s", dup_ro.name,
                                 dup_ro)
                    dup_ro._result = ro.result
                    # self._root_objects[self._root_objects.index(dup_ro)]._result = ro.result
                logger.debug("--------------------------------")

    def save(self):
        if not self._produced:
            logger.fatal("No produce method has been called for %s.", self)
            raise Exception
        self._output_tree.Fill()
        self._output_tree.Write("",ROOT.TObject.kOverwrite)
        self._output_file.Close()
        logger.info("Output file closed: %s.", self._output_file)
