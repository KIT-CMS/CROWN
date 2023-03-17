FriendTree Generation
===========================

CROWN can be used, to generate FriendTrees based on a CROWN ntuple. The concept of FriendTrees is explained here: https://root.cern/manual/trees/#widening-a-ttree-through-friends. They allow to extend an existing ntuple with new quantities. Common use cases are new high level variables like neural network outputs or additional correction factors. 

A FriendTree is generated using a FriendTreeConfiguration. Such a configuration has some major differences, compared to a regular configuration:

1. The input file is a CROWN ntuple, not a ROOT file.
2. Only one scope per user is allowed.
3. No global scope is required
4. No optimization of the configuration is donw (for now)
5. The available inputs have to be specified. The available inputs can be provided by using a CROWN ntuple as input, or a json file. The ntuple can be used for debugging proposes, when running a production, it is recommended to use a json file. The basic structure this quantities map is listed below. Such a json can then be used for multiple eras, sampletypes and scopes.

.. warning::
    One limitation of the current FriendTree generation is that only one ntuple can be used as input, additional friends as input are not supported yet.

.. code-block:: json

    {
        "era_1": {
            "sampletype_1": {
                "scope_1": {
                    "shift_1": [
                        "id_tau_vsJet_VLoose_2",
                        "id_tau_vsJet_VVLoose_2",
                        "jphi_1",
                        "q_2",
                        "id_wgt_tau_vsEle_VLoose_2",
                        "id_tau_vsEle_VTight_2",
                        "extramuon_veto",
                        "id_tau_vsJet_Tight_2",
                        "mass_1",
                        "puweight"
                    ],
                },
            },
        },
    }



The recommended way of producing FriendTrees is to use a workflow tool, that manages the submission of jobs, generation of tarballs and organizing the output. One possible workflow tool choice is KingMaker (https://github.com/KIT-CMS/KingMaker). A more detailed description of the KingMaker workflow can be found in :ref:`KingMaker`.

Writing a FriendTreeConfiguration
---------------------------------

The basic structure of a FriendTreeConfiguration is identical to a regular configuration. When creating a new FriendTree executable, an additional argument has to be provided:

* ``DQUANTITIESMAP`` - The path to the quantities map json file or the crown ntuple root file.

All other parameters are identical to the regular configuration. Setting up producers, outputs and new systematic shifts works the same way as before. The configuration has to be of type ``FriendTreeConfiguration``. When multiple scopes are provided, one executable per scope is generated.