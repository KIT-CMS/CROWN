Changelog
==========

May 2024 - Version 0.4.0

* Switch to ROOT 6.30, supporting now RHEL8 and RHEL9.
* Introduced support for ML inference via `OnnxRuntime`. A generic producer is avialable at this link: https://github.com/KIT-CMS/CROWN/blob/main/include/ml.hxx
* Introduced a CorrectionManager, that is responsible for loading correction files and sharing them among the different producers. This allows to load the corrections only once and share them among the different producers, resulting in a signifiant speedup of the initial loading time. In the course of this implementation, Many functions now have a deprecated version, that does not use the `CorrectionManager`. The old functions will be removed in the next release. A more detailed description can be found in the :ref:`The Correction Manager` page.

Sept. 2023 - Version 0.3.0

* Switched to ROOT 6.28 via LCG 104, resulting in about 20% faster processing times.
* Added support for the generation of friend trees with additional friends as input. For more details, check :ref:`FriendTree Generation`.
* Added option to compile the CROWNlib only, allowing to reuse the same libary for multiple CROWN executables.

Feb. 2023

* Added support for the generation of friend trees. For more details, check :ref:`FriendTree Generation`.
* Added documentation on ntuple and friend production via KingMaker. For more details, check :ref:`KingMaker`.

Jan. 2023

* Added Quantities <-> Shifts mapping to the output files to allow an easier Postprocessing. For more details, check :ref:`Quantity mapping`.
