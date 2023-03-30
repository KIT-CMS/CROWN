KingMaker
===========

KingMaker is a workflow management for producing ntuples with the CROWN framework. The workflow management is based on  law (github.com/riga/law), which is using  luigi (https://github.com/spotify/luigi) as backend.

Setup
-----

.. code-block:: bash

    git clone --recursive git@github.com:KIT-CMS/KingMaker.git
    cd KingMaker
    source setup.sh KingMaker

This should install all required packages and setup the environment. In addition, a ``luigid`` scheduler is started, if not already running. The required port is set to the ```LUIGIPORT``` environment variable.

Management of samples
---------------------



Submission of ntuples
---------------------


Submission of friend trees
--------------------------