CMS Analysis Information
========================

In CMS the relevant information about how to handle physics objects in your analysis 
which includes recommendations for triggers, object identification, calibration, 
correction etc. is provided by the ``Physics Object Groups`` (POGs).

The information is collected in either a CERN Twiki page or wiki docs:

.. admonition:: Egamma POG
   
   | Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPOG

.. admonition:: Muon POG
    
   | Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonPOG
   | Wiki docs: https://muon-wiki.docs.cern.ch/

.. admonition:: Tau POG
    
   | Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/Tau
   | Wiki docs: https://tau-wiki.docs.cern.ch/

.. admonition:: JetMET POG
    
   | Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/JetMET
   | Wiki docs: https://cms-jerc.web.cern.ch/

.. admonition:: BTV (B-Tag & Vertexing) POG
    
   | Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/BtagPOG
   | Wiki docs: https://btv-wiki.docs.cern.ch/

.. admonition:: Luminosity POG
    
    Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM

.. admonition:: Generator Group
    
   | Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/GeneratorMain
   | Wiki docs: https://cms-generators.docs.cern.ch/

There is also additional information specifically for the scale factors and how to 
evaluate them using ``correctionlib``. This tool is used in CROWN. 

.. admonition:: json POG

   | CMS Wiki (corrections): https://cms-analysis-corrections.docs.cern.ch/
   | Correctionlib docs: https://cms-nanoaod.github.io/correctionlib/

Further, common analysis recommendations, which is not specifically relevant for
CROWN but for the subsequent analysis steps like CMS plotting style, naming of systematic
uncertainties or statistical tools, can be found here: https://cms-analysis.docs.cern.ch/