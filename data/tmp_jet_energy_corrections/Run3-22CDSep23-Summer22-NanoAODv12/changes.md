## Changes: 2025-09-23 (Update JERC JSONs with V3 JECs (Fix Asymptotic Discontinuity in L2Relative txt files))

Merge Request: [!1](https://gitlab.cern.ch/cms-analysis-corrections/JME/Run3-22CDSep23-Summer22-NanoAODv12/-/merge_requests/1)

This MR updates the JERC JSON files for `2022_Summer22` from the jet energy correction version V2 (previously) to V3 (current recommendation). 

The new JEC version (V3) is identical as the previous one (V2) except for a bug fix in the `L2Relative` text file. In extremely rare cases, the V2 `L2Relative` JECs exhibited asymptotic behavior in a few high eta bins and a very tiny pT bin, leading to abnormally large corrected jet pT values. The fraction of affected events is exceedingly small (10e-7%), so the majority of analyses will not observe any noticeable effect. More information about this issue can be found in the corresponding [Gitlab issue](https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/153#note_9486099) and the [presentation slides](https://indico.cern.ch/event/1545816/contributions/6507348/attachments/3066299/5423894/cms-jerc-news_13May2025.pdf#page=2). This issue is fixed in the V3 JECs, which are contained in the JERC JSON files of this MR.

The related PR in JECDatabase is: [#206](https://github.com/cms-jet/JECDatabase/pull/206)
