# Objects

## Muons

In nanoAOD, a base selection is already applied on muons (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Muons):

```
 SlimmedMuons passing the following selection are stored in nanoAOD:
    - pT(mu) > 3 GeV
    - pass one of the muonID: 'CutBasedIdLoose' || 'SoftCutBasedId' || 'SoftMvaId' || 'CutBasedIdGlobalHighPt' || 'CutBasedIdTrkHighPt'
```

In insertframeworkname, the following function are available to select muons:

## Electrons

In nanoAOD, a base selection is already applied on Electrons (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Muons). All `SlimmedElectrons` with `pT(e) > 5 GeV` are stored in nanoAOD

## Taus

All Taus, that pass a cut with `pT(tau) > 18 GeV` and at least one tau disciminator are stored in nanoAOD.