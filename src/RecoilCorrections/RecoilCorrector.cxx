#include "../../include/RecoilCorrections/RecoilCorrector.hxx"
#include "../../include/utility/Logger.hxx"

RecoilCorrector::RecoilCorrector(std::string filepath) {
    fileName = filepath;
    TFile *file = new TFile(fileName, "READ");
    if (file->IsZombie()) {
        Logger::get("RecoilCorrector")
            ->debug("file {} is not found...   quitting ", fileName);
        exit(-1);
    }

    TH1D *projH = (TH1D *)file->Get("projH");
    if (projH == NULL) {
        Logger::get("RecoilCorrector")
            ->debug("File should contain histogram with the name projH ");
        Logger::get("RecoilCorrector")
            ->debug("Check content of the file {}", fileName);
        exit(-1);
    }

    TString firstBinStr = projH->GetXaxis()->GetBinLabel(1);
    TString secondBinStr = projH->GetXaxis()->GetBinLabel(2);

    TString paralZStr = firstBinStr;
    TString perpZStr = secondBinStr;
    if (firstBinStr.Contains("Perp")) {
        paralZStr = secondBinStr;
        perpZStr = firstBinStr;
    }
    Logger::get("RecoilCorrector")
        ->debug("Parallel component      (U1) : {}", paralZStr);
    Logger::get("RecoilCorrector")
        ->debug("Perpendicular component (U2) : {}", perpZStr);

    TH1D *ZPtBinsH = (TH1D *)file->Get("ZPtBinsH");
    if (ZPtBinsH == NULL) {
        Logger::get("RecoilCorrector")
            ->debug("File should contain histogram with the name ZPtBinsH");
        Logger::get("RecoilCorrector")
            ->debug("Check content of the file {}", fileName);
        exit(-1);
    }
    int nZPtBins = ZPtBinsH->GetNbinsX();
    float ZPtBins[10];
    TString ZPtStr[10];
    for (int i = 0; i <= nZPtBins; ++i) {
        ZPtBins[i] = ZPtBinsH->GetXaxis()->GetBinLowEdge(i + 1);
        if (i < nZPtBins)
            ZPtStr[i] = ZPtBinsH->GetXaxis()->GetBinLabel(i + 1);
    }

    TH1D *nJetBinsH = (TH1D *)file->Get("nJetBinsH");
    if (nJetBinsH == NULL) {
        Logger::get("RecoilCorrector")
            ->debug("File should contain histogram with the name nJetBinsH");
        Logger::get("RecoilCorrector")
            ->debug("Check content of the file {}", fileName);
        exit(-1);
    }
    int nJetsBins = nJetBinsH->GetNbinsX();
    TString nJetsStr[5];
    for (int i = 0; i < nJetsBins; ++i) {
        nJetsStr[i] = nJetBinsH->GetXaxis()->GetBinLabel(i + 1);
    }
    InitMEtWeights(file, perpZStr, paralZStr, nZPtBins, ZPtBins, ZPtStr,
                   nJetsBins, nJetsStr);
    _epsrel = 5e-4;
    _epsabs = 5e-4;
    _range = 0.95;
}

RecoilCorrector::~RecoilCorrector() {}

void RecoilCorrector::InitMEtWeights(TFile *_file, TString _perpZStr,
                                     TString _paralZStr, int nZPtBins,
                                     float *ZPtBins, TString *_ZPtStr,
                                     int nJetsBins, TString *_nJetsStr) {

    std::vector<float> newZPtBins;
    std::vector<std::string> newZPtStr;
    std::vector<std::string> newNJetsStr;

    std::string newPerpZStr = std::string(_perpZStr);
    std::string newParalZStr = std::string(_paralZStr);

    for (int idx = 0; idx < nZPtBins + 1; ++idx)
        newZPtBins.push_back(ZPtBins[idx]);
    for (int idx = 0; idx < nZPtBins; ++idx)
        newZPtStr.push_back(std::string(_ZPtStr[idx]));

    for (int idx = 0; idx < nJetsBins; ++idx)
        newNJetsStr.push_back(std::string(_nJetsStr[idx]));

    InitMEtWeights(_file, newZPtBins, newPerpZStr, newParalZStr, newZPtStr,
                   newNJetsStr);
}

void RecoilCorrector::InitMEtWeights(
    TFile *_fileMet, const std::vector<float> &ZPtBins,
    const std::string _perpZStr, const std::string _paralZStr,
    const std::vector<std::string> &_ZPtStr,
    const std::vector<std::string> &_nJetsStr) {

    // checking files
    if (_fileMet->IsZombie()) {
        Logger::get("RecoilCorrector")->debug("File {} is not found", fileName);
        Logger::get("RecoilCorrector")->debug("quitting program...");
        exit(-1);
    }

    _nZPtBins = ZPtBins.size() - 1; // the -1 is on purpose!
    _nJetsBins = _nJetsStr.size();
    _ZPtBins = ZPtBins;

    for (int ZPtBin = 0; ZPtBin < _nZPtBins; ++ZPtBin) {
        for (int jetBin = 0; jetBin < _nJetsBins; ++jetBin) {

            TString binStrPerpData =
                _perpZStr + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_data";
            TString binStrParalData = _paralZStr + "_" + _nJetsStr[jetBin] +
                                      _ZPtStr[ZPtBin] + "_data";
            TString binStrPerpMC =
                _perpZStr + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_mc";
            TString binStrParalMC =
                _paralZStr + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_mc";

            _metZParalData[ZPtBin][jetBin] =
                (TF1 *)_fileMet->Get(binStrParalData);
            _metZPerpData[ZPtBin][jetBin] =
                (TF1 *)_fileMet->Get(binStrPerpData);
            _metZParalMC[ZPtBin][jetBin] = (TF1 *)_fileMet->Get(binStrParalMC);
            _metZPerpMC[ZPtBin][jetBin] = (TF1 *)_fileMet->Get(binStrPerpMC);

            // checking functions
            if (_metZParalData[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Function with name {} is not found in file {} "
                            "quitting program...",
                            binStrParalData, fileName);
                exit(-1);
            }
            if (_metZPerpData[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Function with name {} is not found in file {} "
                            "quitting program...",
                            binStrPerpData, fileName);
                exit(-1);
            }

            if (_metZParalMC[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Function with name {} is not found in file {} "
                            "quitting program...",
                            binStrParalMC, fileName);

                exit(-1);
            }
            if (_metZPerpMC[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Function with name {} is not found in file {} "
                            "quitting program...",
                            binStrPerpMC, fileName);
                exit(-1);
            }

            TString binStrPerpDataHist = _perpZStr + "_" + _nJetsStr[jetBin] +
                                         _ZPtStr[ZPtBin] + "_hist_data";
            TString binStrParalDataHist = _paralZStr + "_" + _nJetsStr[jetBin] +
                                          _ZPtStr[ZPtBin] + "_hist_data";
            TString binStrPerpMCHist = _perpZStr + "_" + _nJetsStr[jetBin] +
                                       _ZPtStr[ZPtBin] + "_hist_mc";
            TString binStrParalMCHist = _paralZStr + "_" + _nJetsStr[jetBin] +
                                        _ZPtStr[ZPtBin] + "_hist_mc";

            _metZParalDataHist[ZPtBin][jetBin] =
                ((TH1D *)_fileMet->Get(binStrParalDataHist));
            _metZPerpDataHist[ZPtBin][jetBin] =
                ((TH1D *)_fileMet->Get(binStrPerpDataHist));
            _metZParalMCHist[ZPtBin][jetBin] =
                ((TH1D *)_fileMet->Get(binStrParalMCHist));
            _metZPerpMCHist[ZPtBin][jetBin] =
                ((TH1D *)_fileMet->Get(binStrPerpMCHist));

            // checking histograms
            if (_metZParalDataHist[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Histogram with name {} is not found in file {}... "
                            "quitting program...",
                            binStrParalDataHist, fileName);
                exit(-1);
            }
            if (_metZPerpDataHist[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Histogram with name {} is not found in file {}... "
                            "quitting program...",
                            binStrPerpDataHist, fileName);
                exit(-1);
            }

            if (_metZParalMCHist[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Histogram with name {} is not found in file {}... "
                            "quitting program...",
                            binStrParalMCHist, fileName);
                exit(-1);
            }
            if (_metZPerpMCHist[ZPtBin][jetBin] == NULL) {
                Logger::get("RecoilCorrector")
                    ->debug("Histogram with name {} is not found in file {}... "
                            "quitting program...",
                            binStrPerpMCHist, fileName);
                exit(-1);
            }

            Logger::get("RecoilCorrector")
                ->debug(" {} : {}", _ZPtStr[ZPtBin], _nJetsStr[jetBin]);

            double xminD, xmaxD;

            _metZParalData[ZPtBin][jetBin]->GetRange(xminD, xmaxD);
            _xminMetZParalData[ZPtBin][jetBin] = float(xminD);
            _xmaxMetZParalData[ZPtBin][jetBin] = float(xmaxD);

            _metZPerpData[ZPtBin][jetBin]->GetRange(xminD, xmaxD);
            _xminMetZPerpData[ZPtBin][jetBin] = float(xminD);
            _xmaxMetZPerpData[ZPtBin][jetBin] = float(xmaxD);

            _metZParalMC[ZPtBin][jetBin]->GetRange(xminD, xmaxD);
            _xminMetZParalMC[ZPtBin][jetBin] = float(xminD);
            _xmaxMetZParalMC[ZPtBin][jetBin] = float(xmaxD);

            _metZPerpMC[ZPtBin][jetBin]->GetRange(xminD, xmaxD);
            _xminMetZPerpMC[ZPtBin][jetBin] = float(xminD);
            _xmaxMetZPerpMC[ZPtBin][jetBin] = float(xmaxD);

            _xminMetZParal[ZPtBin][jetBin] =
                TMath::Max(_xminMetZParalData[ZPtBin][jetBin],
                           _xminMetZParalMC[ZPtBin][jetBin]);
            _xmaxMetZParal[ZPtBin][jetBin] =
                TMath::Min(_xmaxMetZParalData[ZPtBin][jetBin],
                           _xmaxMetZParalMC[ZPtBin][jetBin]);

            _xminMetZPerp[ZPtBin][jetBin] =
                TMath::Max(_xminMetZPerpData[ZPtBin][jetBin],
                           _xminMetZPerpMC[ZPtBin][jetBin]);
            _xmaxMetZPerp[ZPtBin][jetBin] =
                TMath::Min(_xmaxMetZPerpData[ZPtBin][jetBin],
                           _xmaxMetZPerpMC[ZPtBin][jetBin]);

            _meanMetZParalData[ZPtBin][jetBin] =
                _metZParalData[ZPtBin][jetBin]->Mean(
                    _xminMetZParalData[ZPtBin][jetBin],
                    _xmaxMetZParalData[ZPtBin][jetBin]);
            _rmsMetZParalData[ZPtBin][jetBin] =
                TMath::Sqrt(_metZParalData[ZPtBin][jetBin]->CentralMoment(
                    2, _xminMetZParalData[ZPtBin][jetBin],
                    _xmaxMetZParalData[ZPtBin][jetBin]));
            _meanMetZPerpData[ZPtBin][jetBin] = 0;
            _rmsMetZPerpData[ZPtBin][jetBin] =
                TMath::Sqrt(_metZPerpData[ZPtBin][jetBin]->CentralMoment(
                    2, _xminMetZPerpData[ZPtBin][jetBin],
                    _xmaxMetZPerpData[ZPtBin][jetBin]));

            _meanMetZParalMC[ZPtBin][jetBin] =
                _metZParalMC[ZPtBin][jetBin]->Mean(
                    _xminMetZParalMC[ZPtBin][jetBin],
                    _xmaxMetZParalMC[ZPtBin][jetBin]);
            _rmsMetZParalMC[ZPtBin][jetBin] =
                TMath::Sqrt(_metZParalMC[ZPtBin][jetBin]->CentralMoment(
                    2, _xminMetZParalMC[ZPtBin][jetBin],
                    _xmaxMetZParalMC[ZPtBin][jetBin]));
            _meanMetZPerpMC[ZPtBin][jetBin] = 0;
            _rmsMetZPerpMC[ZPtBin][jetBin] =
                TMath::Sqrt(_metZPerpMC[ZPtBin][jetBin]->CentralMoment(
                    2, _xminMetZPerpMC[ZPtBin][jetBin],
                    _xmaxMetZPerpMC[ZPtBin][jetBin]));
        }
    }
}

void RecoilCorrector::CorrectWithHist(float MetPx, float MetPy, float genVPx,
                                      float genVPy, float visVPx, float visVPy,
                                      int njets, float &MetCorrPx,
                                      float &MetCorrPy) {

    // input parameters
    // MetPx, MetPy - missing transverse momentum
    // genVPx, genVPy - generated transverse momentum of Z(W)
    // visVPx, visVPy - visible transverse momentum of Z(W)
    // njets - number of jets
    // MetCorrPx, MetCorrPy - corrected missing transverse momentum

    Double_t Zpt = TMath::Sqrt(genVPx * genVPx + genVPy * genVPy);

    Double_t U1 = 0.0;
    Double_t U2 = 0.0;
    Double_t metU1 = 0.0;
    Double_t metU2 = 0.0;

    CalculateU1U2FromMet(MetPx, MetPy, genVPx, genVPy, visVPx, visVPy, U1, U2,
                         metU1, metU2);
    if (Zpt > 1000.0)
        Zpt = 999.0;
    if (njets >= _nJetsBins)
        njets = _nJetsBins - 1;

    int ZptBin = binNumber(Zpt, _ZPtBins);

    TH1D *metZParalDataHist = ((TH1D *)_metZParalDataHist[ZptBin][njets]);
    TH1D *metZPerpDataHist = ((TH1D *)_metZPerpDataHist[ZptBin][njets]);

    TH1D *metZParalMCHist = ((TH1D *)_metZParalMCHist[ZptBin][njets]);
    TH1D *metZPerpMCHist = ((TH1D *)_metZPerpMCHist[ZptBin][njets]);

    if (U1 > _range * _xminMetZParal[ZptBin][njets] &&
        U1 < _range * _xmaxMetZParal[ZptBin][njets]) {

        int nSumProb = 1;
        double q[1];
        double sumProb[1];

        const int ibin = metZParalMCHist->FindBin(U1);
        const double integralToNextBinEdge =
            metZParalMCHist->GetBinContent(ibin) *
            (metZParalMCHist->GetBinLowEdge(ibin + 1) - U1) /
            metZParalMCHist->GetBinWidth(ibin);
        sumProb[0] =
            (metZParalMCHist->Integral(1, ibin) - integralToNextBinEdge) /
            metZParalMCHist->Integral();
        Logger::get("RecoilCorrector")
            ->debug("U1 value: {} bin in MC hist: {}Integral: {}", U1,
                    metZParalMCHist->FindBin(U1), sumProb[0]);

        if (sumProb[0] < 0) {
            Logger::get("RecoilCorrector")
                ->debug("Warning ! ProbSum[0] = {}", sumProb[0]);
            sumProb[0] = 1e-5;
        }
        if (sumProb[0] > 1) {
            Logger::get("RecoilCorrector")
                ->debug("Warning ! ProbSum[0] = {}", sumProb[0]);
            sumProb[0] = 1.0 - 1e-5;
        }
        metZParalDataHist->GetQuantiles(nSumProb, q, sumProb);
        Logger::get("RecoilCorrector")
            ->debug("Parallel component. Detemined probability: {} Projection "
                    "value. old = {}",
                    sumProb[0], U1);
        float U1reco = float(q[0]);
        U1 = U1reco;
        Logger::get("RecoilCorrector")->debug(" new = {}", U1);

    } else {
        Logger::get("RecoilCorrector")
            ->debug("Warning: parallel Met component out of histogram range: "
                    "{}. Correction won't be applied",
                    U1);
        //  float U1reco = rescale(U1,
        //      		   _meanMetZParalData[ZptBin][njets],
        //      		   _meanMetZParalMC[ZptBin][njets],
        //      		   _rmsMetZParalData[ZptBin][njets],
        //      		   _rmsMetZParalMC[ZptBin][njets]);
        //  U1 = U1reco;
    }

    if (std::abs(U2) < _range * _xmaxMetZPerp[ZptBin][njets]) {

        int nSumProb = 1;
        double q[1];
        double sumProb[1];

        Logger::get("RecoilCorrector")
            ->debug("U2 value: {} bin in MC hist: {}", U2,
                    metZParalMCHist->FindBin(U2));
        const double absU2 = std::abs(U2);
        const int signU2 = TMath::Sign(1.0, U2);
        const int ibin = metZPerpMCHist->FindBin(absU2);
        const double integralToNextBinEdge =
            metZPerpMCHist->GetBinContent(ibin) *
            (metZPerpMCHist->GetBinLowEdge(ibin + 1) - absU2) /
            metZPerpMCHist->GetBinWidth(ibin);
        sumProb[0] =
            ((metZPerpMCHist->Integral(1, ibin) - integralToNextBinEdge) /
             metZPerpMCHist->Integral());
        if (sumProb[0] < 0) {
            Logger::get("RecoilCorrector")
                ->debug("Warning ! ProbSum[0] = {}", sumProb[0]);
            sumProb[0] = 1e-5;
        }
        if (sumProb[0] > 1) {
            Logger::get("RecoilCorrector")
                ->debug("Warning ! ProbSum[0] = {}", sumProb[0]);
            sumProb[0] = 1.0 - 1e-5;
        }
        metZPerpDataHist->GetQuantiles(nSumProb, q, sumProb);
        Logger::get("RecoilCorrector")
            ->debug("Perpendicular component. Determined probability: {} "
                    "Projection value. old = {}",
                    sumProb[0], U2);
        float U2reco = float(q[0]) * signU2;
        U2 = U2reco;
        Logger::get("RecoilCorrector")->debug(" new = {}", U2);

    } else {
        Logger::get("RecoilCorrector")
            ->debug("Warning: perpendicular Met component out of histogram "
                    "range: {}. Correction won't be applied",
                    U2);
        //  float U2reco = rescale(U2,
        //      		   _meanMetZPerpData[ZptBin][njets],
        //      		   _meanMetZPerpMC[ZptBin][njets],
        //      		   _rmsMetZPerpData[ZptBin][njets],
        //      		   _rmsMetZPerpMC[ZptBin][njets]);
        //  U2 = U2reco;
    }

    CalculateMetFromU1U2(U1, U2, genVPx, genVPy, visVPx, visVPy, MetCorrPx,
                         MetCorrPy);
}

void RecoilCorrector::CalculateU1U2FromMet(float metPx, float metPy,
                                           float genZPx, float genZPy,
                                           float diLepPx, float diLepPy,
                                           Double_t &U1, Double_t &U2,
                                           Double_t &metU1, Double_t &metU2) {

    auto diLep = ROOT::Math::XYVector(diLepPx, diLepPy);
    auto genZ = ROOT::Math::XYVector(genZPx, genZPy);
    auto met = ROOT::Math::XYVector(metPx, metPy);
    auto hadRec = met + diLep - genZ; // actually, that's: - 1.0 * ( hadrons +
                                      // Z) = Met - neutrinos (if involved)

    auto deltaPhiZHadRec = ROOT::Math::VectorUtil::DeltaPhi(genZ, hadRec);
    auto deltaPhiDiLepMEt = ROOT::Math::VectorUtil::DeltaPhi(diLep, met);

    U1 = hadRec.R() * TMath::Cos(deltaPhiZHadRec);
    U2 = hadRec.R() * TMath::Sin(deltaPhiZHadRec);

    metU1 = met.R() * TMath::Cos(deltaPhiDiLepMEt);
    metU2 = met.R() * TMath::Sin(deltaPhiDiLepMEt);
}

void RecoilCorrector::CalculateMetFromU1U2(float U1, float U2, float genZPx,
                                           float genZPy, float diLepPx,
                                           float diLepPy, float &metPx,
                                           float &metPy) {

    float hadRecPt = TMath::Sqrt(U1 * U1 + U2 * U2);

    float deltaPhiZHadRec = TMath::ATan2(U2, U1);

    float phiZ = TMath::ATan2(genZPy, genZPx);

    float phiHadRec = phiZ + deltaPhiZHadRec;

    float hadRecX = hadRecPt * TMath::Cos(phiHadRec);
    float hadRecY = hadRecPt * TMath::Sin(phiHadRec);

    metPx = hadRecX + genZPx - diLepPx;
    metPy = hadRecY + genZPy - diLepPy;
}