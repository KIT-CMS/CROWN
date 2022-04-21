#include "../../include/RecoilCorrections/MetSystematics.hxx"
#include "../../include/utility/Logger.hxx"

MetSystematic::MetSystematic(std::string filepath) {

    fileName = filepath;
    TFile *file = new TFile(fileName, "READ");
    if (file->IsZombie()) {
        Logger::get("MetSystematics")
            ->debug("file {} is not found...   quitting ", fileName);
        exit(-1);
    }
    TH1D *jetBinsH = (TH1D *)file->Get("nJetBinsH");
    if (jetBinsH == NULL) {
        Logger::get("MetSystematics")
            ->debug("Histogram nJetBinsH should be contained in file {}",
                    fileName);
        Logger::get("MetSystematics")
            ->debug("Check content of the file {}", fileName);
        exit(-1);
    }

    nJetBins = jetBinsH->GetNbinsX();
    std::vector<TString> JetBins;
    for (int i = 0; i < nJetBins; ++i) {
        JetBins.push_back(jetBinsH->GetXaxis()->GetBinLabel(i + 1));
    }

    TString uncType[2] = {"Response", "Resolution"};

    TString histName = "syst";
    TH2D *hist = (TH2D *)file->Get(histName);
    if (hist == NULL) {
        Logger::get("MetSystematics")
            ->debug("Histogram {} should be contained in file {}", histName,
                    fileName);
        Logger::get("MetSystematics")
            ->debug("Check content of the file {}", fileName);
        exit(-1);
    }
    for (int xBin = 0; xBin < 2; ++xBin) {
        for (int yBin = 0; yBin < 3; ++yBin) {
            sysUnc[xBin][yBin] = hist->GetBinContent(xBin + 1, yBin + 1);
            Logger::get("MetSystematics")
                ->debug("Systematics : {} {} = {}", uncType[xBin],
                        JetBins[yBin], sysUnc[xBin][yBin]);
        }
    }

    for (int j = 0; j < nJetBins; ++j) {
        TString histName = JetBins[j];
        responseHist[j] = (TH1D *)file->Get(histName);
        if (responseHist[j] == NULL) {
            Logger::get("MetSystematics")
                ->debug("Histogram {} should be contained in file {}", histName,
                        fileName);
            Logger::get("MetSystematics")
                ->debug("Check content of the file {}", fileName);
            exit(-1);
        }
    }
}

void MetSystematic::ComputeHadRecoilFromMet(float metX, float metY,
                                            float genVPx, float genVPy,
                                            float visVPx, float visVPy,
                                            float &Hparal, float &Hperp) {

    float genVPt = TMath::Sqrt(genVPx * genVPx + genVPy * genVPy);
    float unitX = genVPx / genVPt;
    float unitY = genVPy / genVPt;

    float unitPhi = TMath::ATan2(unitY, unitX);
    float unitPerpX = TMath::Cos(unitPhi + 0.5 * TMath::Pi());
    float unitPerpY = TMath::Sin(unitPhi + 0.5 * TMath::Pi());

    float Hx = -metX - visVPx;
    float Hy = -metY - visVPy;

    Hparal = Hx * unitX + Hy * unitY;
    Hperp = Hx * unitPerpX + Hy * unitPerpY;
}

void MetSystematic::ComputeMetFromHadRecoil(float Hparal, float Hperp,
                                            float genVPx, float genVPy,
                                            float visVPx, float visVPy,
                                            float &metX, float &metY) {

    float genVPt = TMath::Sqrt(genVPx * genVPx + genVPy * genVPy);
    float unitX = genVPx / genVPt;
    float unitY = genVPy / genVPt;

    float unitPhi = TMath::ATan2(unitY, unitX);
    float unitPerpX = TMath::Cos(unitPhi + 0.5 * TMath::Pi());
    float unitPerpY = TMath::Sin(unitPhi + 0.5 * TMath::Pi());

    float det = unitX * unitPerpY - unitY * unitPerpX;
    float Hx = (Hparal * unitPerpY - Hperp * unitY) / det;
    float Hy = (Hperp * unitX - Hparal * unitPerpX) / det;

    metX = -Hx - visVPx;
    metY = -Hy - visVPy;
}

void MetSystematic::ShiftResponseMet(float metPx, float metPy, float genVPx,
                                     float genVPy, float visVPx, float visVPy,
                                     int njets, float sysShift,
                                     float &metShiftPx, float &metShiftPy) {

    float Hparal = 0;
    float Hperp = 0;
    float genVPt = TMath::Sqrt(genVPx * genVPx + genVPy * genVPy);

    // protection against null
    if (genVPt < 1.0) {
        metShiftPx = metPx;
        metShiftPy = metPy;
        return;
    }

    ComputeHadRecoilFromMet(metPx, metPy, genVPx, genVPy, visVPx, visVPy,
                            Hparal, Hperp);

    int jets = njets;
    if (jets > 2)
        jets = 2;
    if (jets < 0) {
        Logger::get("MetSystematics")->debug("Number of jets is negative !");
        exit(-1);
    }

    float mean = -responseHist[jets]->Interpolate(genVPt) * genVPt;
    float shift = sysShift * mean;
    Hparal = Hparal + (shift - mean);

    ComputeMetFromHadRecoil(Hparal, Hperp, genVPx, genVPy, visVPx, visVPy,
                            metShiftPx, metShiftPy);
}

void MetSystematic::ShiftResolutionMet(float metPx, float metPy, float genVPx,
                                       float genVPy, float visVPx, float visVPy,
                                       int njets, float sysShift,
                                       float &metShiftPx, float &metShiftPy) {

    float Hparal = 0;
    float Hperp = 0;
    float genVPt = TMath::Sqrt(genVPx * genVPx + genVPy * genVPy);

    // protection against null
    if (genVPt < 1.0) {
        metShiftPx = metPx;
        metShiftPy = metPy;
        return;
    }

    ComputeHadRecoilFromMet(metPx, metPy, genVPx, genVPy, visVPx, visVPy,
                            Hparal, Hperp);

    int jets = njets;
    if (jets > 2)
        jets = 2;
    if (jets < 0) {
        Logger::get("MetSystematics")->debug("Number of jets is negative !");
        exit(-1);
    }

    float mean = -responseHist[jets]->Interpolate(genVPt) * genVPt;
    Hperp = sysShift * Hperp;
    Hparal = mean + (Hparal - mean) * sysShift;

    ComputeMetFromHadRecoil(Hparal, Hperp, genVPx, genVPy, visVPx, visVPy,
                            metShiftPx, metShiftPy);
}

void MetSystematic::ShiftMet(float metPx, float metPy, float genVPx,
                             float genVPy, float visVPx, float visVPy,
                             int njets, int sysType, float sysShift,
                             float &metShiftPx, float &metShiftPy) {

    metShiftPx = metPx;
    metShiftPy = metPy;

    if (sysType == 0)
        // Responce Shift
        ShiftResponseMet(metPx, metPy, genVPx, genVPy, visVPx, visVPy, njets,
                         sysShift, metShiftPx, metShiftPy);
    else if (sysType == 1)
        // Resolution Shift
        ShiftResolutionMet(metPx, metPy, genVPx, genVPy, visVPx, visVPy, njets,
                           sysShift, metShiftPx, metShiftPy);
    else if (sysType == -1)
        Logger::get("MetSystematics")->debug("No type --> doing nothing");
    else {
        Logger::get("MetSystematics")
            ->debug("Unknown systematic type --> exiting");
        exit(-1);
    }
}

void MetSystematic::ApplyMetSystematic(float metPx, float metPy, float genVPx,
                                       float genVPy, float visVPx, float visVPy,
                                       int njets, int sysType, int sysShift,
                                       float &metShiftPx, float &metShiftPy) {

    int jets = njets;
    if (jets > 2)
        jets = 2;
    if (jets < 0) {
        Logger::get("MetSystematics")->debug("Number of jets is negative !");
        exit(-1);
    }
    if (sysShift == -1) {
        Logger::get("MetSystematics")
            ->debug("sysShift is Nominal, doing nothing !");
        exit(-1);
    }

    int type = 0;
    if (sysType != 0)
        type = 1;

    float scale = 1 + sysUnc[type][jets];
    if (sysShift != 0)
        scale = 1 - sysUnc[type][jets];

    ShiftMet(metPx, metPy, genVPx, genVPy, visVPx, visVPy, njets, type, scale,
             metShiftPx, metShiftPy);
}