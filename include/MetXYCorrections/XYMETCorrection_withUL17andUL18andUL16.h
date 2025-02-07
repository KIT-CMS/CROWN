#ifndef XYMETCorrection_H
#define XYMETCorrection_H

#include <utility>
#include <cmath>
#include "TString.h"
#include "TMath.h"

enum TheRunEra {
    y2016B,
    y2016C,
    y2016D,
    y2016E,
    y2016F,
    y2016G,
    y2016H,
    y2017B,
    y2017C,
    y2017D,
    y2017E,
    y2017F,
    y2018A,
    y2018B,
    y2018C,
    y2018D,
    y2016MC,
    y2017MC,
    y2018MC,
    yUL2016B,
    yUL2016C,
    yUL2016D,
    yUL2016E,
    yUL2016F,
    yUL2016Flate,
    yUL2016G,
    yUL2016H,
    yUL2017B,
    yUL2017C,
    yUL2017D,
    yUL2017E,
    yUL2017F,
    yUL2018A,
    yUL2018B,
    yUL2018C,
    yUL2018D,
    yUL2016MCAPV,
    yUL2016MCnonAPV,
    yUL2017MC,
    yUL2018MC
};

std::pair<double, double> METXYCorr_Met_MetPhi(double uncormet,
                                               double uncormet_phi, int runnb,
                                               TString year, bool isMC, int npv,
                                               bool isUL = false,
                                               bool ispuppi = false) {

    std::pair<double, double> TheXYCorr_Met_MetPhi(uncormet, uncormet_phi);

    if (npv > 100)
        npv = 100;
    int runera = -1;
    bool usemetv2 = false;
    if (isMC && year == "2016" && !isUL)
        runera = y2016MC;
    else if (isMC && year == "2017" && !isUL) {
        runera = y2017MC;
        usemetv2 = true;
    } else if (isMC && year == "2018" && !isUL)
        runera = y2018MC;
    else if (isMC && year == "2016APV" && isUL)
        runera = yUL2016MCAPV;
    else if (isMC && year == "2016nonAPV" && isUL)
        runera = yUL2016MCnonAPV;
    else if (isMC && year == "2017" && isUL)
        runera = yUL2017MC;
    else if (isMC && year == "2018" && isUL)
        runera = yUL2018MC;

    else if (!isMC && runnb >= 272007 && runnb <= 275376 && !isUL)
        runera = y2016B;
    else if (!isMC && runnb >= 275657 && runnb <= 276283 && !isUL)
        runera = y2016C;
    else if (!isMC && runnb >= 276315 && runnb <= 276811 && !isUL)
        runera = y2016D;
    else if (!isMC && runnb >= 276831 && runnb <= 277420 && !isUL)
        runera = y2016E;
    else if (!isMC && runnb >= 277772 && runnb <= 278808 && !isUL)
        runera = y2016F;
    else if (!isMC && runnb >= 278820 && runnb <= 280385 && !isUL)
        runera = y2016G;
    else if (!isMC && runnb >= 280919 && runnb <= 284044 && !isUL)
        runera = y2016H;

    else if (!isMC && runnb >= 297020 && runnb <= 299329 && !isUL) {
        runera = y2017B;
        usemetv2 = true;
    } else if (!isMC && runnb >= 299337 && runnb <= 302029 && !isUL) {
        runera = y2017C;
        usemetv2 = true;
    } else if (!isMC && runnb >= 302030 && runnb <= 303434 && !isUL) {
        runera = y2017D;
        usemetv2 = true;
    } else if (!isMC && runnb >= 303435 && runnb <= 304826 && !isUL) {
        runera = y2017E;
        usemetv2 = true;
    } else if (!isMC && runnb >= 304911 && runnb <= 306462 && !isUL) {
        runera = y2017F;
        usemetv2 = true;
    }

    else if (!isMC && runnb >= 315252 && runnb <= 316995 && !isUL)
        runera = y2018A;
    else if (!isMC && runnb >= 316998 && runnb <= 319312 && !isUL)
        runera = y2018B;
    else if (!isMC && runnb >= 319313 && runnb <= 320393 && !isUL)
        runera = y2018C;
    else if (!isMC && runnb >= 320394 && runnb <= 325273 && !isUL)
        runera = y2018D;

    else if (!isMC && runnb >= 315252 && runnb <= 316995 && isUL)
        runera = yUL2018A;
    else if (!isMC && runnb >= 316998 && runnb <= 319312 && isUL)
        runera = yUL2018B;
    else if (!isMC && runnb >= 319313 && runnb <= 320393 && isUL)
        runera = yUL2018C;
    else if (!isMC && runnb >= 320394 && runnb <= 325273 && isUL)
        runera = yUL2018D;

    else if (!isMC && runnb >= 297020 && runnb <= 299329 && isUL) {
        runera = yUL2017B;
        usemetv2 = false;
    } else if (!isMC && runnb >= 299337 && runnb <= 302029 && isUL) {
        runera = yUL2017C;
        usemetv2 = false;
    } else if (!isMC && runnb >= 302030 && runnb <= 303434 && isUL) {
        runera = yUL2017D;
        usemetv2 = false;
    } else if (!isMC && runnb >= 303435 && runnb <= 304826 && isUL) {
        runera = yUL2017E;
        usemetv2 = false;
    } else if (!isMC && runnb >= 304911 && runnb <= 306462 && isUL) {
        runera = yUL2017F;
        usemetv2 = false;
    }

    else if (!isMC && runnb >= 272007 && runnb <= 275376 && isUL)
        runera = yUL2016B;
    else if (!isMC && runnb >= 275657 && runnb <= 276283 && isUL)
        runera = yUL2016C;
    else if (!isMC && runnb >= 276315 && runnb <= 276811 && isUL)
        runera = yUL2016D;
    else if (!isMC && runnb >= 276831 && runnb <= 277420 && isUL)
        runera = yUL2016E;
    else if (!isMC &&
             ((runnb >= 277772 && runnb <= 278768) || runnb == 278770) && isUL)
        runera = yUL2016F;
    else if (!isMC &&
             ((runnb >= 278801 && runnb <= 278808) || runnb == 278769) && isUL)
        runera = yUL2016Flate;
    else if (!isMC && runnb >= 278820 && runnb <= 280385 && isUL)
        runera = yUL2016G;
    else if (!isMC && runnb >= 280919 && runnb <= 284044 && isUL)
        runera = yUL2016H;

    else {
        // Couldn't find data/MC era => no correction applied
        return TheXYCorr_Met_MetPhi;
    }

    double METxcorr(0.), METycorr(0.);

    if (!usemetv2) { // Current recommendation for 2016 and 2018
        if (!ispuppi) {
            if (runera == y2016B)
                METxcorr = -(-0.0478335 * npv - 0.108032);
            if (runera == y2016B)
                METycorr = -(0.125148 * npv + 0.355672);
            if (runera == y2016C)
                METxcorr = -(-0.0916985 * npv + 0.393247);
            if (runera == y2016C)
                METycorr = -(0.151445 * npv + 0.114491);
            if (runera == y2016D)
                METxcorr = -(-0.0581169 * npv + 0.567316);
            if (runera == y2016D)
                METycorr = -(0.147549 * npv + 0.403088);
            if (runera == y2016E)
                METxcorr = -(-0.065622 * npv + 0.536856);
            if (runera == y2016E)
                METycorr = -(0.188532 * npv + 0.495346);
            if (runera == y2016F)
                METxcorr = -(-0.0313322 * npv + 0.39866);
            if (runera == y2016F)
                METycorr = -(0.16081 * npv + 0.960177);
            if (runera == y2016G)
                METxcorr = -(0.040803 * npv - 0.290384);
            if (runera == y2016G)
                METycorr = -(0.0961935 * npv + 0.666096);
            if (runera == y2016H)
                METxcorr = -(0.0330868 * npv - 0.209534);
            if (runera == y2016H)
                METycorr = -(0.141513 * npv + 0.816732);
            if (runera == y2017B)
                METxcorr = -(-0.259456 * npv + 1.95372);
            if (runera == y2017B)
                METycorr = -(0.353928 * npv - 2.46685);
            if (runera == y2017C)
                METxcorr = -(-0.232763 * npv + 1.08318);
            if (runera == y2017C)
                METycorr = -(0.257719 * npv - 1.1745);
            if (runera == y2017D)
                METxcorr = -(-0.238067 * npv + 1.80541);
            if (runera == y2017D)
                METycorr = -(0.235989 * npv - 1.44354);
            if (runera == y2017E)
                METxcorr = -(-0.212352 * npv + 1.851);
            if (runera == y2017E)
                METycorr = -(0.157759 * npv - 0.478139);
            if (runera == y2017F)
                METxcorr = -(-0.232733 * npv + 2.24134);
            if (runera == y2017F)
                METycorr = -(0.213341 * npv + 0.684588);
            if (runera == y2018A)
                METxcorr = -(0.362865 * npv - 1.94505);
            if (runera == y2018A)
                METycorr = -(0.0709085 * npv - 0.307365);
            if (runera == y2018B)
                METxcorr = -(0.492083 * npv - 2.93552);
            if (runera == y2018B)
                METycorr = -(0.17874 * npv - 0.786844);
            if (runera == y2018C)
                METxcorr = -(0.521349 * npv - 1.44544);
            if (runera == y2018C)
                METycorr = -(0.118956 * npv - 1.96434);
            if (runera == y2018D)
                METxcorr = -(0.531151 * npv - 1.37568);
            if (runera == y2018D)
                METycorr = -(0.0884639 * npv - 1.57089);
            if (runera == y2016MC)
                METxcorr = -(-0.195191 * npv - 0.170948);
            if (runera == y2016MC)
                METycorr = -(-0.0311891 * npv + 0.787627);
            if (runera == y2017MC)
                METxcorr = -(-0.217714 * npv + 0.493361);
            if (runera == y2017MC)
                METycorr = -(0.177058 * npv - 0.336648);
            if (runera == y2018MC)
                METxcorr = -(0.296713 * npv - 0.141506);
            if (runera == y2018MC)
                METycorr = -(0.115685 * npv + 0.0128193);

            // UL2017
            if (runera == yUL2017B)
                METxcorr = -(-0.211161 * npv + 0.419333);
            if (runera == yUL2017B)
                METycorr = -(0.251789 * npv + -1.28089);
            if (runera == yUL2017C)
                METxcorr = -(-0.185184 * npv + -0.164009);
            if (runera == yUL2017C)
                METycorr = -(0.200941 * npv + -0.56853);
            if (runera == yUL2017D)
                METxcorr = -(-0.201606 * npv + 0.426502);
            if (runera == yUL2017D)
                METycorr = -(0.188208 * npv + -0.58313);
            if (runera == yUL2017E)
                METxcorr = -(-0.162472 * npv + 0.176329);
            if (runera == yUL2017E)
                METycorr = -(0.138076 * npv + -0.250239);
            if (runera == yUL2017F)
                METxcorr = -(-0.210639 * npv + 0.72934);
            if (runera == yUL2017F)
                METycorr = -(0.198626 * npv + 1.028);
            if (runera == yUL2017MC)
                METxcorr = -(-0.300155 * npv + 1.90608);
            if (runera == yUL2017MC)
                METycorr = -(0.300213 * npv + -2.02232);

            // UL2018
            if (runera == yUL2018A)
                METxcorr = -(0.263733 * npv + -1.91115);
            if (runera == yUL2018A)
                METycorr = -(0.0431304 * npv + -0.112043);
            if (runera == yUL2018B)
                METxcorr = -(0.400466 * npv + -3.05914);
            if (runera == yUL2018B)
                METycorr = -(0.146125 * npv + -0.533233);
            if (runera == yUL2018C)
                METxcorr = -(0.430911 * npv + -1.42865);
            if (runera == yUL2018C)
                METycorr = -(0.0620083 * npv + -1.46021);
            if (runera == yUL2018D)
                METxcorr = -(0.457327 * npv + -1.56856);
            if (runera == yUL2018D)
                METycorr = -(0.0684071 * npv + -0.928372);
            if (runera == yUL2018MC)
                METxcorr = -(0.183518 * npv + 0.546754);
            if (runera == yUL2018MC)
                METycorr = -(0.192263 * npv + -0.42121);

            // UL2016
            if (runera == yUL2016B)
                METxcorr = -(-0.0214894 * npv + -0.188255);
            if (runera == yUL2016B)
                METycorr = -(0.0876624 * npv + 0.812885);
            if (runera == yUL2016C)
                METxcorr = -(-0.032209 * npv + 0.067288);
            if (runera == yUL2016C)
                METycorr = -(0.113917 * npv + 0.743906);
            if (runera == yUL2016D)
                METxcorr = -(-0.0293663 * npv + 0.21106);
            if (runera == yUL2016D)
                METycorr = -(0.11331 * npv + 0.815787);
            if (runera == yUL2016E)
                METxcorr = -(-0.0132046 * npv + 0.20073);
            if (runera == yUL2016E)
                METycorr = -(0.134809 * npv + 0.679068);
            if (runera == yUL2016F)
                METxcorr = -(-0.0543566 * npv + 0.816597);
            if (runera == yUL2016F)
                METycorr = -(0.114225 * npv + 1.17266);
            if (runera == yUL2016Flate)
                METxcorr = -(0.134616 * npv + -0.89965);
            if (runera == yUL2016Flate)
                METycorr = -(0.0397736 * npv + 1.0385);
            if (runera == yUL2016G)
                METxcorr = -(0.121809 * npv + -0.584893);
            if (runera == yUL2016G)
                METycorr = -(0.0558974 * npv + 0.891234);
            if (runera == yUL2016H)
                METxcorr = -(0.0868828 * npv + -0.703489);
            if (runera == yUL2016H)
                METycorr = -(0.0888774 * npv + 0.902632);
            if (runera == yUL2016MCnonAPV)
                METxcorr = -(-0.153497 * npv + -0.231751);
            if (runera == yUL2016MCnonAPV)
                METycorr = -(0.00731978 * npv + 0.243323);
            if (runera == yUL2016MCAPV)
                METxcorr = -(-0.188743 * npv + 0.136539);
            if (runera == yUL2016MCAPV)
                METycorr = -(0.0127927 * npv + 0.117747);
        }

        // UL2017Puppi
        if (ispuppi) {
            if (runera == yUL2017B)
                METxcorr = -(-0.00382117 * npv + -0.666228);
            if (runera == yUL2017B)
                METycorr = -(0.0109034 * npv + 0.172188);
            if (runera == yUL2017C)
                METxcorr = -(-0.00110699 * npv + -0.747643);
            if (runera == yUL2017C)
                METycorr = -(-0.0012184 * npv + 0.303817);
            if (runera == yUL2017D)
                METxcorr = -(-0.00141442 * npv + -0.721382);
            if (runera == yUL2017D)
                METycorr = -(-0.0011873 * npv + 0.21646);
            if (runera == yUL2017E)
                METxcorr = -(0.00593859 * npv + -0.851999);
            if (runera == yUL2017E)
                METycorr = -(-0.00754254 * npv + 0.245956);
            if (runera == yUL2017F)
                METxcorr = -(0.00765682 * npv + -0.945001);
            if (runera == yUL2017F)
                METycorr = -(-0.0154974 * npv + 0.804176);
            if (runera == yUL2017MC)
                METxcorr = -(-0.0102265 * npv + -0.446416);
            if (runera == yUL2017MC)
                METycorr = -(0.0198663 * npv + 0.243182);

            // UL2018Puppi
            if (runera == yUL2018A)
                METxcorr = -(-0.0073377 * npv + 0.0250294);
            if (runera == yUL2018A)
                METycorr = -(-0.000406059 * npv + 0.0417346);
            if (runera == yUL2018B)
                METxcorr = -(0.00434261 * npv + 0.00892927);
            if (runera == yUL2018B)
                METycorr = -(0.00234695 * npv + 0.20381);
            if (runera == yUL2018C)
                METxcorr = -(0.00198311 * npv + 0.37026);
            if (runera == yUL2018C)
                METycorr = -(-0.016127 * npv + 0.402029);
            if (runera == yUL2018D)
                METxcorr = -(0.00220647 * npv + 0.378141);
            if (runera == yUL2018D)
                METycorr = -(-0.0160244 * npv + 0.471053);
            if (runera == yUL2018MC)
                METxcorr = -(-0.0214557 * npv + 0.969428);
            if (runera == yUL2018MC)
                METycorr = -(0.0167134 * npv + 0.199296);

            // UL2016Puppi
            if (runera == yUL2016B)
                METxcorr = -(-0.00109025 * npv + -0.338093);
            if (runera == yUL2016B)
                METycorr = -(-0.00356058 * npv + 0.128407);
            if (runera == yUL2016C)
                METxcorr = -(-0.00271913 * npv + -0.342268);
            if (runera == yUL2016C)
                METycorr = -(0.00187386 * npv + 0.104);
            if (runera == yUL2016D)
                METxcorr = -(-0.00254194 * npv + -0.305264);
            if (runera == yUL2016D)
                METycorr = -(-0.00177408 * npv + 0.164639);
            if (runera == yUL2016E)
                METxcorr = -(-0.00358835 * npv + -0.225435);
            if (runera == yUL2016E)
                METycorr = -(-0.000444268 * npv + 0.180479);
            if (runera == yUL2016F)
                METxcorr = -(0.0056759 * npv + -0.454101);
            if (runera == yUL2016F)
                METycorr = -(-0.00962707 * npv + 0.35731);
            if (runera == yUL2016Flate)
                METxcorr = -(0.0234421 * npv + -0.371298);
            if (runera == yUL2016Flate)
                METycorr = -(-0.00997438 * npv + 0.0809178);
            if (runera == yUL2016G)
                METxcorr = -(0.0182134 * npv + -0.335786);
            if (runera == yUL2016G)
                METycorr = -(-0.0063338 * npv + 0.093349);
            if (runera == yUL2016H)
                METxcorr = -(0.015702 * npv + -0.340832);
            if (runera == yUL2016H)
                METycorr = -(-0.00544957 * npv + 0.199093);
            if (runera == yUL2016MCnonAPV)
                METxcorr = -(-0.0058341 * npv + -0.395049);
            if (runera == yUL2016MCnonAPV)
                METycorr = -(0.00971595 * npv + -0.101288);
            if (runera == yUL2016MCAPV)
                METxcorr = -(-0.0060447 * npv + -0.4183);
            if (runera == yUL2016MCAPV)
                METycorr = -(0.008331 * npv + -0.0990046);
        }

    } else { // these are the corrections for v2 MET recipe (currently
             // recommended for 2017)
        if (runera == y2016B)
            METxcorr = -(-0.0374977 * npv + 0.00488262);
        if (runera == y2016B)
            METycorr = -(0.107373 * npv + -0.00732239);
        if (runera == y2016C)
            METxcorr = -(-0.0832562 * npv + 0.550742);
        if (runera == y2016C)
            METycorr = -(0.142469 * npv + -0.153718);
        if (runera == y2016D)
            METxcorr = -(-0.0400931 * npv + 0.753734);
        if (runera == y2016D)
            METycorr = -(0.127154 * npv + 0.0175228);
        if (runera == y2016E)
            METxcorr = -(-0.0409231 * npv + 0.755128);
        if (runera == y2016E)
            METycorr = -(0.168407 * npv + 0.126755);
        if (runera == y2016F)
            METxcorr = -(-0.0161259 * npv + 0.516919);
        if (runera == y2016F)
            METycorr = -(0.141176 * npv + 0.544062);
        if (runera == y2016G)
            METxcorr = -(0.0583851 * npv + -0.0987447);
        if (runera == y2016G)
            METycorr = -(0.0641427 * npv + 0.319112);
        if (runera == y2016H)
            METxcorr = -(0.0706267 * npv + -0.13118);
        if (runera == y2016H)
            METycorr = -(0.127481 * npv + 0.370786);
        if (runera == y2017B)
            METxcorr = -(-0.19563 * npv + 1.51859);
        if (runera == y2017B)
            METycorr = -(0.306987 * npv + -1.84713);
        if (runera == y2017C)
            METxcorr = -(-0.161661 * npv + 0.589933);
        if (runera == y2017C)
            METycorr = -(0.233569 * npv + -0.995546);
        if (runera == y2017D)
            METxcorr = -(-0.180911 * npv + 1.23553);
        if (runera == y2017D)
            METycorr = -(0.240155 * npv + -1.27449);
        if (runera == y2017E)
            METxcorr = -(-0.149494 * npv + 0.901305);
        if (runera == y2017E)
            METycorr = -(0.178212 * npv + -0.535537);
        if (runera == y2017F)
            METxcorr = -(-0.165154 * npv + 1.02018);
        if (runera == y2017F)
            METycorr = -(0.253794 * npv + 0.75776);
        if (runera == y2018A)
            METxcorr = -(0.362642 * npv + -1.55094);
        if (runera == y2018A)
            METycorr = -(0.0737842 * npv + -0.677209);
        if (runera == y2018B)
            METxcorr = -(0.485614 * npv + -2.45706);
        if (runera == y2018B)
            METycorr = -(0.181619 * npv + -1.00636);
        if (runera == y2018C)
            METxcorr = -(0.503638 * npv + -1.01281);
        if (runera == y2018C)
            METycorr = -(0.147811 * npv + -1.48941);
        if (runera == y2018D)
            METxcorr = -(0.520265 * npv + -1.20322);
        if (runera == y2018D)
            METycorr = -(0.143919 * npv + -0.979328);
        if (runera == y2016MC)
            METxcorr = -(-0.159469 * npv + -0.407022);
        if (runera == y2016MC)
            METycorr = -(-0.0405812 * npv + 0.570415);
        if (runera == y2017MC)
            METxcorr = -(-0.182569 * npv + 0.276542);
        if (runera == y2017MC)
            METycorr = -(0.155652 * npv + -0.417633);
        if (runera == y2018MC)
            METxcorr = -(0.299448 * npv + -0.13866);
        if (runera == y2018MC)
            METycorr = -(0.118785 * npv + 0.0889588);
    }

    double CorrectedMET_x = uncormet * cos(uncormet_phi) + METxcorr;
    double CorrectedMET_y = uncormet * sin(uncormet_phi) + METycorr;

    double CorrectedMET =
        sqrt(CorrectedMET_x * CorrectedMET_x + CorrectedMET_y * CorrectedMET_y);
    double CorrectedMETPhi;
    if (CorrectedMET_x == 0 && CorrectedMET_y > 0)
        CorrectedMETPhi = TMath::Pi();
    else if (CorrectedMET_x == 0 && CorrectedMET_y < 0)
        CorrectedMETPhi = -TMath::Pi();
    else if (CorrectedMET_x > 0)
        CorrectedMETPhi = TMath::ATan(CorrectedMET_y / CorrectedMET_x);
    else if (CorrectedMET_x < 0 && CorrectedMET_y > 0)
        CorrectedMETPhi =
            TMath::ATan(CorrectedMET_y / CorrectedMET_x) + TMath::Pi();
    else if (CorrectedMET_x < 0 && CorrectedMET_y < 0)
        CorrectedMETPhi =
            TMath::ATan(CorrectedMET_y / CorrectedMET_x) - TMath::Pi();
    else
        CorrectedMETPhi = 0;

    TheXYCorr_Met_MetPhi.first = CorrectedMET;
    TheXYCorr_Met_MetPhi.second = CorrectedMETPhi;
    return TheXYCorr_Met_MetPhi;
}

#endif
