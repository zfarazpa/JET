//
//  @file    chi2Min.C
//  @brief   Chi-squared minimisation for the ATLAS R32 alpha_S extraction
//           using MSHT20 PDFs and TMinuit.
//
//  Provides helper functions to fill systematic uncertainty vectors,
//  compute PDF/scale uncertainty bands, and perform the profile
//  likelihood fit via MIGRAD.
//

#include "common_flag.C"
#include "fcn.C"
#include "getPars.h"
#include "Style.h"

#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"

// ============================================================
//  Internal constants
// ============================================================

/// Number of bins in the UP/DOWN systematic arrays.
static constexpr unsigned int kNSys     = 56;
/// Number of MSHT20 PDF eigenvectors (pairs → 64 histograms total).
static constexpr unsigned int kNPDFMSHT = 64;
/// Number of CT18 PDF histograms (29 pairs).
static constexpr unsigned int kNPDFCT18 = 58;
/// Number of MC replica PDFs used in the replica method.
static constexpr unsigned int kNPDFRep  = 50;
/// Total number of fit parameters (1 alpha_S + 87 nuisance parameters).
static constexpr int           kNPar     = 88;
/// CT18 90% CL → 68% CL conversion factor.
static constexpr double        kCT18Norm = 1.645;

// ============================================================
//  Helper: select histogram and nuisance-parameter arrays
//          by pT region.
// ============================================================

namespace {

/// Return the data/theory/NP histogram corresponding to @p region (0..4).
/// Index: 0 = "incl", 1 = "pt1", ..., 4 = "pt4".
int regionIndex(const std::string& region)
{
    if (region == "incl") return 0;
    if (region == "pt1")  return 1;
    if (region == "pt2")  return 2;
    if (region == "pt3")  return 3;
    if (region == "pt4")  return 4;
    return -1;
}

/// Convenience: return the pT-region label for on-plot text.
const char* regionLabel(const std::string& region)
{
    if (region == "incl") return "p_{T3} > 60 GeV";
    if (region == "pt1")  return "p_{T3} > 0.05 x H_{T2}";
    if (region == "pt2")  return "p_{T3} > 0.1 x H_{T2}";
    if (region == "pt3")  return "p_{T3} > 0.2 x H_{T2}";
    if (region == "pt4")  return "p_{T3} > 0.3 x H_{T2}";
    return "";
}

/// y-axis range for the pre-fit main panel, keyed by region.
std::pair<double,double> prefitYRange(const std::string& region)
{
    if (region == "incl") return {0.10, 1.10};
    if (region == "pt1")  return {0.10, 0.90};
    if (region == "pt2")  return {0.03, 0.90};
    if (region == "pt3")  return {0.01, 0.50};
    if (region == "pt4")  return {0.001,0.12};
    return {0.0, 1.0};
}

} // anonymous namespace

// ============================================================
//  Uncertainty filling helpers
// ============================================================

/// @brief Compute symmetrised bin error from UP/DOWN variations and
///        push it onto the corresponding auxiliary vector.
///
/// For each systematic source s, the symmetrised value is
/// 0.5 * (|UP| + |DOWN|).
inline double symm(double up, double dn)
{
    return 0.5 * (std::fabs(up) + std::fabs(dn));
}

// ============================================================

void FillPoints(TH1D* histoData)
{
    std::cout << "in FillPoints function" << std::endl;
    for (int i = 0; i < histoData->GetNbinsX(); ++i)
        PointsAux->push_back(histoData->GetXaxis()->GetBinCenter(i + 1));
}

// ============================================================

void FillStatTheo(TH1D* histo, unsigned int nmin, unsigned int nmax)
{
    std::cout << "In FillStatTheo" << std::endl;
    for (unsigned int i = nmin; i < nmax; ++i)
        statTheoAux->push_back(histo->GetBinError(i + 1));
}

// ============================================================

void FillData(TH1D* histo, unsigned int nmin, unsigned int nmax)
{
    std::cout << "in FillData function" << std::endl;
    for (unsigned int i = nmin; i < nmax; ++i)
        dataPointsAux->push_back(histo->GetBinContent(i + 1));
}

// ============================================================

void FillDataToys(TH1D* histo, TH1D* histoE,
                  unsigned int nmin, unsigned int nmax, double iseed)
{
    TRandom3 gen(static_cast<UInt_t>(12345. * iseed));
    for (unsigned int i = nmin; i < nmax; ++i) {
        double bin  = histo ->GetBinContent(i + 1);
        double binE = histoE->GetBinContent(i + 1);
        dataPointsAux->push_back(gen.Gaus(bin, binE));
    }
}

// ============================================================

void FillNP(TH1D* npTheo)
{
    for (int i = 0; i < npTheo->GetNbinsX(); ++i)
        npAux->push_back(npTheo->GetBinContent(i + 1));
}

// ============================================================
//  FillVectSyst
//  Reads the 56-entry UP/DOWN systematic arrays for bins [nmin, nmax)
//  and pushes symmetrised values onto the per-source auxiliary vectors.
// ============================================================

void FillVectSyst(TH1D* hup[kNSys], TH1D* hdown[kNSys],
                  unsigned int nmin, unsigned int nmax)
{
    std::cout << "in FillVectSyst" << std::endl;

    for (unsigned int k = nmin; k < nmax; ++k) {
        const int b = k + 1;  // ROOT bin index

        // Stat and model uncertainties
        const double stat    = hup[0]->GetBinContent(b);
        const double ddncUP  = hup[1]->GetBinContent(b);
        const double ddncDN  = hdown[1]->GetBinContent(b);
        const double mod1UP  = hup[2]->GetBinContent(b),  mod1DN  = hdown[2]->GetBinContent(b);
        const double mod2UP  = hup[3]->GetBinContent(b),  mod2DN  = hdown[3]->GetBinContent(b);
        const double mod3UP  = hup[4]->GetBinContent(b),  mod3DN  = hdown[4]->GetBinContent(b);

        // JES 1-35
        double jesUP[35], jesDN[35];
        for (int j = 0; j < 35; ++j) {
            jesUP[j] = hup  [5 + j]->GetBinContent(b);
            jesDN[j] = hdown[5 + j]->GetBinContent(b);
        }

        // JER 1-13
        double jerUP[13], jerDN[13];
        for (int j = 0; j < 13; ++j) {
            jerUP[j] = hup  [40 + j]->GetBinContent(b);
            jerDN[j] = hdown[40 + j]->GetBinContent(b);
        }

        const double prwUP  = hup  [53]->GetBinContent(b);
        const double prwDN  = hdown[53]->GetBinContent(b);
        const double tileUP = hup  [55]->GetBinContent(b);
        const double tileDN = hdown[55]->GetBinContent(b);

        // Push onto auxiliary vectors
        statDataAux ->push_back(std::fabs(stat));
        ddncDataAux ->push_back(symm(ddncUP,  ddncDN));
        model1DataAux->push_back(symm(mod1UP,  mod1DN));
        model2DataAux->push_back(symm(mod2UP,  mod2DN));
        model3DataAux->push_back(symm(mod3UP,  mod3DN));

        std::vector<std::vector<double>*> jesAux = {
            jes1DataAux,  jes2DataAux,  jes3DataAux,  jes4DataAux,  jes5DataAux,
            jes6DataAux,  jes7DataAux,  jes8DataAux,  jes9DataAux,  jes10DataAux,
            jes11DataAux, jes12DataAux, jes13DataAux, jes14DataAux, jes15DataAux,
            jes16DataAux, jes17DataAux, jes18DataAux, jes19DataAux, jes20DataAux,
            jes21DataAux, jes22DataAux, jes23DataAux, jes24DataAux, jes25DataAux,
            jes26DataAux, jes27DataAux, jes28DataAux, jes29DataAux, jes30DataAux,
            jes31DataAux, jes32DataAux, jes33DataAux, jes34DataAux, jes35DataAux
        };
        for (int j = 0; j < 35; ++j)
            jesAux[j]->push_back(symm(jesUP[j], jesDN[j]));

        std::vector<std::vector<double>*> jerAux = {
            jer1DataAux,  jer2DataAux,  jer3DataAux,  jer4DataAux,  jer5DataAux,
            jer6DataAux,  jer7DataAux,  jer8DataAux,  jer9DataAux,  jer10DataAux,
            jer11DataAux, jer12DataAux, jer13DataAux
        };
        for (int j = 0; j < 13; ++j)
            jerAux[j]->push_back(symm(jerUP[j], jerDN[j]));

        prwDataAux ->push_back(symm(prwUP,  prwDN));
        tileDataAux->push_back(symm(tileUP, tileDN));
    }
}

// ============================================================
//  SetStatSysData
//  Compute the total uncertainty per bin (stat ⊕ all syst)
//  and write it as the bin error of hdata.
// ============================================================

void SetStatSysData(TH1D* hdata, TH1D* hup[51], TH1D* hdown[51])
{
    for (int k = 0; k < hdata->GetNbinsX(); ++k) {
        const int b = k + 1;

        const double stat   = std::sqrt(std::pow(hup[0]->GetBinContent(b), 2) +
                                        std::pow(hup[1]->GetBinContent(b), 2));
        const double ddnc   = symm(hup[2]->GetBinContent(b),  hdown[2]->GetBinContent(b));
        const double model1 = symm(hup[3]->GetBinContent(b),  hdown[3]->GetBinContent(b));
        const double model2 = symm(hup[4]->GetBinContent(b),  hdown[4]->GetBinContent(b));
        const double model3 = symm(hup[5]->GetBinContent(b),  hdown[5]->GetBinContent(b));

        // JES 1-30
        double jesSum2 = 0.0;
        for (int j = 0; j < 30; ++j) {
            double v = symm(hup[6+j]->GetBinContent(b), hdown[6+j]->GetBinContent(b));
            jesSum2 += v * v;
        }

        // JER 1-13
        double jerSum2 = 0.0;
        for (int j = 0; j < 13; ++j) {
            double v = symm(hup[36+j]->GetBinContent(b), hdown[36+j]->GetBinContent(b));
            jerSum2 += v * v;
        }

        const double prw  = symm(hup[49]->GetBinContent(b), hdown[49]->GetBinContent(b));
        const double tile = symm(hup[50]->GetBinContent(b), hdown[50]->GetBinContent(b));

        const double total = std::sqrt(stat*stat + ddnc*ddnc +
                                       model1*model1 + model2*model2 + model3*model3 +
                                       jesSum2 + jerSum2 +
                                       prw*prw + tile*tile);
        hdata->SetBinError(b, total);
    }
}

// ============================================================
//  FillTheoPrefit (TF1 k-factor variant)
// ============================================================

void FillTheoPrefit(TH1D* htheo,
                    const std::vector<std::vector<double>>& myparam,
                    double aS, TF1* /*f_myKNNLO*/, TH1D* npcorr)
{
    std::cout << "in FillTheoPrefit function" << std::endl;
    std::cout << "nbins theoprefit: " << htheo->GetNbinsX() << std::endl;

    for (int i = 0; i < htheo->GetNbinsX(); ++i) {
        const double p0  = myparam[i][0];
        const double p1  = myparam[i][1];
        const double p2  = myparam[i][2];
        const double bin = p0*std::pow(aS,1) + p1*std::pow(aS,2) + p2*std::pow(aS,3);
        // NP correction currently forced to unity; restore the line below to use it:
        // const double np  = npcorr->GetBinContent(i+1);
        const double np  = 1.0;
        htheo->SetBinContent(i + 1, np * bin);
        htheo->SetBinError  (i + 1, 0.01 * bin * np);
    }
}

// ============================================================
//  FillTheoPrefit (scalar k-factor variant)
// ============================================================

void FillTheoPrefit(TH1D* htheo,
                    const std::vector<std::vector<double>>& myparam,
                    double aS, double /*myKNNLO*/, TH1D* npcorr)
{
    for (int i = 0; i < htheo->GetNbinsX(); ++i) {
        const double p0  = myparam[i][0];
        const double p1  = myparam[i][1];
        const double p2  = myparam[i][2];
        const double bin = p0*std::pow(aS,1) + p1*std::pow(aS,2) + p2*std::pow(aS,3);
        const double np  = npcorr->GetBinContent(i + 1);
        htheo->SetBinContent(i + 1, np * bin);
        htheo->SetBinError  (i + 1, 0.01 * bin * np);
    }
}

// ============================================================
//  PDF uncertainty filling — MSHT20 eigenvector method
// ============================================================

void FillPDFUncertMSHT(TH1D* histo, TH1D* h_pdf[kNPDFMSHT],
                        unsigned int nmin, unsigned int nmax)
{
    std::cout << "in FillPDFUncertMSHT function" << std::endl;
    std::cout << "nmin: " << nmin << ", nmax: " << nmax << std::endl;

    // Auxiliary vectors for the 32 MSHT eigenvector pairs
    std::vector<std::vector<double>*> pdfAux = {
        pdf1TheoAux,  pdf2TheoAux,  pdf3TheoAux,  pdf4TheoAux,
        pdf5TheoAux,  pdf6TheoAux,  pdf7TheoAux,  pdf8TheoAux,
        pdf9TheoAux,  pdf10TheoAux, pdf11TheoAux, pdf12TheoAux,
        pdf13TheoAux, pdf14TheoAux, pdf15TheoAux, pdf16TheoAux,
        pdf17TheoAux, pdf18TheoAux, pdf19TheoAux, pdf20TheoAux,
        pdf21TheoAux, pdf22TheoAux, pdf23TheoAux, pdf24TheoAux,
        pdf25TheoAux, pdf26TheoAux, pdf27TheoAux, pdf28TheoAux,
        pdf29TheoAux, pdf30TheoAux, pdf31TheoAux, pdf32TheoAux
    };

    for (unsigned int i = nmin; i < nmax; ++i) {
        const double central = histo->GetBinContent(i + 1);
        for (unsigned int pair = 0; pair < 32; ++pair) {
            const double binPlus  = h_pdf[2*pair  ]->GetBinContent(i + 1) / central;
            const double binMinus = h_pdf[2*pair+1]->GetBinContent(i + 1) / central;
            pdfAux[pair]->push_back(0.5 * std::fabs(binPlus - binMinus));
        }
    }
}

// ============================================================
//  PDF uncertainty filling — CT18 eigenvector method
// ============================================================

void FillPDFUncertCT18(TH1D* histo, TH1D* h_pdf[kNPDFCT18],
                        unsigned int nmin, unsigned int nmax)
{
    std::cout << "in FillPDFUncertCT18 function" << std::endl;

    std::vector<std::vector<double>*> pdfAux = {
        pdf1TheoAux,  pdf2TheoAux,  pdf3TheoAux,  pdf4TheoAux,
        pdf5TheoAux,  pdf6TheoAux,  pdf7TheoAux,  pdf8TheoAux,
        pdf9TheoAux,  pdf10TheoAux, pdf11TheoAux, pdf12TheoAux,
        pdf13TheoAux, pdf14TheoAux, pdf15TheoAux, pdf16TheoAux,
        pdf17TheoAux, pdf18TheoAux, pdf19TheoAux, pdf20TheoAux,
        pdf21TheoAux, pdf22TheoAux, pdf23TheoAux, pdf24TheoAux,
        pdf25TheoAux, pdf26TheoAux, pdf27TheoAux, pdf28TheoAux,
        pdf29TheoAux
    };

    for (unsigned int i = nmin; i < nmax; ++i) {
        const double central = histo->GetBinContent(i + 1);
        for (unsigned int pair = 0; pair < 29; ++pair) {
            const double binPlus  = h_pdf[2*pair  ]->GetBinContent(i + 1) / central;
            const double binMinus = h_pdf[2*pair+1]->GetBinContent(i + 1) / central;
            pdfAux[pair]->push_back(0.5 * std::fabs(binPlus - binMinus) / kCT18Norm);
        }
    }
}

// ============================================================
//  PDF uncertainty filling — MC replica method
// ============================================================

void FillPDFUncert(TH1D* histo, TH1D* h_pdf[kNPDFRep],
                   unsigned int nmin, unsigned int nmax)
{
    std::cout << "in FillPDFUncert function" << std::endl;

    // 25 replica pairs → 25 eigenvector-equivalent uncertainties
    std::vector<std::vector<double>*> pdfAux = {
        pdf1TheoAux,  pdf2TheoAux,  pdf3TheoAux,  pdf4TheoAux,  pdf5TheoAux,
        pdf6TheoAux,  pdf7TheoAux,  pdf8TheoAux,  pdf9TheoAux,  pdf10TheoAux,
        pdf11TheoAux, pdf12TheoAux, pdf13TheoAux, pdf14TheoAux, pdf15TheoAux,
        pdf16TheoAux, pdf17TheoAux, pdf18TheoAux, pdf19TheoAux, pdf20TheoAux,
        pdf21TheoAux, pdf22TheoAux, pdf23TheoAux, pdf24TheoAux, pdf25TheoAux
    };

    for (unsigned int i = nmin; i < nmax; ++i) {
        const double central = histo->GetBinContent(i + 1);
        for (unsigned int pair = 0; pair < 25; ++pair) {
            const double rPlus  = h_pdf[2*pair  ]->GetBinContent(i + 1) / central;
            const double rMinus = h_pdf[2*pair+1]->GetBinContent(i + 1) / central;
            pdfAux[pair]->push_back(0.5 * (std::fabs(rPlus  - 1.0) +
                                           std::fabs(rMinus - 1.0)));
        }
    }
}

// ============================================================
//  getTheoUncert — pre-smoothed band from two histograms
// ============================================================

void getTheoUncert(TH1D* h_up, TH1D* h_down,
                   TH1D* h_rel_up, TH1D* h_rel_down)
{
    std::cout << "in getTheoUncert function" << std::endl;
    for (int i = 0; i < h_up->GetNbinsX(); ++i) {
        h_up  ->SetBinContent(i+1, 1.0 + h_rel_up  ->GetBinContent(i+1));
        h_down->SetBinContent(i+1, 1.0 + h_rel_down->GetBinContent(i+1));
        h_up  ->SetBinError(i+1, 0.0);
        h_down->SetBinError(i+1, 0.0);
    }
}

// ============================================================
//  getTheoUncert — envelope from 6 scale variations
// ============================================================

void getTheoUncert(TH1D* h_up, TH1D* h_down,
                   TH1D* theo, TH1D* theo_var[6], double /*knnlo*/,
                   const TString& /*region*/)
{
    std::cout << "in getTheoUncert function (6-variation envelope)" << std::endl;
    for (int i = 0; i < theo->GetNbinsX(); ++i) {
        const double bin = theo->GetBinContent(i + 1);
        std::vector<double> deltas(6);
        for (int k = 0; k < 6; ++k)
            deltas[k] = theo_var[k]->GetBinContent(i + 1) - bin;
        std::sort(deltas.begin(), deltas.end());
        h_up  ->SetBinContent(i+1, deltas.back() / bin);
        h_down->SetBinContent(i+1, deltas.front() / bin);
        h_up  ->SetBinError(i+1, 0.0);
        h_down->SetBinError(i+1, 0.0);
    }
}

// ============================================================
//  getPDFUncertMSHT — Hessian uncertainty from 64 eigenvectors
// ============================================================

void getPDFUncertMSHT(TH1D* h_up, TH1D* h_down,
                      TH1D* theo, TH1D* pdf_var[kNPDFMSHT], double /*knnlo*/)
{
    std::cout << "in getPDFUncertMSHT function" << std::endl;
    for (int i = 0; i < theo->GetNbinsX(); ++i) {
        const double bin = theo->GetBinContent(i + 1);
        double sigPlus2 = 0.0, sigMinus2 = 0.0;
        for (unsigned int pair = 0; pair < 32; ++pair) {
            const double vP = pdf_var[2*pair  ]->GetBinContent(i+1) - bin;
            const double vM = pdf_var[2*pair+1]->GetBinContent(i+1) - bin;
            if (vP > 0) sigPlus2  += vP * vP;
            if (vP < 0) sigMinus2 += vP * vP;
            if (vM > 0) sigPlus2  += vM * vM;
            if (vM < 0) sigMinus2 += vM * vM;
        }
        const double sig = 0.5 * (std::sqrt(sigPlus2) + std::sqrt(sigMinus2));
        h_up  ->SetBinContent(i+1, 1.0 + sig / bin);
        h_down->SetBinContent(i+1, 1.0 - sig / bin);
        h_up  ->SetBinError(i+1, 0.0);
        h_down->SetBinError(i+1, 0.0);
    }
}

// ============================================================
//  getPDFUncert — MC replica method (50 replicas)
// ============================================================

void getPDFUncert(TH1D* h_up, TH1D* h_down,
                  TH1D* theo, TH1D* pdf_var[kNPDFRep], double /*knnlo*/)
{
    std::cout << "in getPDFUncert function" << std::endl;
    for (int i = 0; i < theo->GetNbinsX(); ++i) {
        const double bin = theo->GetBinContent(i + 1);
        double sigPlus2 = 0.0, sigMinus2 = 0.0;
        for (unsigned int j = 0; j < kNPDFRep; ++j) {
            const double delta = pdf_var[j]->GetBinContent(i+1) - bin;
            if (delta > 0) sigPlus2  += delta * delta;
            else           sigMinus2 += delta * delta;
        }
        const double sig = 0.5 * (std::sqrt(sigPlus2) + std::sqrt(sigMinus2));
        h_up  ->SetBinContent(i+1, 1.0 + sig / bin);
        h_down->SetBinContent(i+1, 1.0 - sig / bin);
        h_up  ->SetBinError(i+1, 0.0);
        h_down->SetBinError(i+1, 0.0);
    }
}

// ============================================================
//  getPDFUncertCT18 — CT18 90% CL rescaled to 68% CL
// ============================================================

void getPDFUncertCT18(TH1D* h_up, TH1D* h_down,
                      TH1D* theo, TH1D* pdf_var[kNPDFCT18], double /*knnlo*/)
{
    std::cout << "in getPDFUncertCT18 function" << std::endl;
    for (int i = 0; i < theo->GetNbinsX(); ++i) {
        const double bin = theo->GetBinContent(i + 1);
        double sigma2 = 0.0;
        for (unsigned int pair = 0; pair < 29; ++pair) {
            const double half = 0.5 * (pdf_var[2*pair  ]->GetBinContent(i+1) -
                                       pdf_var[2*pair+1]->GetBinContent(i+1));
            sigma2 += half * half;
        }
        const double sig = std::sqrt(sigma2) / kCT18Norm;
        h_up  ->SetBinContent(i+1, 1.0 + sig / bin);
        h_down->SetBinContent(i+1, 1.0 - sig / bin);
        h_up  ->SetBinError(i+1, 0.0);
        h_down->SetBinError(i+1, 0.0);
    }
}

// ============================================================
//  getRatio
// ============================================================

void getRatio(TH1D* hratio, TH1D* htheo, TH1D* hdata)
{
    std::cout << "in getRatio function" << std::endl;
    for (int i = 0; i < hratio->GetNbinsX(); ++i) {
        const double data  = hdata->GetBinContent(i+1);
        const double dataE = hdata->GetBinError  (i+1);
        const double theo  = htheo->GetBinContent(i+1);
        hratio->SetBinContent(i+1, data  / theo);
        hratio->SetBinError  (i+1, dataE / theo);
    }
}

// ============================================================
//  chi2MinMSHT
//  Main function: load inputs, fill systematic vectors,
//  run MIGRAD, produce pre-fit and nuisance-parameter plots.
// ============================================================

std::vector<double> chi2MinMSHT(TString variation, std::string pTregion,
                                 double KNNLO, int xlow, int xhigh,
                                 std::vector<double>& NuisPar)
{
    NuisPar.clear();
    std::cout << "\nin chi2Min function" << std::endl;
    std::cout << "pTregion: " << pTregion << std::endl;

    // ---------------------------------------------------------
    //  K-factor function
    // ---------------------------------------------------------
    TFile* funcFile = new TFile("./Kfactor_func_allbins.root");
    const int ri = regionIndex(pTregion);
    char funcName[64];
    std::snprintf(funcName, sizeof(funcName), "myexp_pt%d", ri == 0 ? 0 : ri);
    TF1* kfactor_func = static_cast<TF1*>(funcFile->Get(funcName));

    // ---------------------------------------------------------
    //  Theory parameters
    // ---------------------------------------------------------
    const TString fileName = "R32_histograms_PDF_MSHT-HT.root";
    const std::string var_str(variation.Data());
    std::vector<std::vector<double>> mytheoParam =
        getParsMSHT(var_str, pTregion, fileName);

    for (const auto& p : mytheoParam)
        paramAux->push_back(p);

    std::cout << "mytheoParam.size: "    << mytheoParam.size()    << std::endl;
    std::cout << "mytheoParam[0].size: " << mytheoParam[0].size() << std::endl;
    std::cout << "mytheoParam[0]: "
              << mytheoParam[0][0] << ", " << mytheoParam[0][1] << ", "
              << mytheoParam[0][2] << ", " << mytheoParam[0][3] << std::endl;

    pTflag->push_back(pTregion);
    knnlo ->push_back(KNNLO);
    knnlo_func->push_back(kfactor_func);

    if (xhigh > static_cast<int>(mytheoParam.size()))
        xhigh = static_cast<int>(mytheoParam.size());
    const int lbin = xlow, hbin = xhigh;
    std::cout << "xlow: " << lbin << "\nxhigh: " << hbin << std::endl;
    blow ->push_back(xlow);
    bhigh->push_back(xhigh);

    // ---------------------------------------------------------
    //  Open ROOT files
    // ---------------------------------------------------------
    TFile* fileData = TFile::Open("./HEPData-105630-v1.root");
    TFile* fileTheo = TFile::Open(fileName);
    std::cout << "Opening fileData: " << fileData->GetName() << std::endl;
    std::cout << "Opening fileTheo: " << fileTheo->GetName() << std::endl;
    fileDataAux->push_back(fileData);
    fileTheoAux->push_back(fileTheo);

    // ---------------------------------------------------------
    //  Theory and data histograms (5 pT regions)
    // ---------------------------------------------------------
    std::cout << "defining theory histograms" << std::endl;
    TH1D* histoTheo[5], *histoData[5], *npTheo[5];
    for (int i = 0; i < 5; ++i) {
        char name[128];
        std::snprintf(name, sizeof(name), "h_R32_pt%d_muHTp_msht_as118", i);
        histoTheo[i] = static_cast<TH1D*>(fileTheo->Get(name));

        std::snprintf(name, sizeof(name), "Binned_R32__pt3_%d/Hist1D_y1", i);
        histoData[i] = static_cast<TH1D*>(fileData->Get(name));

        std::snprintf(name, sizeof(name), "h_np_pt%d", i);
        npTheo[i] = static_cast<TH1D*>(fileTheo->Get(name));
    }

    // ---------------------------------------------------------
    //  Systematic variation histograms (56 UP/DOWN per pT region)
    // ---------------------------------------------------------
    std::cout << "Defining UP/DOWN systematic variations" << std::endl;
    TH1D* sysUP  [5][kNSys];
    TH1D* sysDOWN[5][kNSys];

    for (unsigned int j = 0; j < kNSys; ++j) {
        char nUP[64], nDOWN[64];
        std::snprintf(nUP,   sizeof(nUP),   "Hist1D_y1_e%dplus",  j+1);
        std::snprintf(nDOWN, sizeof(nDOWN), "Hist1D_y1_e%dminus", j+1);
        for (int r = 0; r < 5; ++r) {
            char path[128];
            std::snprintf(path, sizeof(path), "Binned_R32__pt3_%d/%s", r, nUP);
            sysUP  [r][j] = static_cast<TH1D*>(fileData->Get(path));
            std::snprintf(path, sizeof(path), "Binned_R32__pt3_%d/%s", r, nDOWN);
            sysDOWN[r][j] = static_cast<TH1D*>(fileData->Get(path));
        }
    }

    // ---------------------------------------------------------
    //  Fill auxiliary vectors for the selected pT region
    // ---------------------------------------------------------
    std::cout << "Filling data and syst histograms" << std::endl;
    FillPoints   (histoData[ri]);
    FillStatTheo (histoTheo[ri], 0, histoData[ri]->GetNbinsX());
    FillData     (histoData[ri], 0, histoData[ri]->GetNbinsX());
    FillVectSyst (sysUP[ri], sysDOWN[ri], 0, histoData[ri]->GetNbinsX());
    FillNP       (npTheo[ri]);

    // ---------------------------------------------------------
    //  Scale uncertainty histograms
    // ---------------------------------------------------------
    TH1D* h_theo_up[5], *h_theo_dn[5];
    for (int i = 0; i < 5; ++i) {
        char name[128];
        std::snprintf(name, sizeof(name), "h_R32_scale_uncert_smoothed_up_pt%d",   i);
        h_theo_up[i] = static_cast<TH1D*>(funcFile->Get(name));
        std::snprintf(name, sizeof(name), "h_R32_scale_uncert_smoothed_down_pt%d", i);
        h_theo_dn[i] = static_cast<TH1D*>(funcFile->Get(name));
    }

    // ---------------------------------------------------------
    //  PDF uncertainty histograms (64 MSHT eigenvectors)
    // ---------------------------------------------------------
    TH1D* h_pdf[5][kNPDFMSHT];
    for (unsigned int j = 0; j < kNPDFMSHT; ++j) {
        for (int r = 0; r < 5; ++r) {
            char name[128];
            std::snprintf(name, sizeof(name),
                          "h_R32_pt%d_muHTp_msht_pdf%d", r, j+1);
            h_pdf[r][j] = static_cast<TH1D*>(fileTheo->Get(name));
        }
    }
    FillPDFUncertMSHT(histoTheo[ri], h_pdf[ri], 0, histoData[ri]->GetNbinsX());

    // ---------------------------------------------------------
    //  Pre-fit theory clones (central, αS up, αS down)
    // ---------------------------------------------------------
    TH1D* theo_prefit      = static_cast<TH1D*>(histoData[ri]->Clone("theo_prefit"));
    TH1D* theo_prefit_asUP = static_cast<TH1D*>(histoData[ri]->Clone("theo_prefit_asUP"));
    TH1D* theo_prefit_asDN = static_cast<TH1D*>(histoData[ri]->Clone("theo_prefit_asDOWN"));
    TH1D* data_prefit      = static_cast<TH1D*>(histoData[ri]->Clone("data_prefit"));

    FillTheoPrefit(theo_prefit,      mytheoParam, 0.118, kfactor_func, npTheo[ri]);
    FillTheoPrefit(theo_prefit_asUP, mytheoParam, 0.128, kfactor_func, npTheo[ri]);
    FillTheoPrefit(theo_prefit_asDN, mytheoParam, 0.108, kfactor_func, npTheo[ri]);
    SetStatSysData(data_prefit, sysUP[ri], sysDOWN[ri]);

    std::cout << "theo_prefit nbins:      " << theo_prefit     ->GetNbinsX() << std::endl;
    std::cout << "theo_prefit asUP nbins: " << theo_prefit_asUP->GetNbinsX() << std::endl;
    std::cout << "theo_prefit asDN nbins: " << theo_prefit_asDN->GetNbinsX() << std::endl;

    // ---------------------------------------------------------
    //  Uncertainty band histograms
    // ---------------------------------------------------------
    TH1D* theo_uncert_up   = static_cast<TH1D*>(theo_prefit->Clone("theo_uncert_up"));
    TH1D* theo_uncert_down = static_cast<TH1D*>(theo_prefit->Clone("theo_uncert_down"));
    TH1D* pdf_uncert_up    = static_cast<TH1D*>(theo_prefit->Clone("pdf_uncert_up"));
    TH1D* pdf_uncert_down  = static_cast<TH1D*>(theo_prefit->Clone("pdf_uncert_down"));

    getTheoUncert  (theo_uncert_up, theo_uncert_down, h_theo_up[ri], h_theo_dn[ri]);
    getPDFUncertMSHT(pdf_uncert_up, pdf_uncert_down,  histoTheo[ri], h_pdf[ri], KNNLO);

    for (int i = 0; i < theo_uncert_up->GetNbinsX(); ++i) {
        scalesUPTheoAux  ->push_back(theo_uncert_up  ->GetBinContent(i+1));
        scalesDOWNTheoAux->push_back(theo_uncert_down->GetBinContent(i+1));
    }

    // ---------------------------------------------------------
    //  TMinuit fit
    // ---------------------------------------------------------
    TMinuit minuit(kNPar);
    minuit.SetFCN(fcn_MSHT);
    minuit.SetGraphicsMode(false);
    minuit.SetPrintLevel(-1);

    // Parameter 0: alpha_S
    minuit.DefineParameter(0, " aS ", 0.1153, 1e-5, 0.0, 0.0);

    // Parameters 1..kNPar-1: nuisance parameters lambda_k
    for (int k = 1; k < kNPar; ++k) {
        const std::string name = "lambda" + std::to_string(k);
        minuit.DefineParameter(k, name.c_str(), 0.0, 1e-5, 0.0, 0.0);
    }

    minuit.SetErrorDef(1.0);
    minuit.SetMaxIterations(100000000);

    std::cout << "Performing minimisation" << std::endl;
    minuit.Migrad();

    double outpar[kNPar], errpar[kNPar];
    for (int i = 0; i < kNPar; ++i)
        minuit.GetParameter(i, outpar[i], errpar[i]);

    const double alphaValue = outpar[0];
    const double alphaError = errpar[0];
    const double chi2       = minuit.fAmin;

    // ---------------------------------------------------------
    //  Nuisance parameter plot
    // ---------------------------------------------------------
    const TString nuisNames[kNPar-1] = {
        "DDNC",
        "Model_Sherpa_Eff", "Model_Sherpa_Resp", "Model_Sherpa_Fake",
        "JES_EffectiveNP_Detector1", "JES_EffectiveNP_Detector2",
        "JES_EffectiveNP_Mixed1",    "JES_EffectiveNP_Mixed2",   "JES_EffectiveNP_Mixed3",
        "JES_EffectiveNP_Model1",    "JES_EffectiveNP_Model2",
        "JES_EffectiveNP_Model3",    "JES_EffectiveNP_Model4",
        "JES_EffectiveNP_Stat1",     "JES_EffectiveNP_Stat2",
        "JES_EffectiveNP_Stat3",     "JES_EffectiveNP_Stat4",
        "JES_EffectiveNP_Stat5",     "JES_EffectiveNP_Stat6",
        "JES_EtaIntercalibration_Model",
        "JES_EtaIntercalibration_NonClosure_2018data",
        "JES_EtaIntercalibration_NonClosure_highE",
        "JES_EtaIntercalibration_NonClosure_negEta",
        "JES_EtaIntercalibration_NonClosure_posEta",
        "JES_EtaIntercalibration_TotalStat",
        "JES_Flavor_Composition_prop",  "JES_Flavor_Response_prop",
        "JES_Flavor_PerJet_GenShower",  "JES_Flavor_PerJet_GenShower_HF",
        "JES_Flavor_PerJet_Hadronisation", "JES_Flavor_PerJet_Hadronisation_HF",
        "JES_Flavor_PerJet_Shower",     "JES_Flavor_PerJet_Shower_HF",
        "JES_Pileup_OffsetMu",  "JES_Pileup_OffsetNPV",
        "JES_Pileup_PtTerm",    "JES_Pileup_RhoTopology",
        "JES_PunchThrough_MC16","JES_SingleParticle_HighPt",
        "JER_DataVsMC_MC16",
        "JER_EffectiveNP1",  "JER_EffectiveNP10", "JER_EffectiveNP11",
        "JER_EffectiveNP12restTerm",
        "JER_EffectiveNP2",  "JER_EffectiveNP3",  "JER_EffectiveNP4",
        "JER_EffectiveNP5",  "JER_EffectiveNP6",
        "JER_EffectiveNP7_1","JER_EffectiveNP8",  "JER_EffectiveNP7_2",
        "PURW", "Tile",
        "PDF_eigen1",  "PDF_eigen2",  "PDF_eigen3",  "PDF_eigen4",
        "PDF_eigen5",  "PDF_eigen6",  "PDF_eigen7",  "PDF_eigen8",
        "PDF_eigen9",  "PDF_eigen10", "PDF_eigen11", "PDF_eigen12",
        "PDF_eigen13", "PDF_eigen14", "PDF_eigen15", "PDF_eigen16",
        "PDF_eigen17", "PDF_eigen18", "PDF_eigen19", "PDF_eigen20",
        "PDF_eigen21", "PDF_eigen22", "PDF_eigen23", "PDF_eigen24",
        "PDF_eigen25", "PDF_eigen26", "PDF_eigen27", "PDF_eigen28",
        "PDF_eigen29", "PDF_eigen30", "PDF_eigen31", "PDF_eigen32"
    };

    TH1D* nuisHist   = new TH1D("nuisParameters", "nuisParameters", kNPar-1, 0, kNPar-1);
    TH1D* sigPlus1   = new TH1D("sigmaPlus1",  "", kNPar-1, 0, kNPar-1);
    TH1D* sigMinus1  = new TH1D("sigmaMinus1", "", kNPar-1, 0, kNPar-1);
    TH1D* sigPlus2   = new TH1D("sigmaPlus2",  "", kNPar-1, 0, kNPar-1);
    TH1D* sigMinus2  = new TH1D("sigmaMinus2", "", kNPar-1, 0, kNPar-1);

    for (int k = 1; k < kNPar; ++k) {
        nuisHist ->SetBinContent(k, outpar[k]);
        nuisHist ->SetBinError  (k, errpar[k]);
        sigPlus1 ->SetBinContent(k,  1.0);
        sigMinus1->SetBinContent(k, -1.0);
        sigPlus2 ->SetBinContent(k,  2.0);
        sigMinus2->SetBinContent(k, -2.0);
        NuisPar.push_back(outpar[k]);
        nuisHist->GetXaxis()->SetBinLabel(k, nuisNames[k-1]);
    }

    // Style sigma bands
    for (auto* h : {sigPlus1, sigMinus1}) {
        h->SetFillColor(kGreen); h->SetLineColor(kGreen); h->SetLineWidth(0);
    }
    for (auto* h : {sigPlus2, sigMinus2}) {
        h->SetFillColor(kYellow); h->SetLineColor(kYellow); h->SetLineWidth(0);
    }

    TCanvas* cNuis = new TCanvas("cNuis", "cNuis", 1150, 650);
    cNuis->cd();
    gPad->SetBottomMargin(0.33);
    nuisHist->GetXaxis()->SetLabelSize(0.025);
    nuisHist->GetXaxis()->SetTitleSize(0);
    nuisHist->GetYaxis()->SetTitle("Fit Result");
    nuisHist->GetYaxis()->SetRangeUser(-7.5, 7.5);
    nuisHist->Draw("");
    sigPlus2 ->Draw("same"); sigMinus2->Draw("same");
    sigPlus1 ->Draw("same"); sigMinus1->Draw("same");
    nuisHist ->Draw("same");

    TLegend* leg1 = new TLegend(0.60, 0.75, 0.85, 0.90);
    leg1->AddEntry(nuisHist,  "NP Values");
    leg1->AddEntry(sigPlus1,  "#pm 1 #sigma contour", "F");
    leg1->AddEntry(sigPlus2,  "#pm 2 #sigma contour", "F");
    leg1->SetLineColor(0); leg1->SetLineWidth(0); leg1->SetFillColor(0);
    leg1->Draw();

    ATLAS_INTERNAL_LABEL(0.20, 0.85, 1);
    myText(0.20, 0.77, 1, regionLabel(pTregion), 0.05);
    myText(0.20, 0.45, 1, "MSHT20 PDFs", 0.04);
    nuisHist->Draw("sameaxis");

    const TString nuisStem = "nuisParameters_" + TString(pTregion) + "_" + variation + "_MSHT";
    cNuis->SaveAs(nuisStem + ".pdf");
    cNuis->SaveAs(nuisStem + ".png");

    // ---------------------------------------------------------
    //  Pre-fit canvas
    // ---------------------------------------------------------
    TCanvas* cPre = new TCanvas("cPre", "cPre", 650, 650);
    cPre->SetFillStyle(4000);
    cPre->cd();

    TPad* cPre_1 = new TPad("cPre_1", "", 0.0, 0.23, 0.95, 0.99);
    cPre_1->SetFillStyle(4000); cPre_1->Draw();
    TPad* cPre_2 = new TPad("cPre_2", "", 0.0, 0.02, 0.95, 0.36);
    cPre_2->SetFillStyle(4000); cPre_2->Draw();

    cPre_1->cd();
    gPad->SetLogx();
    auto [ylo, yhi] = prefitYRange(pTregion);
    theo_prefit->GetYaxis()->SetTitle("R_{32}");
    theo_prefit->GetYaxis()->SetTitleSize(0.10);
    theo_prefit->GetYaxis()->SetTitleOffset(0.5);
    theo_prefit->GetYaxis()->SetRangeUser(ylo, yhi);
    theo_prefit->GetXaxis()->SetTitle("H_{T2} (GeV)");
    theo_prefit->GetXaxis()->SetLabelSize(0);
    theo_prefit->GetXaxis()->SetTitleSize(0);

    data_prefit->SetMarkerStyle(20);
    theo_prefit     ->SetLineColor(2);
    theo_prefit_asUP->SetLineColor(3);
    theo_prefit_asDN->SetLineColor(4);

    const double xupper = theo_prefit->GetXaxis()->GetBinUpEdge(hbin);
    theo_prefit->GetXaxis()->SetRangeUser(250., 0.9 * xupper);
    theo_prefit->Draw("histo");
    theo_prefit_asUP->Draw("histosame");
    theo_prefit_asDN->Draw("histosame");
    data_prefit->Draw("pex0same");

    myMarkerText(0.60, 0.89, 1, 20, "Data (stat #oplus sys)", 1.2, 0.04, 0.04);
    myText(0.60, 0.80, 1, "NNLO #otimes NP", 0.035);
    myLine(0.60, 0.75, 0.03, 4, 1, "#alpha_{S}(M_{Z}) = 0.108 - Min. Value", 0.03, 0.03);
    myLine(0.60, 0.70, 0.03, 2, 1, "#alpha_{S}(M_{Z}) = 0.118",              0.03, 0.03);
    myLine(0.60, 0.65, 0.03, 3, 1, "#alpha_{S}(M_{Z}) = 0.128 - Max. Value", 0.03, 0.03);
    myText(0.20, 0.76, 1, "MSHT20 PDFs",       0.04);
    myText(0.20, 0.70, 1, regionLabel(pTregion), 0.05);
    ATLAS_INTERNAL_LABEL(0.20, 0.85, 1);

    cPre_2->cd();
    gPad->SetLogx();
    gPad->SetBottomMargin(0.33);

    TH1D* ratio_pre = static_cast<TH1D*>(data_prefit->Clone("ratio_pre"));
    getRatio(ratio_pre, theo_prefit, data_prefit);
    ratio_pre->GetYaxis()->SetRangeUser(0.82, 1.23);
    ratio_pre->GetYaxis()->SetLabelSize(0.10);
    ratio_pre->GetYaxis()->SetTitleSize(0.10);
    ratio_pre->GetYaxis()->SetTitleOffset(0.5);
    ratio_pre->GetYaxis()->SetTitle("data/theory");
    ratio_pre->GetXaxis()->SetTitle("H_{T2} [GeV]");
    ratio_pre->GetXaxis()->SetTitleOffset(0.9);
    ratio_pre->GetXaxis()->SetLabelSize(0.10);
    ratio_pre->GetXaxis()->SetTitleSize(0.16);
    ratio_pre->GetXaxis()->SetMoreLogLabels();
    ratio_pre->GetXaxis()->SetNoExponent();
    ratio_pre->SetMarkerStyle(20);
    ratio_pre->SetMarkerColor(1);
    ratio_pre->SetLineColor(1);
    ratio_pre->GetXaxis()->SetRangeUser(250., xupper);
    ratio_pre->Draw("pex0");

    theo_uncert_up->SetLineColor(2); theo_uncert_up->SetFillColor(2);
    theo_uncert_up->SetFillStyle(3001);
    theo_uncert_up->Draw("histosame");
    theo_uncert_down->SetFillColor(10);
    theo_uncert_down->Draw("histosame");

    pdf_uncert_up->SetLineColor(4); pdf_uncert_up->SetFillColor(4);
    pdf_uncert_up->SetFillStyle(3001);
    pdf_uncert_up->Draw("histosame");
    pdf_uncert_down->SetFillColor(2); pdf_uncert_down->SetFillStyle(3001);
    pdf_uncert_down->Draw("histosame");
    theo_uncert_down->Draw("histosame");
    ratio_pre->Draw("pex0same");

    myBoxText(0.27, 0.85, 0.05, 2, 3001, "Theo. uncertainty (#mu_{R},#mu_{F})", 0.08, 0.08);
    myBoxText(0.65, 0.85, 0.05, 4, 3001, "Theo. uncertainty (PDFs)",            0.08, 0.08);
    gPad->RedrawAxis();

    TLine* ll_pre = new TLine(250., 1.0, xupper, 1.0);
    ll_pre->Draw("same");

    const TString preStem = (KNNLO != 1.0 ? "prefit_" : "prefit_noKnnlo_")
                            + TString(pTregion) + "_" + variation + "_MSHT";
    cPre->SaveAs(preStem + ".pdf");
    cPre->SaveAs(preStem + ".png");

    // ---------------------------------------------------------
    //  Correlation matrix
    // ---------------------------------------------------------
    const TString allNames[kNPar] = {
        "aS",
        "DDNC", "Model_ShLund", "Model_ShAHADIC",
        "JES_EffectiveNP_Detector1", "JES_EffectiveNP_Detector2",
        "JES_EffectiveNP_Mixed1",    "JES_EffectiveNP_Mixed2",   "JES_EffectiveNP_Mixed3",
        "JES_EffectiveNP_Model1",    "JES_EffectiveNP_Model2",
        "JES_EffectiveNP_Model3",    "JES_EffectiveNP_Model4",
        "JES_EffectiveNP_Stat1",     "JES_EffectiveNP_Stat2",
        "JES_EffectiveNP_Stat3",     "JES_EffectiveNP_Stat4",
        "JES_EffectiveNP_Stat5",     "JES_EffectiveNP_Stat6",
        "JES_EtaIntercalibration_Model",
        "JES_EtaIntercalibration_NonClosure_2018data",
        "JES_EtaIntercalibration_NonClosure_highE",
        "JES_EtaIntercalibration_NonClosure_negEta",
        "JES_EtaIntercalibration_NonClosure_posEta",
        "JES_EtaIntercalibration_TotalStat",
        "JES_Flavor_Composition_prop",  "JES_Flavor_Response_prop",
        "JES_Flavor_PerJet_GenShower",  "JES_Flavor_PerJet_GenShower_HT",
        "JES_Flavor_PerJet_Hadronisation", "JES_Flavor_PerJet_Hadronisation_HF",
        "JES_Flavor_PerJet_Shower",     "JES_Flavor_PerJet_Shower_HT",
        "JES_Pileup_OffsetMu",  "JES_Pileup_OffsetNPV",
        "JES_Pileup_PtTerm",    "JES_Pileup_RhoTopology",
        "JES_PunchThrough_MC16","JES_SingleParticle_HighPt",
        "JER_DataVsMC_MC16",
        "JER_EffectiveNP1",  "JER_EffectiveNP10", "JER_EffectiveNP11",
        "JER_EffectiveNP12restTerm",
        "JER_EffectiveNP2",  "JER_EffectiveNP3",  "JER_EffectiveNP4",
        "JER_EffectiveNP5",  "JER_EffectiveNP6",
        "JER_EffectiveNP7_1","JER_EffectiveNP8",  "JER_EffectiveNP7_2",
        "PURW", "Tile",
        "PDF_eigen1",  "PDF_eigen2",  "PDF_eigen3",  "PDF_eigen4",
        "PDF_eigen5",  "PDF_eigen6",  "PDF_eigen7",  "PDF_eigen8",
        "PDF_eigen9",  "PDF_eigen10", "PDF_eigen11", "PDF_eigen12",
        "PDF_eigen13", "PDF_eigen14", "PDF_eigen15", "PDF_eigen16",
        "PDF_eigen17", "PDF_eigen18", "PDF_eigen19", "PDF_eigen20",
        "PDF_eigen21", "PDF_eigen22", "PDF_eigen23", "PDF_eigen24",
        "PDF_eigen25", "PDF_eigen26", "PDF_eigen27", "PDF_eigen28",
        "PDF_eigen29", "PDF_eigen30", "PDF_eigen31", "PDF_eigen32"
    };

    Double_t mat[kNPar][kNPar];
    minuit.mnemat(&mat[0][0], kNPar);

    std::vector<double> sigmas(kNPar);
    for (int k = 0; k < kNPar; ++k)
        sigmas[k] = std::sqrt(mat[k][k]);

    TH2D* covMatrix = new TH2D("covMatrix", "covMatrix",
                               kNPar, 0, kNPar, kNPar, 0, kNPar);
    for (int i = 0; i < kNPar; ++i) {
        covMatrix->GetXaxis()->SetBinLabel(i+1, allNames[i]);
        covMatrix->GetYaxis()->SetBinLabel(i+1, allNames[i]);
        for (int j = 0; j < kNPar; ++j)
            covMatrix->SetBinContent(i+1, j+1, mat[i][j] / (sigmas[i] * sigmas[j]));
    }

    TCanvas* c2 = new TCanvas("c2", "c2", 1250, 950);
    c2->cd();
    gPad->SetTopMargin(0.02);
    gPad->SetRightMargin(0.09);
    gPad->SetBottomMargin(0.22);
    gPad->SetLeftMargin(0.135);
    covMatrix->GetXaxis()->SetLabelSize(0.015);
    covMatrix->GetYaxis()->SetLabelSize(0.015);
    covMatrix->GetZaxis()->SetLabelSize(0.030);
    covMatrix->GetZaxis()->SetRangeUser(-1, 1);
    covMatrix->Draw("colz");

    const TString corrStem = "corrMatrix_" + TString(pTregion) + "_" + variation + "_MSHT";
    c2->SaveAs(corrStem + ".pdf");
    c2->SaveAs(corrStem + ".png");

    // ---------------------------------------------------------
    //  Clean up ROOT files and all auxiliary vectors
    // ---------------------------------------------------------
    fileData->Close();
    fileTheo->Close();
    fileDataAux->clear();
    fileTheoAux->clear();
    blow->clear(); bhigh->clear();
    pTflag->clear(); knnlo->clear(); knnlo_func->clear();
    paramAux->clear(); npAux->clear();
    PointsAux->clear(); dataPointsAux->clear();
    statDataAux->clear(); statTheoAux->clear();
    ddncDataAux->clear();
    model1DataAux->clear(); model2DataAux->clear(); model3DataAux->clear();

    for (auto* v : {jes1DataAux,  jes2DataAux,  jes3DataAux,  jes4DataAux,  jes5DataAux,
                    jes6DataAux,  jes7DataAux,  jes8DataAux,  jes9DataAux,  jes10DataAux,
                    jes11DataAux, jes12DataAux, jes13DataAux, jes14DataAux, jes15DataAux,
                    jes16DataAux, jes17DataAux, jes18DataAux, jes19DataAux, jes20DataAux,
                    jes21DataAux, jes22DataAux, jes23DataAux, jes24DataAux, jes25DataAux,
                    jes26DataAux, jes27DataAux, jes28DataAux, jes29DataAux, jes30DataAux,
                    jes31DataAux, jes32DataAux, jes33DataAux, jes34DataAux, jes35DataAux})
        v->clear();

    for (auto* v : {jer1DataAux,  jer2DataAux,  jer3DataAux,  jer4DataAux,  jer5DataAux,
                    jer6DataAux,  jer7DataAux,  jer8DataAux,  jer9DataAux,  jer10DataAux,
                    jer11DataAux, jer12DataAux, jer13DataAux})
        v->clear();

    prwDataAux->clear(); tileDataAux->clear();
    scalesUPTheoAux->clear(); scalesDOWNTheoAux->clear();

    for (auto* v : {pdf1TheoAux,  pdf2TheoAux,  pdf3TheoAux,  pdf4TheoAux,
                    pdf5TheoAux,  pdf6TheoAux,  pdf7TheoAux,  pdf8TheoAux,
                    pdf9TheoAux,  pdf10TheoAux, pdf11TheoAux, pdf12TheoAux,
                    pdf13TheoAux, pdf14TheoAux, pdf15TheoAux, pdf16TheoAux,
                    pdf17TheoAux, pdf18TheoAux, pdf19TheoAux, pdf20TheoAux,
                    pdf21TheoAux, pdf22TheoAux, pdf23TheoAux, pdf24TheoAux,
                    pdf25TheoAux, pdf26TheoAux, pdf27TheoAux, pdf28TheoAux,
                    pdf29TheoAux, pdf30TheoAux, pdf31TheoAux, pdf32TheoAux,
                    pdf33TheoAux, pdf34TheoAux, pdf35TheoAux, pdf36TheoAux,
                    pdf37TheoAux, pdf38TheoAux, pdf39TheoAux, pdf40TheoAux,
                    pdf41TheoAux, pdf42TheoAux, pdf43TheoAux, pdf44TheoAux,
                    pdf45TheoAux, pdf46TheoAux, pdf47TheoAux, pdf48TheoAux,
                    pdf49TheoAux, pdf50TheoAux})
        v->clear();

    return {alphaValue, alphaError, chi2};
}
