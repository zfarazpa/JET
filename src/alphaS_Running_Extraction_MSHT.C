/// @file theFitter_MSHT_Running.C
/// @brief alpha_S extraction via chi2 minimisation using MSHT PDFs.


#include "chi2Min.C"

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Result returned by a single chi2 minimisation.
struct FitResult {
  double alphaS   = 0.0;
  double sysError = 0.0;
  double chi2     = 0.0;
};

/// Kinematic range for one Q bin.
struct QBinConfig {
  int    nmin = 0;
  int    nmax = 0;
  double KNNLO = 1.0;
};

/// How scale variations are handled.
enum class ScaleVariationMode {
  Separate,   ///< Run each variation independently (doScaleVarSep = true)
  Zero,       ///< Fill theory shifts with zeros (neither flag set)
  // Combined mode could be added here in the future
};


// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Scale labels used for theory-uncertainty variations.
static const std::vector<std::string> kScaleNames = {
  "muRU_muFU",
  "muR0_muFU",
  "muRU_muF0",
  "muRL_muFL",
  "muR0_muFL",
  "muRL_muF0",
};

static constexpr int kNQ = 6;  ///< Number of Q bins.


// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Unpack the raw vector returned by chi2MinMSHT into a FitResult.
/// Returns a zeroed FitResult and prints an error if the vector is too short.
FitResult unpackFitResult(const std::vector<double>& raw,
                          const std::string& label)
{
  if (raw.size() < 3) {
    std::cerr << "[ERROR] chi2MinMSHT returned " << raw.size()
              << " values (expected >= 3) for \"" << label
              << "\". Returning zeroed result.\n";
    return {};
  }
  return {raw[0], raw[1], raw[2]};
}

/// Print a FitResult to stdout with a descriptive label.
void printFitResult(const std::string& label, const FitResult& r)
{
  std::cout << "[FitResult] " << label
            << "  alphaS=" << r.alphaS
            << "  sysError=" << r.sysError
            << "  chi2=" << r.chi2
            << '\n';
}

/// Print a section banner to stdout.
void printBanner(const std::string& text)
{
  const std::string line(40, '/');
  std::cout << '\n' << line << '\n' << text << '\n' << line << '\n';
}


// ---------------------------------------------------------------------------
// Default Q-bin configuration
// ---------------------------------------------------------------------------

/// Build the default array of QBinConfig objects (before region overrides).
std::vector<QBinConfig> buildDefaultBins()
{
  std::vector<QBinConfig> bins(kNQ);
  bins[0] = {  0,  5, 1.0 };
  bins[1] = {  5, 10, 1.0 };
  bins[2] = { 10, 15, 1.0 };
  bins[3] = { 15, 20, 1.0 };
  bins[4] = { 20, 25, 1.0 };
  bins[5] = { 25, 30, 1.0 };
  return bins;
}

/// Apply region-specific overrides to the bin configuration.
void applyRegionOverrides(const std::string& region,
                          std::vector<QBinConfig>& bins)
{
  if (region == "pt3") {
    bins[5].nmin = 20;
    bins[5].nmax = 25;
  } else if (region == "pt4") {
    bins[4].nmin = 20;  bins[4].nmax = 24;
    bins[5].nmin = 20;  bins[5].nmax = 24;
  }
  // Add further regions here without touching the rest of the code.
}


// ---------------------------------------------------------------------------
// Fitting
// ---------------------------------------------------------------------------

/// Run the nominal fit for one Q bin and return its FitResult.
/// @param[out] nuisanceParameters  Populated by chi2MinMSHT.
FitResult runNominalFit(const std::string& region,
                        const QBinConfig& bin,
                        std::vector<double>& nuisanceParameters)
{
  std::vector<double> raw = chi2MinMSHT(
      "muR0_muF0", region, bin.KNNLO,
      bin.nmin, bin.nmax, nuisanceParameters);

  return unpackFitResult(raw, "muR0_muF0 (nominal)");
}

/// Run all six scale variations and return the sorted vector of
/// (variation alphaS - nominal alphaS) shifts.
std::vector<double> computeTheoryShifts(const std::string& region,
                                        const QBinConfig& bin,
                                        std::vector<double>& nuisanceParameters,
                                        const FitResult& nominal)
{
  std::vector<double> shifts;
  shifts.reserve(kScaleNames.size());

  for (const auto& scaleName : kScaleNames) {
    std::vector<double> raw = chi2MinMSHT(
        scaleName, region, bin.KNNLO,
        bin.nmin, bin.nmax, nuisanceParameters);

    FitResult r = unpackFitResult(raw, scaleName);
    printFitResult(scaleName, r);
    shifts.push_back(r.alphaS - nominal.alphaS);
  }

  std::sort(shifts.begin(), shifts.end());
  return shifts;
}


// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

/// Write alpha_S central value and uncertainties to a ROOT file.
void writeAlphaSOutput(const std::string& filename,
                       const FitResult& nominal,
                       const std::vector<double>& theoryShifts)
{
  if (theoryShifts.empty()) {
    std::cerr << "[ERROR] theoryShifts is empty; skipping output for "
              << filename << '\n';
    return;
  }

  TFile output(filename.c_str(), "RECREATE");
  if (!output.IsOpen()) {
    std::cerr << "[ERROR] Could not open output file: " << filename << '\n';
    return;
  }

  auto makeHist = [](const char* name, double content) {
    TH1D h(name, "", 1, 0.0, 1.0);
    h.SetBinContent(1, content);
    h.SetBinError(1, 0.0);
    return h;
  };

  makeHist("h_as",        nominal.alphaS        ).Write();
  makeHist("h_as_up",     nominal.sysError      ).Write();
  makeHist("h_as_down",  -nominal.sysError      ).Write();
  makeHist("h_as_th_up",  theoryShifts.back()   ).Write();
  makeHist("h_as_th_down",theoryShifts.front()  ).Write();

  output.Close();
  std::cout << "[Output] Wrote " << filename << '\n';
}


// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

void theFitter_MSHT_Running(const std::string& region)
{
  // ── Configuration ──────────────────────────────────────────────────────
  const ScaleVariationMode scaleMode = ScaleVariationMode::Separate;

  std::vector<QBinConfig> bins = buildDefaultBins();
  applyRegionOverrides(region, bins);

  // ── Nominal fits ────────────────────────────────────────────────────────
  printBanner("Nominal Fits");

  std::vector<FitResult>           nominalResults(kNQ);
  std::vector<std::vector<double>> nuisanceParameters(kNQ);

  for (int iQ = 0; iQ < kNQ; ++iQ) {
    std::string label = "Q" + std::to_string(iQ);
    std::cout << "\n[Fit] Running nominal fit for " << label << '\n';

    nominalResults[iQ] = runNominalFit(region, bins[iQ], nuisanceParameters[iQ]);
    printFitResult(label + " nominal", nominalResults[iQ]);
  }

  // ── Theory (scale-variation) shifts ─────────────────────────────────────
  std::vector<std::vector<double>> theoryShifts(kNQ);

  if (scaleMode == ScaleVariationMode::Separate) {
    printBanner("Scale Variations");

    for (int iQ = 0; iQ < kNQ; ++iQ) {
      std::string label = "Q" + std::to_string(iQ);
      printBanner(label + " bin");
      printFitResult("muR0_muF0", nominalResults[iQ]);

      theoryShifts[iQ] = computeTheoryShifts(
          region, bins[iQ], nuisanceParameters[iQ], nominalResults[iQ]);
    }
  } else {
    // Zero out theory shifts for ScaleVariationMode::Zero
    for (int iQ = 0; iQ < kNQ; ++iQ)
      theoryShifts[iQ].assign(kScaleNames.size(), 0.0);
  }

  // ── Summary ─────────────────────────────────────────────────────────────
  printBanner("Theory Error Summary");

  for (int iQ = 0; iQ < kNQ; ++iQ) {
    std::string label = "Q" + std::to_string(iQ);
    printFitResult(label, nominalResults[iQ]);

    const auto& shifts = theoryShifts[iQ];
    if (!shifts.empty()) {
      std::cout << "  Theory uncertainty " << label
                << ": [" << shifts.front()
                << ", " << shifts.back() << "]\n";
    }
  }

  // ── Write ROOT output ────────────────────────────────────────────────────
  printBanner("Writing Output Files");

  for (int iQ = 0; iQ < kNQ; ++iQ) {
    std::string filename =
        "aS_mZ_global_Fit_Q" + std::to_string(iQ)
        + "_" + region + ".root";

    writeAlphaSOutput(filename, nominalResults[iQ], theoryShifts[iQ]);
  }
}
