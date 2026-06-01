//  Collection of ROOT-based utilities for high-energy physics analyses.
//  Includes functions for graph arithmetic, error propagation,
//  uncertainty-band generation, histogram conversion, and
//  ATLAS-style plot formatting.
//
//  Designed to support consistent visualization of experimental
//  data and theoretical predictions within the ATLAS framework.

#ifndef ATLAS_UTILS_H
#define ATLAS_UTILS_H

#include "Rtypes.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"

// ============================================================
//  Label & annotation helpers
// ============================================================

/// @brief Draw the official ATLAS label at normalized pad coordinates.
///
/// Renders "ATLAS" in the standard bold italic typeface.  Coordinates
/// are in NDC (Normalized Device Coordinates), i.e. [0, 1] × [0, 1].
///
/// @param x      Horizontal position in NDC.
/// @param y      Vertical position in NDC.
/// @param color  ROOT color index (default: kBlack).
void ATLAS_LABEL(Double_t x, Double_t y, Color_t color = kBlack);

/// @brief Draw plain text at normalized pad coordinates.
///
/// @param x      Horizontal position in NDC.
/// @param y      Vertical position in NDC.
/// @param color  ROOT color index for the text.
/// @param text   Null-terminated string to render.
void myText(Double_t x, Double_t y, Color_t color, const char* text);

/// @brief Draw a legend entry consisting of a filled box and a text label.
///
/// @param x        Horizontal position in NDC.
/// @param y        Vertical position in NDC.
/// @param boxsize  Size of the colored box in NDC units.
/// @param mcolor   ROOT color index used to fill the box.
/// @param text     Null-terminated label string.
void myBoxText(Double_t x, Double_t y, Double_t boxsize,
               Int_t mcolor, const char* text);

/// @brief Draw a legend entry consisting of a marker symbol and a text label.
///
/// @param x      Horizontal position in NDC.
/// @param y      Vertical position in NDC.
/// @param color  ROOT color index for the marker.
/// @param mstyle ROOT marker style (e.g. kFullCircle).
/// @param text   Null-terminated label string.
/// @param msize  Marker size in ROOT units (default: 2.0).
void myMarkerText(Double_t x, Double_t y, Int_t color, Int_t mstyle,
                  const char* text, Float_t msize = 2.0f);

// ============================================================
//  Graph arithmetic
// ============================================================

/// @brief Compute the ratio of two symmetric-error graphs (g1 / g2).
///
/// Errors are propagated in quadrature assuming the two graphs are
/// uncorrelated.  The caller takes ownership of the returned object.
///
/// @param g1  Numerator graph.   Must not be nullptr.
/// @param g2  Denominator graph. Must not be nullptr.
/// @return    New TGraphErrors representing g1 / g2, or nullptr on failure.
TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1, TGraphErrors* g2);

/// @brief Compute the ratio of two asymmetric-error graphs (g1 / g2).
///
/// High and low errors are propagated independently.  The caller takes
/// ownership of the returned object.
///
/// @param g1  Numerator graph.   Must not be nullptr.
/// @param g2  Denominator graph. Must not be nullptr.
/// @return    New TGraphAsymmErrors representing g1 / g2, or nullptr on failure.
TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,
                                        TGraphAsymmErrors* g2);

// ============================================================
//  Uncertainty band construction
// ============================================================

/// @brief Build an asymmetric uncertainty band from three graphs.
///
/// @p g0 is treated as the central value; @p g1 and @p g2 represent the
/// lower and upper variations respectively.  The resulting band encodes
/// asymmetric errors suitable for shaded-band plots.  The caller takes
/// ownership of the returned object.
///
/// @param g0  Central-value graph.
/// @param g1  Lower-variation graph.
/// @param g2  Upper-variation graph.
/// @return    New TGraphAsymmErrors encoding the band, or nullptr on failure.
TGraphAsymmErrors* myMakeBand(TGraphErrors* g0,
                              TGraphErrors* g1,
                              TGraphErrors* g2);

/// @brief Add a systematic contribution from @p g1 into an existing band @p g2.
///
/// Errors are summed in quadrature so that successive calls accumulate
/// independent uncertainty sources.
///
/// @param g1  Graph whose errors are added to the band.
/// @param g2  Existing band graph that is modified in place.
void myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2);

// ============================================================
//  Histogram conversion
// ============================================================

/// @brief Convert a TH1 histogram to a TGraphErrors.
///
/// Bin centers become x-values, bin contents become y-values, and bin
/// errors become y-errors.  x-errors are set to half the bin width.
/// The caller takes ownership of the returned object.
///
/// @param h1  Source histogram. Must not be nullptr.
/// @return    New TGraphErrors mirroring @p h1, or nullptr on failure.
TGraphErrors* TH1TOTGraph(const TH1* h1);

#endif  // ATLAS_UTILS_H
