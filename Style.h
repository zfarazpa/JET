//
//  @file    Style.h
//  @brief   ATLAS plotting style for ROOT canvases.
//
//  Provides two entry points:
//    - AtlasStyle()    — constructs and returns a configured TStyle object.
//    - SetAtlasStyle() — applies that style to the current ROOT session.
//
//  Covers fonts, margins, tick marks, color palette, stat/title box
//  visibility, and all other graphical settings required for
//  ATLAS publication-quality figures.
//
//  Usage:
//  @code
//    #include "Style.h"
//    SetAtlasStyle();          // call once before creating any canvas
//  @endcode
//
//  Based on the standard ATLAS plotting conventions.
//

#ifndef STYLE_H
#define STYLE_H

#include "TStyle.h"

/// @brief Construct a TStyle object configured to the ATLAS plotting guidelines.
///
/// Configures:
///   - Canvas and pad margins (top, bottom, left, right)
///   - Font family and size for axes, titles, and labels
///   - Tick mark style (inward ticks on all four sides)
///   - Stat and title box visibility (both suppressed)
///   - Marker defaults and color palette
///
/// The caller does @b not take ownership — the returned pointer is managed
/// by ROOT's global style registry (gROOT).
///
/// @return Pointer to the newly created and registered TStyle.
TStyle* AtlasStyle();

/// @brief Apply the ATLAS style to the current ROOT session.
///
/// Calls AtlasStyle() and passes the result to gROOT->SetStyle(), then
/// calls gROOT->ForceStyle() so that any already-loaded objects are
/// immediately updated to reflect the new settings.
///
/// Call this once at the start of a macro or program, before creating
/// any TCanvas or TH1, to ensure all subsequent plots use ATLAS style.
void SetAtlasStyle();

#endif  // STYLE_H
