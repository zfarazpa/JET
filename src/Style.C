//
//  @file    Style.cxx
//  @brief   Implementation of the ATLAS plotting style for ROOT.
//
//  Based on a style file originally developed for the BaBar experiment.
//


#include "chi2Min.C"
#include <iostream>
#include "../include/Style.h"
#include "TROOT.h"
#include "TStyle.h"

// ============================================================

void SetAtlasStyle()
{
    std::cout << "\nApplying ATLAS style settings...\n" << std::endl;

    static TStyle* atlasStyle = nullptr;
    if (!atlasStyle)
        atlasStyle = AtlasStyle();

    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();
}

// ============================================================

TStyle* AtlasStyle()
{
    TStyle* atlasStyle = new TStyle("ATLAS", "Atlas style");

    // ---------------------------------------------------------
    //  Colors — plain black on white throughout
    // ---------------------------------------------------------
    const Int_t kWhiteIdx = 0;
    atlasStyle->SetFrameBorderMode(kWhiteIdx);
    atlasStyle->SetFrameFillColor(kWhiteIdx);
    atlasStyle->SetCanvasBorderMode(kWhiteIdx);
    atlasStyle->SetCanvasColor(kWhiteIdx);
    atlasStyle->SetPadBorderMode(kWhiteIdx);
    atlasStyle->SetPadColor(kWhiteIdx);
    atlasStyle->SetStatColor(kWhiteIdx);
    // Note: SetFillColor is intentionally omitted — applying it globally
    // forces a white fill on every object, which is usually undesirable.

    // ---------------------------------------------------------
    //  Paper & pad sizes
    // ---------------------------------------------------------
    atlasStyle->SetPaperSize(20, 26);       // cm, A4-like portrait

    // ---------------------------------------------------------
    //  Pad margins
    // ---------------------------------------------------------
    atlasStyle->SetPadTopMargin(0.05);
    atlasStyle->SetPadRightMargin(0.05);
    atlasStyle->SetPadBottomMargin(0.16);
    atlasStyle->SetPadLeftMargin(0.16);

    // ---------------------------------------------------------
    //  Axis title offsets
    // ---------------------------------------------------------
    atlasStyle->SetTitleXOffset(1.4);
    atlasStyle->SetTitleYOffset(1.4);

    // ---------------------------------------------------------
    //  Fonts & sizes
    //
    //  Font 42 = Helvetica (plain).
    //  Font 72 = Helvetica italic (uncomment to use instead).
    // ---------------------------------------------------------
    const Int_t    kFont  = 42;
    const Double_t kTSize = 0.05;

    atlasStyle->SetTextFont(kFont);
    atlasStyle->SetTextSize(kTSize);

    for (const char* axis : {"x", "y", "z"}) {
        atlasStyle->SetLabelFont(kFont,  axis);
        atlasStyle->SetTitleFont(kFont,  axis);
        atlasStyle->SetLabelSize(kTSize, axis);
        atlasStyle->SetTitleSize(kTSize, axis);
    }

    // ---------------------------------------------------------
    //  Markers & lines
    // ---------------------------------------------------------
    atlasStyle->SetMarkerStyle(20);               // filled circle
    atlasStyle->SetMarkerSize(1.2);
    atlasStyle->SetHistLineWidth(2);
    atlasStyle->SetLineStyleString(2, "[12 12]"); // dashed (PostScript)

    // Error bar style: no caps, symmetric x-errors kept (uncomment
    // SetErrorX to suppress x-error bars entirely).
    // atlasStyle->SetErrorX(0.001);
    atlasStyle->SetEndErrorSize(0.);

    // ---------------------------------------------------------
    //  Histogram decorations — suppress title, stats, and fit box
    // ---------------------------------------------------------
    atlasStyle->SetOptTitle(0);
    atlasStyle->SetOptStat(0);
    atlasStyle->SetOptFit(0);

    // ---------------------------------------------------------
    //  Tick marks on all four sides
    // ---------------------------------------------------------
    atlasStyle->SetPadTickX(1);
    atlasStyle->SetPadTickY(1);

    // ---------------------------------------------------------
    //  Color palette
    // ---------------------------------------------------------
    atlasStyle->SetPalette(1);
    // Uncomment below if a 2-D color plot needs a wider right margin:
    // atlasStyle->SetPadRightMargin(0.16);

    return atlasStyle;
}
