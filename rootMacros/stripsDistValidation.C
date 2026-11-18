
// plotStripsByChamber.C
//
// Plot the fired-strip (mmStrip) distribution per chamber, split into two
// canvases (X view and Y view), each with 4 pads (one per chamber). Only
// strips passing a cut on their raw_q samples are filled.
//
// The cut is applied to the strip's raw_q waveform directly:
//   requireAll = true  -> keep the strip if ALL 12 raw_q samples > rawQCut
//                         (per-strip minimum > cut)
//   requireAll = false -> keep the strip if ANY raw_q sample > rawQCut
//                         (per-strip maximum > cut)
// Either way a passing strip is filled once (we plot the strips that pass).
//
// Usage:
//   root -l 'plotStripsByChamber.C("merged_run1415-1420.root", 60, true)'
//   root -l 'plotStripsByChamber.C("merged_run1415-1420.root", 60, false)'  // ANY

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

void plotStripsByChamber(const char* filename = "../Python/merged_run1415-1420.root",
                         double rawQCut   = 0.0,
                         double maxQCut = 20.0,
                         bool   requireAll = true,      // true=ALL samples>cut, false=ANY
                         const char* treeName = "apv_raw") {

    // plane name -> (chamber 0..3, view: 0 = X, 1 = Y)
    const std::map<std::string, std::pair<int, int>> planeMap = {
        {"ARTEMIS_X",  {0, 0}}, {"APOLLON_Y",  {0, 1}},
        {"BEATRICE_X", {1, 0}}, {"BYRON_Y",    {1, 1}},
        {"CLEO_X",     {2, 0}}, {"CASTOR_Y",   {2, 1}},
        {"DAPHNE_X",   {3, 0}}, {"DIONYSOS_Y", {3, 1}},
    };

    const int kNStrips = 640;

    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << filename << "\n"; return; }
    TTree* t = dynamic_cast<TTree*>(f->Get(treeName));
    if (!t) { std::cerr << "no tree '" << treeName << "'\n"; return; }

    // branches: read the FULL raw_q waveform per strip (int16 -> short)
    std::vector<std::string>*        mmChamber = nullptr;
    std::vector<int>*                mmStrip   = nullptr;
    std::vector<std::vector<short>>* raw_q     = nullptr;
    std::vector<short>* max_q = nullptr;
    t->SetBranchAddress("mmChamber", &mmChamber);
    t->SetBranchAddress("mmStrip",   &mmStrip);
    t->SetBranchAddress("raw_q",     &raw_q);
    t->SetBranchAddress("max_q", &max_q);

    // one histogram per plane
    std::map<std::string, TH1F*> h;
    TH1::AddDirectory(kFALSE);
    const char* mode = requireAll ? "all" : "any";
    for (const auto& [name, cv] : planeMap) {
        h[name] = new TH1F(Form("h_%s", name.c_str()),
                           Form("%s (%s raw_q>%.0f, max_q>%.0f);strip;fired strips",
                                name.c_str(), mode, rawQCut, maxQCut),
                           kNStrips, 0.5, kNStrips + 0.5);
    }

    // fill
    const Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        const size_t m = mmStrip->size();
        for (size_t s = 0; s < m; ++s) {
            auto it = planeMap.find((*mmChamber)[s]);
            if (it == planeMap.end()) continue;

            const std::vector<short>& wf = (*raw_q)[s];   // this strip's raw_q samples
            short maxQ = (*max_q)[s];
            //if (wf.empty()) continue;

            bool pass = false;
            if(rawQCut>0){
                if (requireAll) {
                    // ALL samples > cut  <=>  minimum sample > cut
                    pass = (*std::min_element(wf.begin(), wf.end()) > rawQCut);
                } else {
                    // ANY sample > cut   <=>  maximum sample > cut
                    pass = (*std::max_element(wf.begin(), wf.end()) > rawQCut);
                }
            }else if(maxQCut>0){
                    pass = maxQ>maxQCut;
            }else{
                pass=true;
            }

            if (pass) h[(*mmChamber)[s]]->Fill((*mmStrip)[s]);
        }
        if (i % 20000 == 0)
            std::cout << "  " << i << "/" << nEntries << "\r" << std::flush;
    }
    std::cout << std::endl;

    // draw: X view on one canvas, Y view on the other, chamber -> pad
    gStyle->SetOptStat(0);
    TCanvas* cX = new TCanvas("cX_strips", "mmStrip per chamber - X view", 1200, 800);
    TCanvas* cY = new TCanvas("cY_strips", "mmStrip per chamber - Y view", 1200, 800);
    cX->Divide(2, 2);
    cY->Divide(2, 2);

    for (const auto& [name, cv] : planeMap) {
        const int chamber = cv.first;
        const int view    = cv.second;
        TCanvas* can = (view == 0) ? cX : cY;
        can->cd(chamber + 1);
        h[name]->SetFillColor((view == 0) ? kAzure + 1 : kOrange + 7);
        h[name]->SetLineColor((view == 0) ? kAzure + 2 : kOrange + 8);
        h[name]->Draw("hist");
    }
    cX->cd(); cX->Update();
    cY->cd(); cY->Update();
}