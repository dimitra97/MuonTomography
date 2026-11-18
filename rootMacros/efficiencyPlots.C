// hitEfficiency.C
//
// Hit efficiency per layer, without tracking:
//
//     efficiency(layer) = (number of trigger events with >= 1 hit in the layer)
//                         / (total number of trigger events)
//
// A "hit" is a strip in the layer with max_q > qMin and tMin <= t_max_q <= tMax.
// Defaults (qMin = 0, no time cut) count any strip present in the raw data.
// Every entry of the apv_raw tree is one trigger event.
//
// Also prints, per chamber, the fraction of triggers with a hit in X AND Y and
// in X OR Y.
//
// Usage:
//   root -l -b -q 'hitEfficiency.C("merged_run1415-1420.root")'
//   root -l -b -q 'hitEfficiency.C("merged_run1415-1420.root", 50, 4, 10)'   // qMin=50, 4<=t<=10
//
// Output: table on stdout, hitEfficiency.root / hitEfficiency.pdf (bar chart).

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TEfficiency.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

void hitEfficiency(const char* rawFile = "merged_run1415-1420.root",
                   double qMin = 0, int tMin = 0, int tMax = 99,
                   const char* outFile = "hitEfficiency.root") {
    const int kN = 8;
    const char* layer[kN] = {"ARTEMIS_X", "BEATRICE_X", "CLEO_X",   "DAPHNE_X",
                             "APOLLON_Y", "BYRON_Y",    "CASTOR_Y", "DIONYSOS_Y"};
    // physical chambers {X layer, Y layer}, top to bottom (from the track fits)
    const int chamber[4][2] = {{0, 4}, {2, 7}, {1, 6}, {3, 5}};

    TFile* f = TFile::Open(rawFile, "READ");
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << rawFile << "\n"; return; }
    TTree* t = dynamic_cast<TTree*>(f->Get("apv_raw"));
    if (!t) { std::cerr << "no apv_raw tree\n"; return; }

    std::vector<std::string>* name = nullptr;
    std::vector<short>*       q    = nullptr;
    std::vector<int>*         tq   = nullptr;
    std::vector<std::vector<short>>* raw  = nullptr;
    t->SetBranchStatus("*", 1);
    t->SetBranchAddress("mmChamber", &name);
    t->SetBranchAddress("max_q", &q);
    t->SetBranchAddress("t_max_q", &tq);
    t->SetBranchAddress("raw_q", &raw);

    std::unordered_map<std::string, int> index;
    for (int l = 0; l < kN; ++l) index[layer[l]] = l;

    const Long64_t nTrig = t->GetEntries();
    long long nHit[kN] = {0};
    long long nAnd[4] = {0}, nOr[4] = {0};

    for (Long64_t i = 0; i < nTrig; ++i) {
        t->GetEntry(i);
        bool hit[kN] = {false};
        
        for (size_t k = 0; k < name->size(); ++k) {
            auto it = index.find((*name)[k]);
            if (it == index.end()) continue;
            if ((*q)[k] <= qMin || (*tq)[k] < tMin || (*tq)[k] > tMax) continue;
            
            const auto& samples = (*raw)[k];
            if (std::any_of(samples.begin(), samples.end(), [](short s) { return s < 20; }))
                continue;
            hit[it->second] = true;
        }
        for (int l = 0; l < kN; ++l) nHit[l] += hit[l];
        for (int c = 0; c < 4; ++c) {
            const bool x = hit[chamber[c][0]], y = hit[chamber[c][1]];
            nAnd[c] += x && y;
            nOr[c]  += x || y;
        }
    }

    auto print = [nTrig](const char* label, long long k) {
        const double e  = double(k) / nTrig;
        const double up = TEfficiency::ClopperPearson(nTrig, k, 0.683, true) - e;
        const double dn = e - TEfficiency::ClopperPearson(nTrig, k, 0.683, false);
        printf("  %-30s %9lld / %-9lld = %6.2f %%  (+%.2f -%.2f)\n", label, k, (long long)nTrig, 100 * e, 100 * up, 100 * dn);
    };

    printf("\nHit efficiency = triggers with >= 1 hit in the layer / all triggers\n");
    printf("hit: max_q > %.0f, %d <= t_max_q <= %d\n\n", qMin, tMin, tMax);
    for (int l = 0; l < kN; ++l) print(layer[l], nHit[l]);
    printf("\nPer chamber:\n");
    for (int c = 0; c < 4; ++c) {
        print(Form("ch%d %s AND %s", c, layer[chamber[c][0]], layer[chamber[c][1]]), nAnd[c]);
        print(Form("ch%d %s OR  %s", c, layer[chamber[c][0]], layer[chamber[c][1]]), nOr[c]);
    }
    printf("\n");

    // bar chart
    gStyle->SetOptStat(0);
    TFile fout(outFile, "RECREATE");
    TH1D h("hHitEff", "Hit efficiency per layer;;triggers with a hit / all triggers", kN, 0, kN);
    for (int l = 0; l < kN; ++l) {
        const double e = double(nHit[l]) / nTrig;
        h.SetBinContent(l + 1, e);
        h.SetBinError(l + 1, std::sqrt(e * (1 - e) / nTrig));
        h.GetXaxis()->SetBinLabel(l + 1, layer[l]);
    }
    h.SetMinimum(0);
    h.SetMaximum(1.05);
    h.SetFillColor(kAzure - 4);
    h.SetBarWidth(0.7);
    h.SetBarOffset(0.15);
    TCanvas c("c", "hit efficiency_maxq>20", 900, 550);
    h.Draw("bar");
    h.Draw("e same");
    TString pdf = outFile;
    pdf.ReplaceAll(".root", ".pdf");
    c.Print(pdf);
    h.Write();
    fout.Close();
    std::cout << "Wrote " << outFile << " and " << pdf << "\n";
}