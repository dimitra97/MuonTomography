// visualizeClusters.C
//
// Visualize Micromegas strip clusters from a per-event `clusters` tree.
// For a given event it fills two canvases (X view and Y view); each canvas has
// one pad per chamber, and every pad shows that layer's clusters as
// charge-vs-strip bars (one colour per cluster) with a marker at each centroid.
//
// Expected `clusters` tree layout (ONE ENTRY PER EVENT), as written by the
// C++ application:
//     event      int
//     nClusters  int
//     chamber    std::vector<int>     (one per cluster)
//     view       std::vector<int>     (0 = X, 1 = Y)
//     size       std::vector<int>     (strips per cluster)
//     charge     std::vector<double>
//     centroid   std::vector<double>
//     position   std::vector<double>
//     strips     std::vector<int>     (flattened across clusters; slice with size)
//     charges    std::vector<double>  (flattened, aligned with strips)
//
// Usage (interactive ROOT):
//     root -l
//     .L visualizeClusters.C
//     visualizeClusters("clusters.root", 27);      // draw a single event
//     visualizeClustersLoop("clusters.root", 0, 20);// step through events 0..19
//
// or in one line:
//     root -l 'visualizeClusters.C("clusters.root", 27)'

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMarker.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {
constexpr int kNChambers = 4;    // 2x2 pad grid per canvas
constexpr int kNStrips   = 640;  // strip axis range

int clusterColor(int j) {
    static const int pal[] = {kRed + 1,   kBlue + 1, kGreen + 2, kMagenta + 1,
                              kOrange + 7, kCyan + 2, kViolet,    kAzure + 1};
    return pal[j % 8];
}
}  // namespace

class ClusterViewer {
 public:
    explicit ClusterViewer(const char* filename) {
        f_ = TFile::Open(filename, "READ");
        if (!f_ || f_->IsZombie()) {
            std::cerr << "ClusterViewer: cannot open " << filename << "\n";
            return;
        }
        t_ = dynamic_cast<TTree*>(f_->Get("clusters"));
        if (!t_) {
            std::cerr << "ClusterViewer: no 'clusters' tree in " << filename << "\n";
            return;
        }
        t_->SetBranchAddress("event",    &b_event_);
        t_->SetBranchAddress("chamber",  &v_chamber_);
        t_->SetBranchAddress("view",     &v_view_);
        t_->SetBranchAddress("size",     &v_size_);
        t_->SetBranchAddress("charge",   &v_charge_);
        t_->SetBranchAddress("centroid", &v_centroid_);
        t_->SetBranchAddress("strips",   &v_strips_);
        t_->SetBranchAddress("charges",  &v_charges_);

        TH1::AddDirectory(kFALSE);   // histograms we create are ours to delete
        gStyle->SetOptStat(0);

        cX_ = new TCanvas("cX_clusters", "Clusters - X view", 1200, 800);
        cY_ = new TCanvas("cY_clusters", "Clusters - Y view", 1200, 800);
        cX_->Divide(1, 4);
        cY_->Divide(1, 4);
    }

    bool     ok()      const { return f_ && t_; }
    Long64_t entries() const { return t_ ? t_->GetEntries() : 0; }

    bool draw(Long64_t entry) {
        if (!ok() || entry < 0 || entry >= t_->GetEntries()) {
            std::cerr << "draw: entry " << entry << " out of range\n";
            return false;
        }
        // clear canvases BEFORE deleting old objects (removes dangling refs),
        // then delete last event's objects, then rebuild the pad grid.
        cX_->Clear();
        cY_->Clear();
        clearTrash();
        cX_->Divide(1, 4);
        cY_->Divide(1, 4);

        t_->GetEntry(entry);

        // --- unpack the flattened clusters of this event ---
        struct Cl {
            int chamber, view;
            double centroid;
            std::vector<int>    strips;
            std::vector<double> charges;
        };
        std::vector<Cl> cls;
        double ymax[kNChambers][2] = {{0}};

        const int nc = v_size_ ? static_cast<int>(v_size_->size()) : 0;
        if(nc<4){
            return false;

        }
        int offset = 0;
        for (int j = 0; j < nc; ++j) {
            Cl c;
            c.chamber  = (*v_chamber_)[j];
            c.view     = (*v_view_)[j];
            c.centroid = (*v_centroid_)[j];
            const int n = (*v_size_)[j];
            double localMax = 0.0;
            for (int k = 0; k < n; ++k) {
                const int    s = (*v_strips_)[offset + k];
                const double q = (*v_charges_)[offset + k];
                c.strips.push_back(s);
                c.charges.push_back(q);
                localMax = std::max(localMax, q);
            }
            offset += n;
            if (c.chamber >= 0 && c.chamber < kNChambers && (c.view == 0 || c.view == 1))
                ymax[c.chamber][c.view] = std::max(ymax[c.chamber][c.view], localMax);
            cls.push_back(std::move(c));
        }

        // --- draw an axis frame in every pad ---
        for (int c = 0; c < kNChambers; ++c) {
            for (int vw = 0; vw < 2; ++vw) {
                TCanvas* can = (vw == 0) ? cX_ : cY_;
                can->cd(c + 1);
                const double ym = (ymax[c][vw] > 0) ? ymax[c][vw] * 1.25 : 10.0;
                TString title = Form("Event %d  |  Chamber %d  %s;strip;charge",
                                     b_event_, c, (vw == 0 ? "X" : "Y"));
                gPad->DrawFrame(0.5, 0.0, kNStrips + 0.5, ym, title);
            }
        }

        // --- draw each cluster as a coloured histogram + centroid marker ---
        int colorIdx[kNChambers][2] = {{0}};
        for (const auto& c : cls) {
            if (c.chamber < 0 || c.chamber >= kNChambers ||
                (c.view != 0 && c.view != 1))
                continue;
            TCanvas* can = (c.view == 0) ? cX_ : cY_;
            can->cd(c.chamber + 1);

            const int idx = colorIdx[c.chamber][c.view]++;
            const int col = clusterColor(idx);

            TH1F* h = new TH1F(Form("h_e%d_c%d_v%d_%d", b_event_, c.chamber, c.view, idx),
                               "", kNStrips, 0.5, kNStrips + 0.5);
            double peak = 0.0;
            for (size_t k = 0; k < c.strips.size(); ++k) {
                h->SetBinContent(c.strips[k], c.charges[k]);   // bin i == strip i
                peak = std::max(peak, c.charges[k]);
            }
            h->SetFillStyle(1001);
            h->Draw("hist same");
            trash_.push_back(h);

    

            TLine* l = new TLine(c.centroid, 0.0, c.centroid, peak);
            l->SetLineStyle(2);
            l->SetLineColor(col);
            l->Draw();
            trash_.push_back(l);
        }

        cX_->Update();
        cY_->Update();
        return true;
    }

 private:
    void clearTrash() {
        for (TObject* o : trash_) delete o;
        trash_.clear();
    }

    TFile*   f_  = nullptr;
    TTree*   t_  = nullptr;
    TCanvas* cX_ = nullptr;
    TCanvas* cY_ = nullptr;

    int b_event_ = 0;
    std::vector<int>*    v_chamber_  = nullptr;
    std::vector<int>*    v_view_     = nullptr;
    std::vector<int>*    v_size_     = nullptr;
    std::vector<int>*    v_strips_   = nullptr;
    std::vector<double>* v_charge_   = nullptr;
    std::vector<double>* v_centroid_ = nullptr;
    std::vector<double>* v_charges_  = nullptr;

    std::vector<TObject*> trash_;   // per-event objects, deleted between events
};

// ---- single event ----
void visualizeClusters(const char* filename = "clusters.root", Long64_t event = 0) {
    ClusterViewer* viewer = new ClusterViewer(filename);   // kept alive after macro returns
    if (!viewer->ok()) return;
    viewer->draw(event);
}

// ---- step through a range of events (Enter = next, q+Enter = quit) ----
void visualizeClustersLoop(const char* filename = "clusters.root",
                           Long64_t startEvent = 0, Long64_t nEvents = -1) {
    ClusterViewer* viewer = new ClusterViewer(filename);
    if (!viewer->ok()) return;
    const Long64_t total = viewer->entries();
    const Long64_t last  = (nEvents < 0) ? total
                                         : std::min<Long64_t>(total, startEvent + nEvents);
    for (Long64_t e = startEvent; e < last; ++e) {
       bool draw = viewer->draw(e);
       if(!draw){
        continue;
       }
        gSystem->ProcessEvents();
        std::cout << "Event " << e << "/" << (last - 1)
                  << "   [Enter]=next, q+Enter=quit: " << std::flush;
        std::string s;
        std::getline(std::cin, s);
        if (s == "q" || s == "Q") break;
    }
}