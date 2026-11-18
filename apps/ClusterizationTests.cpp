
#include <iostream>
#include <vector>

#include "DataReader.hpp"
#include "Clusterization.hpp"
#include "Cluster.hpp"


#include <TFile.h>
#include <TTree.h>


namespace MuonTomography::ClusterApplication {

    using HitsType = std::vector<std::vector<DataReader::Hit>>;
    using HitsPerLayer = Clusterization::HitsPerLayer<DataReader::Hit>;
    using ClusterMapVec = Clusterization::ClusterMapVec;
    


    int run(int argc, char* argv[]){

        if(argc < 4){
            std::cout << "Usage: " << argv[0]
                      << " <input.root> <treeName> <nEvents> [output.root]\n";
            return 1;
        }

        const std::string inPath = argv[1];
        const std::string inputTreeName = argv[2];
        const unsigned int nEvents = static_cast<unsigned int>(std::stoul(argv[3]));
        const std::string outPath = (argc > 4) ? argv[4] : "clusters.root";

        //  read the raw hits, one std::vector<Hit> per event 
        DataReader::Config datacfg;
        datacfg.inputFile = inPath;
        datacfg.treeName = inputTreeName;
        DataReader reader(datacfg);

        HitsType totalHits = reader.readEvents(nEvents);
        //for evert event sort the hits per chamber's layer
        std::vector<HitsPerLayer> sortedHitsAllEvents;
        sortedHitsAllEvents.reserve(nEvents);

        for(std::size_t i=0; i < totalHits.size(); i++){
            HitsPerLayer layerHits = Clusterization::sortHitsPerChamber(totalHits[i]);
            //we can either clusterize every layer of this particular event on the fly or 
            //sort->create the container with the sorted hits for all the events and call the clusterization
             std::cout<<"after sorting: /n"<<std::endl;
            for(const auto& [chamber, hits] : layerHits){
                std::cout<<"at chamber/view:"<<chamber.first<<"/"<<static_cast<int>(chamber.second)<<std::endl;
                for(const auto& hit : hits){
                 std::cout<<"hit : "<<"chamber:" << hit.chamber<<"view: "<< static_cast<int>(hit.view)<<"strip"<<hit.strip<<std::endl;
                }
            }
            sortedHitsAllEvents.emplace_back(layerHits); 
           
        }

         // configure and build the clusterizer 
        Clusterization::Config cfg;
        cfg.maxGap          = 3;
        cfg.maxChargeCutOff = 20.;
        cfg.pitch           = 0.6;                 
        Clusterization clusterizer(cfg);

        ClusterMapVec clusters = clusterizer.clusterizeEvents(sortedHitsAllEvents);

        TFile fout(outPath.c_str(), "RECREATE");
        TTree tree("clusters", "Micromegas strip clusters");

        int b_event;      // event index
        int b_nClusters;  // number of clusters in this event

    // one element per cluster in the event:
    std::vector<int>    b_chamber, b_view, b_size;
    std::vector<double> b_charge, b_centroid, b_position;

    // strips/charges of ALL clusters in the event, concatenated.
    // Cluster j occupies the slice [offset_j, offset_j + b_size[j]) of these,
    // where offset_j = sum(b_size[0..j-1]).
    std::vector<int>    b_strips;
    std::vector<double> b_charges;

    tree.Branch("event",     &b_event);
    tree.Branch("nClusters", &b_nClusters);
    tree.Branch("chamber",   &b_chamber);
    tree.Branch("view",      &b_view);
    tree.Branch("size",      &b_size);
    tree.Branch("charge",    &b_charge);
    tree.Branch("centroid",  &b_centroid);
    tree.Branch("position",  &b_position);
    tree.Branch("strips",    &b_strips);
    tree.Branch("charges",   &b_charges);

    // append one cluster into the per-event buffers
    auto addCluster = [&](int view, const Cluster& cl) {
        b_chamber.push_back(cl.chamber);
        b_view.push_back(view);
        b_size.push_back(static_cast<int>(cl.hits.size()));
        b_charge.push_back(cl.charge);
        b_centroid.push_back(cl.centroid);
        b_position.push_back(cl.position);
        for (const auto& h : cl.hits) {            // adapt to your Cluster's hit storage
            b_strips.push_back(static_cast<int>(h.strip));
            b_charges.push_back(static_cast<double>(h.charge));
        }
    };

    //  event loop to fill the tree per event ----
    long long nClustersTotal = 0;

    for (std::size_t event = 0; event < clusters.size(); ++event) {

        b_event = static_cast<int>(event);
        b_chamber.clear(); b_view.clear();    b_size.clear();
        b_charge.clear();  b_centroid.clear(); b_position.clear();
        b_strips.clear();  b_charges.clear();

        // gather every cluster of every chamber (X then Y) into the buffers
        for (const auto& [chamber, cc] : clusters[event]) {
            for (const auto& cl : cc.x) addCluster(0, cl);
            for (const auto& cl : cc.y) addCluster(1, cl);
        }

        b_nClusters = static_cast<int>(b_chamber.size());
        nClustersTotal += b_nClusters;

        tree.Fill();   // <-- one entry for the whole event

        if (event % 20000 == 0)
            std::cout << "  " << event << "/" << clusters.size()
                      << "  (" << nClustersTotal << " clusters)\r" << std::flush;
    }

    tree.Write();
    fout.Close();
    std::cout << "\nWrote " << nClustersTotal << " clusters over "
              << clusters.size() << " events to " << outPath << ":clusters\n";
    return 0;
}

}  // namespace MuonTomography::ClusterApplication

int main(int argc, char* argv[]) {
    return MuonTomography::ClusterApplication::run(argc, argv);
}


