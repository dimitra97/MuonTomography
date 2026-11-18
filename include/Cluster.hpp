#pragma once
#include <vector>
#include <string>

#include "DataReader.hpp"

namespace MuonTomography{

struct Cluster{

    using View = DataReader::Hit::View;

    //the strips belonging to the cluster
    std::vector<DataReader::Hit> hits;
    //the chamber the cluster belongs to (0,1,2 or 3)
    int chamber;
    //the X or Y view of the cluster
    View view;
    //the average charge of the cluster
    double averageCharge;
    //the average time of the cluster
    double averageTime;
    //centroid from the charge weighting strip 
    double centroid;
    //cluster position in the chamber 
    double position;
    //the total charge of the cluster
    double charge;

};

struct ChamberClusters {
    std::vector<Cluster> x;   
    std::vector<Cluster> y;   
};

}