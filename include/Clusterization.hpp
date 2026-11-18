#pragma once

#include <vector>
#include <map>
#include <concepts>

#include "DataReader.hpp"   // needed for DataReader::Hit::View
#include "Cluster.hpp"      // Cluster + ChamberClusters

namespace MuonTomography {

// A concept to constrain the type of hits the clusterization can process.
template <typename T>
concept StripLike = requires(const T& h) {
    { h.chamber } -> std::convertible_to<int>;
    { h.view }    -> std::convertible_to<int>;
    { h.strip }   -> std::convertible_to<int>;
    { h.charge }  -> std::convertible_to<double>;
    { h.chargeDist } -> std::convertible_to<std::vector<short>>;
};

class Clusterization {
 public:
    using View = DataReader::Hit::View;

    using ClusterMap    = std::map<int, ChamberClusters>;   // chamber -> {x, y}
    using ClusterMapVec = std::vector<ClusterMap>;          // one per event

    template <typename dataType>
    using HitsPerLayer = std::map<std::pair<int, View>, std::vector<dataType>>;

    struct Config {
        unsigned int maxGap{1};
        double maxChargeCutOff{1};
        double pitch{0.6}; //in mm
    };

    explicit Clusterization(const Config& config) : m_config(config) {}
    ~Clusterization() = default;

    template <typename dataType>
        requires StripLike<dataType>
    ClusterMapVec clusterizeEvents(
        const std::vector<HitsPerLayer<dataType>>& events) const;
    
    template <typename dataType>
        requires StripLike<dataType>
    ClusterMap clusterize(const HitsPerLayer<dataType>& hitsPerLayer) const;

    
    template <typename dataType>
        requires StripLike<dataType>
    static HitsPerLayer<dataType> sortHitsPerChamber(const std::vector<dataType>& hits);

 private:
    Config m_config;

    
};

}  // namespace MuonTomography

#include "Clusterization.ipp"