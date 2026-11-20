
#pragma once

#include <algorithm>
#include <iostream>

#include "Clusterization.hpp"

namespace MuonTomography {

// -------- group one event's hits by (chamber, view) --------
template <typename dataType>
    requires StripLike<dataType>
Clusterization::HitsPerLayer<dataType>
Clusterization::sortHitsPerChamber(const std::vector<dataType>& hits) {
    HitsPerLayer<dataType> sortedHits;
    for (const auto& h : hits) {
        const int  chamber = static_cast<int>(h.chamber);
        const View layer   = static_cast<View>(h.view);
        sortedHits[{chamber, layer}].push_back(h);   // was: sortedHits[{chamber, view}]
    }
    return sortedHits;
}

// -------- strip selection --------
template <typename dataType>
    requires StripLike<dataType>
bool Clusterization::acceptStrip(const dataType& hit) const {
    // time window: the peak of a muon signal sits at a fixed delay after the
    // trigger, noise peaks at random times
    const int t = static_cast<int>(hit.time);
    if (t < m_config.tMin || t > m_config.tMax) return false;

    // raw_q cut: if a sample is below the cut-off consider it a noise strip
    if (m_config.useRawQCut &&
        std::ranges::any_of(hit.chargeDist, [&](short q) { return q < m_config.maxChargeCutOff; }))
        return false;
    return true;
}

// -------- build clusters from one event's grouped hits --------
template <typename dataType>
    requires StripLike<dataType>
Clusterization::ClusterMap
Clusterization::clusterize(const HitsPerLayer<dataType>& hitsPerLayer) const {
    ClusterMap result;

    for (const auto& [layer, layerHitsConst] : hitsPerLayer) {
        const int  chamber = layer.first;
        const View view    = layer.second;
        if (layerHitsConst.empty()) continue;

        // hitsPerLayer is const -> take a mutable copy before sorting
        std::vector<dataType> layerHits = layerHitsConst;
        std::sort(layerHits.begin(), layerHits.end(),
                  [](const dataType& a, const dataType& b) {
                      return static_cast<int>(a.strip) < static_cast<int>(b.strip);
                  });

        std::vector<Cluster>  clusters;
        std::vector<dataType> clusterHits;   // hits of the current cluster

        // finalize the current run of strips into a Cluster 
        auto dumpCluster = [&]() -> bool {
            if (clusterHits.empty()) return false;

            double totCharge = 0.0;
            double centroid  = 0.0;    // numerator; normalized below
            double maxq      = 0.0;
            for (const auto& cs : clusterHits) {
                const double q = static_cast<double>(cs.charge);
                maxq       = std::max(maxq, q);
                centroid  += static_cast<int>(cs.strip) * q;
                totCharge += q;
            }


            // charge-weighted centroid 
            centroid = (totCharge > 0.0)
                         ? centroid / totCharge
                         : static_cast<double>(static_cast<int>(clusterHits.front().strip));

            Cluster cl;
            cl.chamber  = chamber;
            cl.view     = view;
            cl.hits     = clusterHits;
            cl.charge = totCharge;    
            cl.centroid = centroid;
            cl.position = m_config.pitch * centroid;

            std::cout << "Cluster created for chamber/layer: " << cl.chamber << "/"
                      << static_cast<int>(cl.view) << " | hits=" << cl.hits.size()
                      << " charge=" << cl.charge << " centroid=" << cl.centroid
                      << " position=" << cl.position << std::endl;

            clusters.push_back(std::move(cl));
            clusterHits.clear();
            return true;
        };

        // drop noise strips first, then group neighbours (gap <= maxGap)
        std::erase_if(layerHits, [&](const dataType& h) { return !acceptStrip(h); });
        if (layerHits.empty()) continue;

        int last = static_cast<int>(layerHits[0].strip);
        clusterHits.push_back(layerHits[0]);
        for (std::size_t i = 1; i < layerHits.size(); ++i) {
            const int strip = static_cast<int>(layerHits[i].strip);
            const unsigned int dist = static_cast<unsigned int>(strip - last);
            if (dist > m_config.maxGap) {
                dumpCluster();                        // flush the finished cluster
                clusterHits.push_back(layerHits[i]);  // start a new one
                last = strip;
                continue;
            }
            clusterHits.push_back(layerHits[i]);
            last = strip;
        }
        dumpCluster();  // last cluster

        if (clusters.empty()) continue;

        if (view == View::X) result[chamber].x = std::move(clusters);
        else                 result[chamber].y = std::move(clusters);
    }
    return result;
}

// -------- cluster a set of pre-grouped events --------
template <typename dataType>
    requires StripLike<dataType>
Clusterization::ClusterMapVec
Clusterization::clusterizeEvents(
    const std::vector<HitsPerLayer<dataType>>& events) const {
    ClusterMapVec result;
    result.reserve(events.size());
    for (const auto& hitsInEvent : events) {
        ClusterMap clusters = clusterize(hitsInEvent);   // compute once
        if (!clusters.empty()) {
            result.push_back(std::move(clusters));       // was: clusterize() again
        }
    }
    return result;
}

}  // namespace MuonTomography