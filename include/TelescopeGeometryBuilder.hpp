#pragma once

#include<vector>


#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"






//This class builds the geometry of a telescope detector using the ACTS geometry interface (Blueprint)
//The geometry construction is configurable 
namespace MuonTomography{

class TelescopeGeometryBuilder {

    public:

    struct Config{

        //the surfaces (planes) of the telescope detector
        std::vector<std::shared_ptr<Acts::Surface>> surfaces{};
        //the bounds of the full telescope
        std::shared_ptr<Acts::VolumeBounds> bounds=nullptr;
        //the transform of the detector
        Acts::Transform3 transform = Acts::Transform3::Identity();
        //the name of the final volume
        std::string name="telescope";

    };

    explicit TelescopeGeometryBuilder(const Config& cfg):m_cfg(cfg){}
    ~TelescopeGeometryBuilder() = default;


    //it constructs the final tracking geometry
    std::unique_ptr<Acts::TrackingGeometry> buildTelescopeTrackingGeometry(
      const Acts::GeometryContext& gctx);


    private:

    Config m_cfg;

      //It constructs the node of the tracking volume with the bounds and the surfaces from the config
    std::unique_ptr<Acts::TrackingVolume> buildTelescopeVolume() const;


};

}