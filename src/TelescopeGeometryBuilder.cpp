#include "TelescopeGeometryBuilder.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Definitions/Units.hpp"

using namespace Acts::UnitLiterals; 

namespace MuonTomography{

    std::unique_ptr<Acts::TrackingGeometry> 
    TelescopeGeometryBuilder::buildTelescopeTrackingGeometry(const Acts::GeometryContext& gctx){

    Acts::Blueprint::Config cfg;
    cfg.envelope[Acts::AxisDirection::AxisZ] = {20_mm, 20_mm};
    cfg.envelope[Acts::AxisDirection::AxisR] = {0_mm, 20_mm};

    auto blueprint = std::make_unique<Acts::Blueprint>(cfg);

    std::unique_ptr<Acts::TrackingVolume> telescopeVolume = buildTelescopeVolume();

    auto& rootNode = blueprint->addStaticVolume(std::move(telescopeVolume));

    Acts::BlueprintOptions blueprintOptions;
    std::unique_ptr<Acts::TrackingGeometry> trackingGeometry = blueprint->construct(blueprintOptions, gctx);

    return trackingGeometry;

    }

    std::unique_ptr<Acts::TrackingVolume>
    TelescopeGeometryBuilder::buildTelescopeVolume() const{

        auto boundsType = m_cfg.bounds->type();
        if(boundsType != Acts::VolumeBounds::eCuboid){
            throw std::invalid_argument("Telescope volume builder: Only cuboidal volume is supported for the telescope detector");
        }

        std::unique_ptr<Acts::TrackingVolume> trackingVolume = std::make_unique<Acts::TrackingVolume>(m_cfg.transform, m_cfg.bounds, m_cfg.name);
        std::vector<std::shared_ptr<Acts::Surface>> surfaces = m_cfg.surfaces;

        for(const auto& surface : surfaces){
            trackingVolume->addSurface(surface);
        }
        return trackingVolume;

    }
}

