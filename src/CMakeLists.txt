add_library(mtomo
    DataReader.cpp
    TelescopeGeometryBuilder.cpp
)
add_library(MTomo::mtomo ALIAS mtomo)

target_include_directories(mtomo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(mtomo PUBLIC cxx_std_20)

# ROOT (for DataReader)
target_link_libraries(mtomo PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::TreePlayer ROOT::Hist)

# ACTS (optional)
if(MTOMO_USE_ACTS)
    target_link_libraries(mtomo PUBLIC Acts::Core)
    target_compile_definitions(mtomo PUBLIC MTOMO_HAS_ACTS)
endif()

