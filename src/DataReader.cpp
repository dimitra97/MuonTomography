#include "DataReader.hpp"

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>




namespace MuonTomography {

namespace{
//a map that maps the chambers layers (X and Y) names with the chamber number
//WARNING: This is data specific and should be changed if the chamber names are different!
const std::map<std::pair<std::string,std::string>, int> chamberMap = {
    {{"ARTEMIS_X","APOLLON_Y"}, 0},
    {{"CLEO_X","DIONYSOS_Y"}, 1},
    {{"BEATRICE_X","CASTOR_Y"}, 2},
    {{"DAPHNE_X","BYRON_Y"}, 3}
};

}

struct DataReader::ReaderImpl {
    std::unique_ptr<TFile> file;
  TTreeReader reader;

  // Only the branches clustering needs. The view (X/Y) is taken from the plane
  // name, so mmReadout / mmLayer are not read here.
  std::unique_ptr<TTreeReaderValue<std::vector<std::string>>> mmChamber;
  std::unique_ptr<TTreeReaderValue<std::vector<std::int32_t>>> mmStrip;
  std::unique_ptr<TTreeReaderValue<std::vector<std::int16_t>>> maxQ;
  std::unique_ptr<TTreeReaderValue<std::vector<std::int32_t>>> maxT;
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<short>>>> rawQ;
};

DataReader::DataReader(const Config& config): m_config(config) {
    //open the input file and read the number of entries
    //this is data specific and should be changed if the input file format is different
    TFile* inputFile = TFile::Open(m_config.inputFile.c_str());
    if(!inputFile || inputFile->IsZombie()) {
        delete inputFile;
        throw std::runtime_error("Could not open input file: " + m_config.inputFile);
    }
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get(m_config.treeName.c_str()));
    if(!tree) {
        inputFile->Close();
        delete inputFile;
        throw std::runtime_error("Could not find tree: " + m_config.treeName + " in input file: " + m_config.inputFile);
    }
    m_entries = tree->GetEntries();
    if(m_entries <= 0) {
        throw std::runtime_error("No entries found in tree: " + m_config.treeName + " in input file: " + m_config.inputFile);
    }
    m_readerImpl = std::make_unique<ReaderImpl>();
    m_readerImpl->file.reset(inputFile);
    m_readerImpl->reader.SetTree(tree);
    m_readerImpl->mmChamber = std::make_unique<TTreeReaderValue<std::vector<std::string>>>(m_readerImpl->reader, "mmChamber");
    m_readerImpl->mmStrip = std::make_unique<TTreeReaderValue<std::vector<std::int32_t>>>(m_readerImpl->reader, "mmStrip");
    m_readerImpl->maxQ = std::make_unique<TTreeReaderValue<std::vector<std::int16_t>>>(m_readerImpl->reader, "max_q");
    m_readerImpl->maxT = std::make_unique<TTreeReaderValue<std::vector<std::int32_t>>>(m_readerImpl->reader, "t_max_q");
    m_readerImpl->rawQ = std::make_unique<TTreeReaderValue<std::vector<std::vector<short>>>>(m_readerImpl->reader, "raw_q");
}

DataReader::~DataReader() {
    if(m_readerImpl) {
        m_readerImpl->file->Close();
    }

}


std::vector<std::vector<DataReader::Hit>> DataReader::readEvents(std::size_t maxEvents) const {
    if(maxEvents > static_cast<std::size_t>(m_entries)) {
        throw std::out_of_range("maxEvents is greater than the number of entries in the input file");
    }

    std::vector<std::vector<Hit>> hitsPerEvent;
    hitsPerEvent.reserve(maxEvents);

    for(std::size_t i = 0; i < maxEvents; ++i) {
        hitsPerEvent.push_back(readEvent(i));
        std::cout<<"finished event i: "<<i<<std::endl;
    }


    return hitsPerEvent;
}

std::vector<DataReader::Hit> DataReader::readEvent(unsigned int i) const {
    if(i >= m_entries) {
        throw std::out_of_range("Event index out of range: " + std::to_string(i));
    }

    if (m_readerImpl->reader.SetEntry(i) != TTreeReader::kEntryValid) {
        throw std::runtime_error("Could not read event: " + std::to_string(i));
    }
    const auto& chambers = **m_readerImpl->mmChamber;
    const auto& strips = **m_readerImpl->mmStrip;
    const auto& charges = **m_readerImpl->maxQ;
    const auto& times = **m_readerImpl->maxT;
    const auto& rawQ = **m_readerImpl->rawQ;

    if(chambers.size() != strips.size() || chambers.size() != charges.size() || chambers.size() != times.size()) {
        throw std::runtime_error("Inconsistent data sizes in event: " + std::to_string(i));
    }

    std::vector<Hit> hits;
    hits.reserve(chambers.size());
    for(size_t j = 0; j < chambers.size(); ++j) {

        const auto& chamberName = chambers[j];
        const auto& strip = strips[j];
        const auto& charge = charges[j];
        const auto& time = times[j];
        const auto& chargeDist = rawQ[j];

        // Determine the chamber number and view from the chamber name
        int chamberNumber = -1;
        Hit::View view;
        for(const auto& entry : chamberMap) {
            if(entry.first.first == chamberName || entry.first.second == chamberName) {
                chamberNumber = entry.second;
                view = (entry.first.first == chamberName) ? Hit::X : Hit::Y;
                break;
            }
        }
        if(chamberNumber == -1) {
            throw std::runtime_error("Unknown chamber name: " + chamberName);
        }

        std::cout<<"pushing back hit: strip"<<strip<<",charge:"<<charge<<"chamberNumber: "<<chamberNumber<<"view: "<<view<<std::endl;

        hits.push_back({strip, charge, chargeDist, static_cast<short int>(time), chamberNumber, view});
    }

    return hits;

}

}//namespace 