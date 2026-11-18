
#pragma once

#include <string>
#include <vector>

#include <map>
#include <unordered_set>
#include <stdexcept>
#include <memory>
#include <concepts>


namespace MuonTomography {

//This is a class that prepares the data from the raw input data set
//It creates a vector of hits for each event that hold the strips informaiton (chamber layer, charge, time etc)
class DataReader {

    public:

    //The struct for the hit information 
    struct Hit {

        //enum for X or Y layer of the chamber
        enum View { X, Y };
        // the strip number of the latyer
        int strip;
        //the charge of the strip
        short int charge;
         //the raw charge distribution on the strip
        std::vector<short> chargeDist;
        
        //the time information of the strip
        short int time;
        //the chamber the strip belongs to
        int chamber;
        //the X or Y view 
        View view;
       
     
    };

    struct Config {
        //the path of the input file with the raw data
        std::string inputFile{""};
        //the name of the tree in the input file of the raw data
        std::string treeName{""};

    };

    //default constructor
    DataReader() = default;

      ///constructor
    ///@param config the configuration for the data reader
    explicit DataReader(const Config& config);

    //destructor
    ~DataReader();


   ///Reads the raw data from the input file and returns a vector of hits for the first maxEvents events
   ///It calls the readEvent function for each event and stores the hits in a vector of vectors
   ///@param maxEvents the maximum number of events to read
   ///@return a vector of hits for each event
    std::vector<std::vector<Hit>> readEvents(std::size_t maxEvents) const;


   private:

   Config m_config;

   struct ReaderImpl;

   std::unique_ptr<ReaderImpl> m_readerImpl;

   int m_entries{0};

  ///Reads the raw data from the input file and returns a vector of hits for the ith event
  ///@param i the index of the event to read
  ///@return a vector of hits for this event
  std::vector<Hit> readEvent(unsigned int i) const;
    
};

}//namespace