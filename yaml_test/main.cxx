#include <iostream>
#include <fstream>
#include "yaml-cpp/yaml.h"

int main() {

   //std::ifstream fin("c12_input.conf");
   //YAML::Parser parser(fin);
   YAML::Node aNode = YAML::LoadFile("c12_input.conf");
   //std::vector<YAML::Node> nodes = YAML::LoadAllFromFile("c12_input.conf");
   //YAML::Dump(aNode);


   //double HRS_length = doc["HRS_length"].as<double>();
   //std::cout << HRS_length << std::endl;
   //for(int i = 0; i<nodes.size() ; i++) {
   //   std::cout << i << std::endl;
   //   YAML::Node aNode = nodes[i];
   for(YAML::const_iterator it=aNode.begin();it!=aNode.end();++it) {
      if( it->second.IsScalar() ) {
         std::cout << " " << it->first.as<std::string>() << " =  " << it->second.as<std::string>() << "\n";
      }
      if( it->second.IsMap() ) {
         for(YAML::const_iterator it2=it->second.begin();it2!=it->second.end();++it2) {
            std::cout << " " << it2->first.as<std::string>() << " =  " << it2->second.as<std::string>() << "\n";
         }
      }
   }
   //}

   std::cout << aNode["debug"]["MMMS"] <<std::endl;


   //for(YAML::Iterator it=doc.begin();it!=doc.end();++it) {
   //   std::string key, value;
   //   it.first() >> key;
   //   it.second() >> value;
   //   std::cout << "Key: " << key << ", value: " << value << std::endl;

   //}


   std::cout << "Hello World!" << std::endl;
   YAML::Node lineup = YAML::Load("{1B: Prince Fielder, 2B: Rickie Weeks, LF: Ryan Braun}");
   for(YAML::const_iterator it=lineup.begin();it!=lineup.end();++it) {
      std::cout << "Playing at " << it->first.as<std::string>() << " is " << it->second.as<std::string>() << "\n";
   }
   return 0;
}

