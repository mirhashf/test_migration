#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <cassert>
#include <boost/algorithm/string.hpp>

struct ChainBlock {
  int x, y, size;
};

bool CompareChainBlock(const int x, const ChainBlock& cb) {
  return cb.x > x;
}

std::vector<std::string> split(std::string str, std::string delim) { 
  size_t start = 0;
  size_t end; 
  std::vector<std::string> v; 

  while( (end = str.find(delim, start)) != std::string::npos ) { 
    v.push_back(str.substr(start, end-start)); 
    start = end + delim.length();
    std::cerr << "start = " << start << " end = " << end << " delim = " << delim << " pushed " << v[v.size() - 1] << " npos " << std::string::npos <<  std::endl;
  } 
  v.push_back(str.substr(start)); 
  return v; 
}

int lift_over_location(std::vector<ChainBlock>& chainBlocks, int location) {
    std::vector<ChainBlock>::iterator it = std::upper_bound(chainBlocks.begin(), chainBlocks.end(), location, CompareChainBlock);
    size_t offset = std::distance(chainBlocks.begin(), it);
    //std::cout << "block found: " << it->x << " " << it->y << " " << it->size << " offset = " << offset << std::endl;
    offset--;
    const ChainBlock& foundBlock = chainBlocks[offset];
    if (foundBlock.x <= location && foundBlock.x + foundBlock.size > location) {
      return foundBlock.y + location - foundBlock.x;
    }
    return -1;
}

int main(int argc, char** argv) {
  std::string chain_file_name = argv[1];
  std::string line;

  std::ifstream chain_ifs(argv[1]);

  std::unordered_map<std::string, std::pair<std::string, std::vector<ChainBlock>* > > chrToChainMap;
  while (true) {
    std::getline(chain_ifs, line);
    boost::algorithm::trim(line);
    if (line.empty()) {
      break;
    }
    //std::cout << "chain: " << line << std::endl;
    std::stringstream ss(line);

    std::string chain, tName, tStrand, qName, qStrand;
    int score, tSize, tStart, tEnd, qSize, qStart, qEnd, id;
    ss >> chain >> score >> tName >> tSize >> tStrand >> tStart >> tEnd >> qName >> qSize >> qStrand >> qStart >> qEnd >> id;

    //std::cout << "Chain " << id << " maps " << tName << ":" << tStart << "-" << tEnd << " to " << qName << ":" << qStart << "-" << qEnd << std::endl;

    int ungapped_sum = 0;
    int tSum = 0;
    int qSum = 0;
    std::vector<ChainBlock>* chainBlocks = new std::vector<ChainBlock>();
    ChainBlock cb;
    while (true) {
      std::getline(chain_ifs, line);
      boost::algorithm::trim(line);
      if (line.empty()) {
        break;
      }
      ss.str(line);
      ss.clear();

      int size = 0, dq = 0, dt = 0;
      ss >> size;
      ss >> dt;
      if ((ss.rdstate() & std::ifstream::failbit) || (ss.rdstate() & std::ifstream::eofbit)) {
        dq = 0; dt = 0;
      } else {
        ss >> dq;
      }
      cb.x = qSum; cb.y = tSum; cb.size = size;
      chainBlocks->push_back(cb);
      //std::cout << size << " " << dt << " " << dq << std::endl;
      ungapped_sum += size;
      tSum += dt + size;
      qSum += dq + size;
    }
    cb.x = qSum; cb.y = tSum; cb.size = 0;
    chainBlocks->push_back(cb);
    
    //std::cout << ungapped_sum << " " << tSum << " " << qSum << std::endl;

    chrToChainMap[qName] = std::pair<std::string, std::vector<ChainBlock>* >(tName, chainBlocks);
  }

  std::cerr << "Will start converting now " << std::endl;
  while(true) {
    std::getline(std::cin, line);
    std::cerr << "line is: " << line << std::endl;
    if (line.empty()) break;

    bool skip = false;
    if (line.find("@rand_") != 0) {
      std::vector<std::string> name_fields = split(line, ":");
      std::vector<std::string> first_parts = split(name_fields[0], "_");

      std::string chr = first_parts[0];
      for (int i = 1; i < first_parts.size() - 7; i++) {
        chr += "_" + first_parts[i];
      }

      int loc1 = atoi(first_parts[first_parts.size() - 7].c_str());
      int loc2 = atoi(first_parts[first_parts.size() - 6].c_str());

      std::cerr << "converting " << chr.substr(1) << " " << loc1 << " " << loc2 << std::endl;

      std::unordered_map<std::string, std::pair<std::string, std::vector<ChainBlock>* > >::iterator it = chrToChainMap.find(chr.substr(1));
      if (it != chrToChainMap.end()) {

        std::string chrNew = it->second.first;
        int loc1_new = lift_over_location(*(it->second.second), loc1);
        int loc2_new = lift_over_location(*(it->second.second), loc2);

        if (loc1_new == -1 || loc2_new == -1) skip = true;

        if (!skip) {
          std::cout << "@" << chrNew << "_" << loc1_new << "_" << loc2_new;
          for (int i = first_parts.size() - 5; i < first_parts.size(); i++) {
            std::cout << "_" << first_parts[i];
          }
          for (int i = 1; i < name_fields.size(); i++) {
            std::cout << ":" << name_fields[i];
          }
          std::cout << std::endl;
        } else {
          std::cerr << "Skipping " << line << " since couldn't lift over" << std::endl;
        }
      } else {
        skip = true;
        std::cerr << "Skipping " << line << " since chromosome couldn't be found" << std::endl;
      }
    } else {
      std::cout << line << std::endl;
    }
    for (int i = 0; i < 3; i++) {
      std::getline(std::cin, line);
      if (!skip) std::cout << line << std::endl;
    }
  }  
}
