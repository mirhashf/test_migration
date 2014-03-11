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

  //std::cout << "Will start converting now " << std::endl;
  while(true) {
    std::getline(std::cin, line);
    //std::cout << "line is: " << line << std::endl;
    if (line.empty()) break;

    std::stringstream ss(line);
    std::string chr;
    int loc;
    ss >> chr >> loc;
    //std::cout << "chr = " << chr << " loc = " << loc << std::endl;
    std::string dstChr = chrToChainMap[chr].first;
    std::vector<ChainBlock>* chainBlocks = chrToChainMap[chr].second;

    std::vector<ChainBlock>::iterator it = std::upper_bound(chainBlocks->begin(), chainBlocks->end(), loc, CompareChainBlock);
    size_t offset = std::distance(chainBlocks->begin(), it);
    //std::cout << "block found: " << it->x << " " << it->y << " " << it->size << " offset = " << offset << std::endl;
    offset--;
    const ChainBlock& foundBlock = (*chainBlocks)[offset];
    if (foundBlock.x <= loc && foundBlock.x + foundBlock.size > loc) {
      std::cout << chr << " " << loc << " maps to " << dstChr << " " << foundBlock.y + loc - foundBlock.x << std::endl;
    } else {
      std::cout << chr << " " << loc << " outside of " << foundBlock.x << ", " << foundBlock.x + foundBlock.size << std::endl;
    }
    std::cout << std::endl;
  }  
}
