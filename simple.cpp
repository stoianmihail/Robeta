#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <cmath>
#include <cassert>
#define NORMAL

#ifndef NORMAL  
  #include "../util.h"
#endif
  
using namespace std;

typedef vector<uint32_t> VU;
typedef vector<double> VD;
typedef pair<uint32_t, uint32_t> PUU;
typedef vector<pair<uint32_t, uint32_t>> VUU;
typedef unordered_map<uint32_t, uint32_t> table;

class Robeta {
private:
  uint32_t computeLimit(const vector<double>& values, const double threshold, const int32_t dir = 1) {
    return lower_bound(values.begin(), values.end(), threshold, [&dir](double elem, double x) {
      return (dir == 1) ? (elem < x) : (elem > x);
    }) - values.begin();
  }

  uint32_t computeRho() { 
    uint32_t limitRho = std::min(globalOccupied, static_cast<uint32_t>(log2(memoryLimit * (1u << 20) / 4)));
    
    cerr << "memLimit=" << memoryLimit << " globalOccupied=" << globalOccupied << endl;

    VD ratios;
    double custom;
    VU lastPrefix(globalOccupied + 1), count(globalOccupied + 1, 1);
    for (unsigned index = 0; index != n; ++index) {
      for (unsigned rho = 1; rho <= globalOccupied; ++rho) {
        uint32_t currPrefix = (keys[index] - globalMin) >> (globalOccupied - rho);
        if (currPrefix == lastPrefix[rho])
          continue;
        count[rho]++;
        lastPrefix[rho] = currPrefix;
      }
    }
    for (unsigned rho = 1; rho <= globalOccupied; ++rho) {
      double ratio = 1.0 * count[rho] / (1ull << rho);
      cerr << "rho=" << rho << " count=" << count[rho] << " power=" << (1ull << rho) << " ratio=" << ratio << endl;
      ratios.push_back(ratio);
      custom += ratio;
    }
    custom /= globalOccupied;
    cerr << "custom=" << custom << " place=" << computeLimit(ratios, custom, -1) << endl;

    uint32_t acceptedRho = computeLimit(ratios, custom, -1);
    return (acceptedRho <= limitRho) ? acceptedRho : limitRho;
  }
#if 0
  void buildRadix(uint32_t rho, const VU& keys) {
    uint32_t n = keys.size();
    uint32_t globalMin = keys.front();
    uint32_t globalMax = keys.back();
    uint32_t globalOccupied = 32 - __builtin_clz(globalMax - globalMin);
    uint32_t globalShiftWith = globalOccupied - rho;
    
    uint32_t numBuckets = 0, numSingles = 0, numBelows = 0;
    auto analyzeSlice = [&globalMin, &globalShiftWith, &keys](uint32_t a, uint32_t b) -> uint32_t {
      // [a, b[
      // TODO: what happens for a == b - 1?
      // TODO: repair this assert, since in for 'shiftWith - 1'
      assert(globalShiftWith);
      
      // TODO: do not forget "- min" for each key!
      uint32_t mask = (1u << globalShiftWith) - 1;
      uint32_t localMin = (keys[a] - globalMin) & mask;
      uint32_t localMax = (keys[b - 1] - globalMin) & mask;
      uint32_t localOccupied = 32 - __builtin_clz(localMax - localMin);
      
      // Find the best candidate for beta
      VD ratios;
      double custom = 0;
      for (unsigned cand = 1; cand <= localOccupied; ++cand) {
        uint32_t localShiftWith = localOccupied - cand;
        // TODO: do not forget "- globalMin" for each key!
        uint32_t count = 0, lastPrefix = 0;
        for (unsigned index = a; index != b; ++index) {
        uint32_t currPrefix = (((keys[index] - globalMin) & mask) - localMin) >> localShiftWith;
        count += (currPrefix != lastPrefix);
        lastPrefix = currPrefix;
      }
      count++;
      
      double ratio = 1.0 * count / (1u << cand);
      //cerr << "cand=" << cand << " ratio=" << ratio << endl;
      ratios.push_back(ratio);
      custom += ratio;
    }
    custom /= (b - a);
    //cerr << "custom=" << custom << endl;
    
    uint32_t beta = computeLimit(ratios, custom, -1);
    return beta;
  };
  
  uint32_t lastPrefix = 0, lastPtr = 0;
  double betaMean = 0, betaMax = 0;
  for (unsigned index = 0; index != n; ++index) {
    uint32_t currPrefix = (keys[index] - globalMin) >> globalShiftWith;
    if (currPrefix != lastPrefix) {
      if ((index - lastPtr) == 1) {
        numSingles++;
      } else if (1 ? ((index - lastPtr) > 4) : true) {
        uint32_t ret = analyzeSlice(lastPtr, index);
        //cerr << "slice [" << lastPtr << ":" << index << "] receives (omega=" << ret.first << ", beta=" << ret.second << ")" << endl;  
        betaMean += ret;
        if (ret > betaMax)
          betaMax = ret;
      } else {
        // Put in a SSE register?
        numBelows++;
      }
      lastPrefix = currPrefix;
      lastPtr = index;
      numBuckets++;
    }
  }
  if ((n - lastPtr) == 1) {
    numSingles++;
  } else if (1 ? ((n - lastPtr) > 4) : true) {
    uint32_t ret = analyzeSlice(lastPtr, n);
    // cerr << "slice [" << lastPtr << ":" << n << "] receives (omega=" << ret.first << ", beta=" << ret.second << ")" << endl;  
    betaMean += ret;
    if (ret > betaMax)
      betaMax = ret;
  } else {
    numBelows++;
  }
  numBuckets++;
  
  betaMean /= (numBuckets - numSingles - numBelows);
  cerr << "numBuckets=" << (numBuckets - numSingles - numBelows) << " numSingles=" << numSingles << " numBelows=" << numBelows << " betaMean=" << betaMean << " betaMax=" << betaMax << endl;
}
#endif
  void buildSimpleRadix() {
    radixHint.resize((1u << globalRho) + 1);
    globalShiftWith = globalOccupied - globalRho;
    
    ofstream output("debug.txt");
    output << globalRho << endl;
    uint32_t lastPrefix = 0, lastPtr = 0;
    for (unsigned index = 0; index != n; ++index) {
      uint32_t currPrefix = (keys[index] - globalMin) >> globalShiftWith;
      if (currPrefix == lastPrefix)
        continue;
      for (unsigned prefix = lastPrefix + 1; prefix <= currPrefix; ++prefix)
        radixHint[prefix] = index;
      output << lastPrefix << "," << (index - lastPtr) << endl;
      lastPrefix = currPrefix;
      lastPtr = index;
    }
    for (; lastPrefix != (1u << globalRho); ++lastPrefix)
      radixHint[lastPrefix + 1] = n;
    output << lastPrefix << ", " << (n - lastPtr) << endl;
  }

public:
  string filename;
  VU keys, radixHint;
  const uint32_t memoryLimit;
  uint32_t n, globalMin, globalMax, globalOccupied, globalRho, globalShiftWith;

  Robeta(const uint32_t memoryLimit, const string filename) : memoryLimit(memoryLimit), filename(filename) {
    readKeys();
    n = keys.size();
    globalMin = keys.front();
    globalMax = keys.back();
    globalOccupied = 32 - __builtin_clz(globalMax - globalMin);
    globalRho = computeRho();
    buildSimpleRadix();
    
    cerr << "globalMin=" << globalMin << " globalMax=" << globalMax << " globalRho=" << globalRho << endl;
    // debug();
  }

  void readKeys() {
#ifndef NORMAL
    keys = util::load_data<uint32_t>(filename);
#else
    ifstream input(filename);
    input >> n;
    keys.resize(n);
    for (unsigned index = 0; index != n; ++index)
      input >> keys[index];
#endif
    if (!is_sorted(keys.begin(), keys.end()))
      sort(keys.begin(), keys.end());
	  auto last = std::unique(keys.begin(), keys.end());
	  keys.erase(last, keys.end());
  }
  
  inline uint32_t lookup(const uint32_t x) const {
    if (x <= globalMin)
      return 0;
    uint32_t prefix = (x - globalMin) >> globalShiftWith;
    uint32_t begin = radixHint[prefix], end = radixHint[prefix + 1];

    uint32_t index;
    switch (end - begin) {
      case 0:
        index = begin;
        break;
      case 1:
        index = (x <= keys[begin]) ? begin : end;
        break;
      default:
        index = lower_bound(keys.begin() + begin, keys.begin() + end, x, [](uint32_t key, uint32_t x) {
          return key < x;
        }) - keys.begin();
    }
    return index;
  }

  ~Robeta() {
    keys.clear();
    radixHint.clear();
  }
};

int main(int argc, char** argv) {
	if (argc != 3) {
    cerr << "Check the parameters" << endl;
    exit(1);
  }
  uint32_t memoryLimit = atoi(argv[1]);
  string filename(argv[2]);
  Robeta robeta(memoryLimit, filename);
}
