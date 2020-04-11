#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <string>
#include <cmath>
#include <cassert>

using namespace std;

typedef vector<uint32_t> VU;
typedef vector<double> VD;
typedef pair<uint32_t, uint32_t> PUU;
typedef vector<pair<uint32_t, uint32_t>> VUU;
typedef unordered_map<uint32_t, uint32_t> table;

uint32_t computeLimit(vector<double>& values, double threshold, int dir = 1) {
	unsigned index = 1;
	for (auto elem : values) {
		if (dir == 1 && elem > threshold)
			break;
    if (dir == -1 && elem < threshold)
      break;
		++index;
	}
	return index - 1;
}

uint32_t findRho(VU& keys) {
  sort(keys.begin(), keys.end());
	auto last = std::unique(keys.begin(), keys.end());
	keys.erase(last, keys.end()); 
    
	uint32_t n = keys.size();
	uint32_t min = keys.front();
	uint32_t max = keys.back();
	uint32_t occupied = 32 - __builtin_clz(max - min);
	
  VD ratios;
  double custom = 0;
	for (unsigned rho = 1; rho <= occupied; ++rho) {
		uint32_t shiftWith = occupied - rho;
		uint32_t lastPrefix = 0, count = 0;
    for (unsigned ptr = 0; ptr != n; ++ptr) {
      uint32_t currPrefix = (keys[ptr] - min) >> shiftWith;
      count += (currPrefix != lastPrefix);
      lastPrefix = currPrefix;
    }
    count++;
    
		double ratio = 1.0 * count / (1u << rho);
    ratios.push_back(ratio);
    custom += ratio;
  }
	custom /= (occupied - 1);
  return computeLimit(ratios, custom, -1);
}

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
      } else if (0 ? ((index - lastPtr) > 4) : true) {
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
  } else if (0 ? ((n - lastPtr) > 4) : true) {
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

int main(int argc, char** argv) {
	unsigned N;
	cin >> N;
	VU keys(N);
	for (unsigned index = 0; index != N; ++index)
		cin >> keys[index];
	
	VU radixHint;
  uint32_t rho = findRho(keys);
  cerr << "N=" << N << " selected rho=" << rho << endl;
	buildRadix(rho, keys); 
}
