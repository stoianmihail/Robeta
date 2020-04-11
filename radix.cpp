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
typedef pair<uint32_t, uint32_t> PUU;
typedef vector<pair<uint32_t, uint32_t>> VUU;
typedef unordered_map<uint32_t, uint32_t> table;

double computeIndex(const string& name, unsigned rho, table& count) {
	VUU pairs(count.size());
	unsigned index = 0;
	for (auto [key, freq] : count)
		pairs[index++] = make_pair(freq, key);
	sort(pairs.begin(), pairs.end());
	
	unsigned total = (1u << rho);
	double sum = 0, mean = 0;
	for (auto [freq, key] : pairs)
		sum += freq;
	mean = sum / total;
		
  // TODO: what is count.size() == 2, and one of them takes it all?
  
	if (name == "gini") {
		unsigned ptr = total - count.size() + 1;
		double div = 1.0 / total, factorSum = 0;
		
		for (auto [freq, key] : pairs)
			factorSum += (ptr++) * freq;
		return div * (2 * factorSum / sum - 1) - 1;
	} else if (name == "theil") {
		double theil = 0;
		for (auto [freq, key] : pairs)
			theil += freq / mean * log(freq / mean);
		theil /= total;
		return theil;
	} else if (name == "sqrt") {
		double atkinson = 0;
		for (auto [freq, key] : pairs)
			atkinson += sqrt(freq);
		atkinson /= total;
		atkinson *= atkinson;
		atkinson /= mean;
		return 1 - atkinson;
	}
}

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
	
	unsigned occupied = 32 - __builtin_clz(max - min);
	
	cerr << "diff=" << (max-min) << " occupied=" << occupied << endl; 
	
	double giniLogMean = 0, theilLogMean, giniMean = 0, theilMean = 0, atkinsonLogMean = 0, harGini = 0, harAtkinson = 0, prop = 0, lower = 0, custom = 0, geoCustom = 0;
	unsigned num = 0;
	vector<double> gs, ts, ats, rs;
	for (unsigned rho = 1; rho != occupied; ++rho) {
		double factor = 1.0 / (1u << rho);
    uint32_t shiftWith = occupied - rho;
		table count;
		for (unsigned ptr = 0; ptr != n; ++ptr) {
			uint32_t currPrefix = (keys[ptr] - min) >> shiftWith;
			++count[currPrefix];
		}
		double ratio = 1.0 * count.size() / (1u << rho);
    custom += ratio;
    geoCustom += log(ratio);
		rs.push_back(ratio);
    
    double gini = computeIndex("gini", rho, count);
		double theil = computeIndex("theil", rho, count);
		double atkinson = computeIndex("sqrt", rho, count);
		cerr << "rho=" << rho << " buckets=" << count.size() << " vs " << (1u << rho) << " ratio = " << (1.0 * count.size() / (1u << rho)) << " --- gini=" << gini << " -- theil=" << theil <<  " --- atkinson(0.5)=" << atkinson << endl;
		giniLogMean += log(gini);
		theilLogMean += log(theil);
		atkinsonLogMean += log(atkinson);
		//giniMean += gini;
		//theilMean += theil;
		gs.push_back(gini);
		ts.push_back(theil);
		ats.push_back(atkinson);
        harGini += 1.0 / gini;
        harAtkinson += 1.0 / atkinson;
    prop += factor * (1 - gini);
		num++;
	}
	
	harGini = num / harGini;
  harAtkinson = num / harAtkinson;
  
  
  custom /= num;
  
  cerr << "try" << endl;
  double currSum = 0;
  unsigned currPtr = 1;
  for (auto elem : rs) {
    currSum += elem;
    cerr << (currPtr++) << ": " << (currSum / custom / num) << endl;
  }
  cerr << endl;
    
  geoCustom = exp(geoCustom / num);
  cerr << "custom =" << custom << " stelle=" << computeLimit(rs, custom, -1) << endl;
  //cerr << "geoCustom=" << geoCustom << " hai=" << computeLimit(rs, geoCustom, -1) << endl;
  
  lower = 1 - 1.0 / (1u << occupied);
  prop /= lower;
  
  cerr << "eigene Mean=" << prop << " and its place=" << computeLimit(rs, prop, -1) << endl;
  
	giniLogMean = exp(giniLogMean / num);
	theilLogMean = exp(theilLogMean / num);
	atkinsonLogMean = exp(atkinsonLogMean / num);
	cerr << "giniLogMean=" << giniLogMean << " -- theilLogMean=" << theilLogMean << " -- atkinson(0.5)LogMean=" << atkinsonLogMean << endl;
	cerr << "gini=" << computeLimit(gs, giniLogMean) << " theil=" << computeLimit(ts, theilLogMean) << " atkinson=" << computeLimit(ats, atkinsonLogMean) << endl;
	cerr << "harGini=" << harGini << " with " << computeLimit(gs, harGini) << "  harAtkinson=" << harAtkinson << " with " << computeLimit(ats, harAtkinson) << endl;
    
  return computeLimit(ats, harAtkinson);
}

void buildRadix(uint32_t rho, const VU& keys) {
  uint32_t n = keys.size();
	uint32_t globalMin = keys.front();
	uint32_t globalMax = keys.back();
  uint32_t globalOccupied = 32 - __builtin_clz(globalMax - globalMin);
	uint32_t globalShiftWith = globalOccupied - rho;
  
  uint32_t numBucket = 0;
  auto analyzeSlice = [&](uint32_t a, uint32_t b) -> PUU {
    // [a, b[
    // TODO: what happens for a == b - 1?
    // TODO: repair this assert, since in for 'shiftWith - 1'
    assert(globalShiftWith);
    
    // TODO: do not forget "- min" for each key!
    uint32_t mask = (1u << globalShiftWith) - 1;
    uint32_t localMin = (keys[a] - globalMin) & mask;
    uint32_t localMax = (keys[b - 1] - globalMin) & mask;
    uint32_t localOccupied = 32 - __builtin_clz(localMax - localMin);
    
    // This is the number of 0s between the prefix and the block
    /*     !---!
     * 1100100000010
     * 1100100000011
     * 1100100001011
     * 1100100001100
     * 1100100011110
     *     !---!
    */
    uint32_t alpha = globalShiftWith - localOccupied;
    
    // Find the best candidate for beta
    vector<double> store;
    double har = 0;
    for (unsigned cand = 1; cand != localOccupied; ++cand) {
      uint32_t localShiftWith = localOccupied - cand;
      table count;
      // TODO: do not forget "- globalMin" for each key!
      for (unsigned index = a; index != b; ++index) {
        uint32_t prefix = ((keys[index] - globalMin - localMin) & mask) >> localShiftWith;
        count[prefix]++;
      }
      double distributionIndex = computeIndex("gini", cand, count);
      
      if (!distributionIndex) {
        cerr << "nein!!!! " << cand << endl; 
      }
      cerr << "cand=" << cand << " localGini=" << distributionIndex << endl;
      store.push_back(distributionIndex);
      har += 1.0 / distributionIndex;
    }
    har = (b - a) / har;
    uint32_t beta = computeLimit(store, har);
    cerr << "localHar=" << har << endl;
    return make_pair(alpha, beta);
  };
  
  uint32_t lastPrefix = 0, lastPtr = 0;
  for (unsigned index = 0; index != n; ++index) {
    uint32_t currPrefix = (keys[index] - globalMin) >> globalShiftWith;
    if (currPrefix != lastPrefix) {
      PUU ret = analyzeSlice(lastPtr, index);
      cerr << "slice [" << lastPtr << ":" << index << "] receives (alfa=" << ret.first << ", beta=" << ret.second << ")" << endl;  
      lastPrefix = currPrefix;
      lastPtr = index;
      numBucket++;
    }
  }
  PUU ret = analyzeSlice(lastPtr, n);
  cerr << "slice [" << lastPtr << ":" << n << "] receives (alfa=" << ret.first << ", beta=" << ret.second << ")" << endl;  
}

int main(int argc, char** argv) {
	unsigned N;
	cin >> N;
	VU keys(N);
	for (unsigned index = 0; index != N; ++index)
		cin >> keys[index];
	
	VU radixHint;
  uint32_t rho = findRho(keys);
  cerr << "Selected rho=" << rho << endl;
	///buildRadix(rho, keys); 
}
