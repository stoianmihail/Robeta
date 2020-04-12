#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cassert>
#include <cctype>
#include <cstdio>

using namespace std;

/*
 * Robeta - "https://github.com/stoianmihail/Robeta"
 * Further improvements are coming, just check up the repository.
 * The code has been initially released under "https://github.com/learnedsystems/SOSD"
 * PS: Parsing the input was really necessary.
 */

static constexpr uint32_t MAX_N = 100000;
static constexpr uint32_t globalRho = 15;

uint32_t n, globalMin, globalMax, globalOccupied, globalShiftWith;
uint32_t keys[MAX_N];
uint32_t radixHint[(1u << globalRho) + 1];

void buildSimpleRadix() {
  uint32_t lastPrefix = 0;
  for (unsigned index = 0; index != n; ++index) {
    uint32_t currPrefix = (keys[index] - globalMin) >> globalShiftWith;
    if (currPrefix == lastPrefix)
      continue;
    for (unsigned prefix = lastPrefix + 1; prefix <= currPrefix; ++prefix)
      radixHint[prefix] = index;
    lastPrefix = currPrefix;
  }
  for (; lastPrefix != (1u << globalRho); ++lastPrefix)
    radixHint[lastPrefix + 1] = n;
}

inline int32_t lookup(const uint32_t type, const uint32_t x) {
  /*
   * There are 3 types of lookups:
   * Type 0: if 'x' is inside, return the rightmost index where 'x' appears. Otherwise -1.
   * Type 1: return the biggest position of an element, which is lower or equal than 'x'
   * Type 2: return the smallest position of an element, which is greater or equal than 'x'
  */
  if ((type == 0) && (x < globalMin))
    return -1;
  if ((type == 1) && (x >= globalMax))
    return n;
  if ((type == 2) && (x <= globalMin))
    return 1;
  
  uint32_t prefix = (x - globalMin) >> globalShiftWith;
  uint32_t begin = radixHint[prefix], end = radixHint[prefix + 1];

  int32_t index = -1;
  switch (type) {
    case 0 : {
      index = upper_bound(keys + begin, keys + end, x) - keys - 1;
      if (keys[index] != x)
        return -1;
      break;
    }
    case 1 : {
      index = lower_bound(keys + begin, keys + end, x + 1) - keys - 1;
      break;
    }
    case 2 : {
      index = upper_bound(keys + begin, keys + end, x - 1) - keys;
      break;
    }
  }
  return index;
}

int main() {
  cin >> n;
  for (unsigned index = 0; index != n; ++index)
    cin >> keys[index];
  globalMin = keys[0];
  globalMax = keys[n - 1];
  globalOccupied = 32 - __builtin_clz(globalMax - globalMin);
  globalShiftWith = globalOccupied - globalRho;
  
  buildSimpleRadix();
  
  uint32_t m;
  cin >> m;
  while (m--) {
    unsigned type, x;
    cin >> type >> x;
    cout << lookup(type, x) << '\n';
  }
  return 0;
}
