#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <iterator>
#include <functional>

using namespace std;

typedef pair<double, double> Coord;

template<class K, 
struct treap {
  Coord coord;
  int prior;
  treap *l, *r;
  treap() {}
  treap(Coord coord, int prior) : coord(coord), prior(prior), l(nullptr), r(nullptr) {}
};

struct item {
  double key;
  int prior;
  treap* segment;
  item *l, *r;
  item() {}
  item (int key, int prior) : key(key), prior(prior), l(nullptr), r(nullptr) {}
};
typedef item * pitem;

bool search(pitem t, int key) {
  if (!t)
    return false;
  if (key == t->key) {
    return true;
  } else if (key < t->key) {
    return search(t->l, key);
  } else {
    return search(t->r, key);
  }
}

void split (pitem t, int key, pitem & l, pitem & r) {
    if (!t)
        l = r = NULL;
    else if (key < t->key)
        split (t->l, key, l, t->l),  r = t;
    else
        split (t->r, key, t->r, r),  l = t;
}

void insert (pitem & t, pitem it) {
    if (!t)
        t = it;
    else if (it->prior > t->prior)
        split (t, it->key, it->l, it->r),  t = it;
    else
        insert (it->key < t->key ? t->l : t->r, it);
}

void merge (pitem & t, pitem l, pitem r) {
    if (!l || !r)
        t = l ? l : r;
    else if (l->prior > r->prior)
        merge (l->r, l->r, r),  t = l;
    else
        merge (r->l, l, r->l),  t = r;
}

void erase (pitem & t, int key) {
    if (t->key == key)
        merge (t, t->l, t->r);
    else
        erase (key < t->key ? t->l : t->r, key);
}

pitem unite (pitem l, pitem r) {
    if (!l || !r)  return l ? l : r;
    if (l->prior < r->prior)  swap (l, r);
    pitem lt, rt;
    split (r, l->key, lt, rt);
    l->l = unite (l->l, lt);
    l->r = unite (l->r, rt);
    return l;
}

vector<Coord> generator(size_t size) {
  vector<Coord> data(size);
  
  // First create an instance of an engine.
  random_device rnd_device;
  
  // Specify the engine and distribution.
  mt19937 mersenne_engine {rnd_device()};
  
  // Generates random integers
  uniform_int_distribution<uint64_t> dist {0, numeric_limits<uint64_t>::max()};

  auto gen = [&dist, &mersenne_engine]() {
    return dist(mersenne_engine);
  };

  for (unsigned index = 0; index != size; ++index)
    data[index] = make_pair(gen(), index);
  sort(data.begin(), data.end());
  return data;
}

static vector<Coord> buildCdf(const vector<Coord>& data) {
  vector<Coord> cdf;
  
  // Determine for each distinct key, where its first occurrence appeared
  uint32_t pos = 0, refresh = 0;
  double last = data.front().first;
  for (auto d : data) {
    if (d.first != last) {
      cdf.push_back({last, refresh});
      refresh = pos;
      last = d.first;
    }
    pos++;
  }
  // And add the last point
  cdf.push_back({last, refresh});
  return cdf;
}

static vector<Coord> tautString(const vector<Coord>& data, double epsilon) {
  // Add the first point
  vector<Coord> result;
  if (data.empty())
    return result;
  result.push_back(data.front());
  if (data.size() == 1)
    return result;

  // Taut the string
  vector<Coord>::const_iterator iter = data.begin(), limit = data.end();
  Coord upperLimit,lowerLimit,last=result.back();

  for (++iter; iter!=limit; ++iter) {
    // Add the new bounds
    Coord u = *iter, l = *iter, b = result.back();
    u.second += epsilon; l.second -= epsilon;

    // Check if we cut the error corridor
    if ((last != b) && ((cmpDevs(b, upperLimit, *iter) < 0) || (cmpDevs(b, lowerLimit, *iter) > 0))) {
      result.push_back(last);
      b = last;
    }

    // Update the error margins
    if ((last == b) || (cmpDevs(b, upperLimit, u) > 0))
      upperLimit = u;
    if ((last == b) || (cmpDevs(b, lowerLimit, l) < 0))
      lowerLimit = l;

    // And remember the current point
    last = *iter;
  }

  // Add the last point
  result.push_back(*(--data.end()));
  return result;
}

int main(int argc, char** argv) {
  if (!argc) {
    cerr << "Did you forget the size of data?" << endl;
    return -1;
  }
  pitem t = nullptr;
  
  unsigned n = 10;
  item elem;
  vector<item> elems;
  for (unsigned index = 0; index != n; ++index) {
    elems.emplace_back(item(index, rand()));
    insert(t, &elems.back());
  }
  
  for (unsigned index = 0; index != n; ++index) {
    cerr << index << " : " << search(t, index) << endl;
  }
  
  
  pitem l, r;
  split(t, 5, l, r);
  cerr << "after split" << endl;
  
  size_t size = atoi(argv[1]);
  vector<Coord> data = generator(size);

#if 0  
  for (unsigned index = 0; index != size; ++index) {
    cerr << data[index].first << ", " << data[index].second << endl;
  }
#endif

  vector<Coord> cdf = buildCdf(data);
  vector<Coord> spline = tautString(cdf, 2);
  
  unsigned ptr = 0;
  for (unsigned index = 0, limit = spline.size(); index != limit - 1; ++index) {
    treap *segment = nullptr;
    while ((ptr != cdf.size()) && (cdf[index] != spline[index + 1])) {
      insert(segment, cdf[index]);
      ++ptr;
    }
  }
  
  for (auto elem : cdf) {
    if (elem == currKnot) {
      
    }
  }
}
