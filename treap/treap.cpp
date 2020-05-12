#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

struct item {
    int key, prior;
    item * l, * r;
    item() { }
    item (int key, int prior) : key(key), prior(prior), l(NULL), r(NULL) { }
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

int main(void) {
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
  
  
  for (unsigned index = 0; index != n; ++index) {
    cerr << index << " : left=" << search(l, index) << " right=" << search(r, index) << endl;
  }
}
