#ifndef HASHER_H_
#define HASHER_H_

#include "clhash.h"

#include <array>
#include <cassert>
#include <cstdio>
#include <random>

class CLHasher {
public:
  void init(std::mt19937_64 &gen, int size){
    assert(size >= 1);
    random_data.resize(size);
    for (size_t i = 0; i < random_data.size(); i++) {
      random_data[i] = get_random_key_for_clhash(gen(), gen());
    }

  }
  inline void Hash(void *data, size_t len, uint64_t *out) const {
    for (size_t i = 0; i < random_data.size(); i++) {
      out[i] = clhash(random_data[i], (const char *)data, len);
    }
  }
  inline ~CLHasher() {
    for (size_t i = 0; i < random_data.size(); i++) {
      std::free(random_data[i]);
    }
  }

private:
  std::vector<void *> random_data;
};

#endif /* HASHER_H_ */
