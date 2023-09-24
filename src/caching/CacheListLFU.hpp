/*
 * d4
 * Copyright (C) 2020  Univ. Artois & CNRS
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <boost/program_options.hpp>
#include <vector>

#include "BucketManager.hpp"
#include "CacheCleaningManager.hpp"
#include "CachedBucket.hpp"
#include "src/caching/CacheCleaningManager.hpp"
#include "src/caching/cnf/BucketManagerCnf.hpp"
#include "src/hashing/HashString.hpp"
#include "src/specs/SpecManager.hpp"
#include <3rdParty/plf_colony.h>
#include <set>

namespace d4 {
namespace po = boost::program_options;
template <class T> class CacheListLFU : public CacheManager<T> {
private:
  const unsigned SIZE_HASH = 999331;
  struct Node {
    CachedBucket<T> bucket;
    size_t last_use;
    unsigned hits;
    bool del;
  };

  std::vector<std::vector<Node>> hashTable;
  size_t counter = 0;
  size_t max_mem = (1ull << 30) * 6; // 10GB

public:
  /**
   * @brief Construct a new Cache List object
   *
   * @param vm is a map to get the option.
   * @param nbVar is the number of variables.
   * @param specs is a structure to get data about the formula.
   * @param out is the stream where are printed out the logs.
   */
  CacheListLFU(po::variables_map &vm, unsigned nbVar, SpecManager *specs,
               std::ostream &out)
      : CacheManager<T>(vm, nbVar, specs, out) {
    out << "c [CACHE LIST CONSTRUCTOR]\n";
    initHashTable(nbVar);
  } // constructor

  /**
   * @brief Destroy the Cache List object
   */
  ~CacheListLFU() {} // destructor

  /**
   * @brief Add an entry in the cache.
   *
   * @param cb is the bucket we want to add.
   * @param hashValue is the hash value of the bucket.
   * @param val is the value we associate with the bucket.
   */
  inline void pushInHashTable(CachedBucket<T> &cb, unsigned int hashValue,
                              T val) {
    hashTable[hashValue % SIZE_HASH].push_back(Node{cb, counter, 0, false});
    CachedBucket<T> &cbIn = (hashTable[hashValue % SIZE_HASH].back().bucket);
    cbIn.lockedBucket(val);
    this->m_nbCreationBucket++;
    this->m_sumDataSize += cb.szData();
    this->m_nbEntry++;
    this->counter++;
    if (this->m_bucketManager->usedMemory() > max_mem) {

      std::vector<Node *> last_use;
      for (auto &v : hashTable) {
        for (auto &b : v) {
          last_use.push_back(&b);
        }
      }
      if (last_use.size() < 1) {
        return;
      }
      std::sort(last_use.begin(), last_use.end(),
                [](const Node &a, const Node &b) {
                  if (a.hits == b.hits) {
                    return a.last_use < b.last_use;
                  }
                  return a.hits < b.hits;
                });
      for (size_t i = 0; i < last_use.size() / 2; i++) {
        last_use[i].del = true;
        this->m_bucketManager->releaseMemory(last_use[i].bucket.data,
                                             last_use[i].bucket.szData());
      }
      for (auto &v : hashTable) {
        std::erase_if(v, [&](const Node &n) { return n.del; });
      }
      for (auto &v : hashTable) {
        for (auto &i : v) {
          i.hits /= 2;
        }
      }
      std::cout << "Removed: Mem" << this->m_bucketManager->usedMemory()
                << std::endl;
    }

  } // pushinhashtable

  /**
   * @brief Research in the set of buckets if the bucket pointed by i already
   * exist.
   *
   * @param cb is the bucket we are looking for.
   * @param hashValue is the hash value computed from the bucket.
   * @return a valid entry if it is in the cache, null otherwise.
   */
  CachedBucket<T> *bucketAlreadyExist(CachedBucket<T> &cb, unsigned hashValue) {
    char *refData = cb.data;
    std::vector<Node> &listCollision = hashTable[hashValue % SIZE_HASH];
    for (auto &cbi : listCollision) {
      if (!cb.sameHeader(cbi.bucket))
        continue;

      if (!memcmp(refData, cbi.bucket.data, cbi.bucket.szData())) {
        this->m_nbPositiveHit++;
        cbi.last_use = this->counter;
        return &cbi.bucket;
      }
    }

    this->m_nbNegativeHit++;
    return NULL;
  } // bucketAlreadyExist

  /**
   * Create a bucket and store it in the cache.
   *
   * @param varConnected is the set of variable.
   * @param c is the value we want to store.
   */
  inline void createAndStoreBucket(std::vector<Var> &varConnected, T &c) {
    CachedBucket<T> *formulaBucket =
        this->m_bucketManager->collectBuckect(varConnected);
    unsigned int hashValue = computeHash(*formulaBucket);
    pushInHashTable(*formulaBucket, hashValue, c);
  } // createBucket

  /**
   * @brief Init the hashTable.
   *
   * @param maxVar is the number of variable.
   */
  void initHashTable(unsigned maxVar) override {
    this->setInfoFormula(maxVar);
    // init hash tables
    hashTable.clear();
    hashTable.resize(SIZE_HASH, {});
  } // initHashTable

  /**
   * @brief Clean up the cache.
   *
   * @param test is a function that is used to decide if an entry must be
   * removed.
   *
   * @return the number of entry we removed.
   */
  unsigned removeEntry(std::function<bool(CachedBucket<T> &c)> test) {
    unsigned nbRemoveEntry = 0;
    for (auto &list : hashTable) {
      unsigned j = 0;
      for (unsigned i = 0; i < list.size(); i++) {
        CachedBucket<T> &cb = list[i].bucket;

        if (test(cb)) {
          assert((int)cb.szData() > 0);
          this->releaseMemory(cb.data, cb.szData());
          cb.reset();
          nbRemoveEntry++;
        } else
          list[j++] = list[i];
      }
      list.resize(j);
    }
    return nbRemoveEntry;
  } // removeEntry
};

} // namespace d4
