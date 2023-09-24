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

#include "src/hashing/hasher.hpp"
#include <boost/program_options.hpp>
#include <vector>

#include "BucketManager.hpp"
#include "CacheCleaningManager.hpp"
#include "CachedBucket.hpp"
#include "sparsepp/spp.h"
#include "src/caching/CacheCleaningManager.hpp"
#include "src/caching/cnf/BucketManagerCnf.hpp"
#include "src/hashing/HashString.hpp"
#include "src/specs/SpecManager.hpp"
#include <set>

namespace d4 {
namespace po = boost::program_options;
template <class T> class CacheListLRU : public CacheManager<T> {
private:
  const unsigned SIZE_HASH = 1024 * 1024;
  struct Node {
    CachedBucket<T> bucket;
    size_t last_use;
  };

  std::vector<std::vector<Node>> hashTable;
  size_t counter = 0;
  size_t max_mem = (1ull << 30) * 5; // 5GB

public:
  /**
   * @brief Construct a new Cache List object
   *
   * @param vm is a map to get the option.
   * @param nbVar is the number of variables.
   * @param specs is a structure to get data about the formula.
   * @param out is the stream where are printed out the logs.
   */
  CacheListLRU(po::variables_map &vm, unsigned nbVar, SpecManager *specs,
               std::ostream &out)
      : CacheManager<T>(vm, nbVar, specs, out) {
    out << "c [CACHE LIST CONSTRUCTOR]\n";
    max_mem = (1ull << 30ull) * vm["cache-fixed-size"].as<int>();
    initHashTable(nbVar);
  } // constructor

  /**
   * @brief Destroy the Cache List object
   */
  ~CacheListLRU() {} // destructor

  /**
   * @brief Add an entry in the cache.
   *
   * @param cb is the bucket we want to add.
   * @param hashValue is the hash value of the bucket.
   * @param val is the value we associate with the bucket.
   */
  inline void pushInHashTable(CachedBucket<T> &cb, size_t hashValue, T val) {
    hashTable[hashValue % SIZE_HASH].push_back(Node{cb, counter});
    CachedBucket<T> &cbIn = (hashTable[hashValue % SIZE_HASH].back().bucket);
    cbIn.lockedBucket(val);
    this->m_nbCreationBucket++;
    this->m_sumDataSize += cb.szData();
    this->m_nbEntry++;
    this->counter++;
    if (this->m_bucketManager->usedMemory() > max_mem) {

      std::vector<size_t> last_use;
      for (auto &v : hashTable) {
        for (auto &b : v) {
          last_use.push_back(b.last_use);
        }
      }
      if (last_use.size() < 1) {
        return;
      }
      std::sort(last_use.begin(), last_use.end());
      size_t median = last_use[last_use.size() / 3];
      size_t removed = 0;

      for (auto &v : hashTable) {
        for (auto &e : v) {
          if (e.last_use < median) {
            this->m_bucketManager->releaseMemory(e.bucket.data,
                                                 e.bucket.szData());
            removed++;
          }
        }
      }

      for (auto &v : hashTable) {
        std::erase_if(v, [&](const Node &n) {
          if (n.last_use < median) {
            return true;
          }
          return false;
        });
      }
      std::cout << "Removed: " << removed << " Mem"
                << this->m_bucketManager->usedMemory() << std::endl;
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
  CachedBucket<T> *bucketAlreadyExist(CachedBucket<T> &cb, size_t hashValue) {
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

template <class T> class CacheListProbLRU : public CacheManager<T> {
private:
  struct Node {
    uint64_t *key;
    size_t last_use = 0;
    T data;
  };
  struct Hasher {
    size_t operator()(const Node &n)const { return n.key[0]; }
  };
  struct KeyEq {
    size_t key_len;
    bool operator()(const Node &lhs, const Node &rhs)const {
      return !memcmp(lhs.key, rhs.key, key_len * sizeof(uint64_t));
    }
  };
  spp::sparse_hash_set<Node, Hasher, KeyEq> table;
  size_t max_mem = (1ull << 30) * 8; // 8gb
  size_t key_len = 2;
  size_t counter = 0;
  static std::mt19937_64 gen;
  CLHasher hasher;
  CachedBucket<T> tmp;

public:
  /**
   * @brief Construct a new Cache List object
   *
   * @param vm is a map to get the option.
   * @param nbVar is the number of variables.
   * @param specs is a structure to get data about the formula.
   * @param out is the stream where are printed out the logs.
   */
  CacheListProbLRU(po::variables_map &vm, unsigned nbVar, SpecManager *specs,
                   std::ostream &out)
      : CacheManager<T>(vm, nbVar, specs, out), table(0, Hasher{}, KeyEq{2}) {
    out << "c [CACHE LIST CONSTRUCTOR]\n";
    max_mem = (1ull << 30ull) * vm["cache-fixed-size"].as<int>();
    std::cout << "Fixed Cache: " << max_mem << std::endl;
    initHashTable(nbVar);
    hasher.init(gen, key_len);

  } // constructor
  CacheListProbLRU(const CacheListProbLRU &) = delete;

  CacheListProbLRU &operator=(const CacheListProbLRU &) = delete;

  /**
   * @brief Destroy the Cache List object
   */
  ~CacheListProbLRU() {} // destructor

  /**
   * @brief Add an entry in the cache.
   *
   * @param cb is the bucket we want to add.
   * @param hashValue is the hash value of the bucket.
   * @param val is the value we associate with the bucket.
   */
  size_t computeHash(CachedBucket<T> &bucket) final {
    uint64_t *key = new uint64_t[key_len];
    hasher.Hash(bucket.data, bucket.szData(), key);
    this->m_bucketManager->releaseMemory(bucket.data, bucket.szData());

    bucket.data = nullptr;

    return reinterpret_cast<size_t>(key);
  }
  size_t node_size() { return (sizeof(Node) + sizeof(uint64_t) * key_len); }
  size_t estimate_mem() { return table.size() * (node_size() + 8); }
  inline void pushInHashTable(CachedBucket<T> &_cb, size_t hashValue, T val) {
    uint64_t *key = reinterpret_cast<uint64_t *>(hashValue);
    table.insert(Node{key, counter, val});
    this->m_nbCreationBucket++;
    this->m_sumDataSize += node_size();
    this->m_nbEntry++;
    this->counter++;
    if (estimate_mem() > max_mem) {
      std::cout << "Start clean " << table.size() << std::endl;

      std::vector<size_t> last_use;
      for (auto &b : table) {
        last_use.push_back(b.last_use);
      }
      if (last_use.size() <= 2) {
        return;
      }
      std::sort(last_use.begin(), last_use.end());
      size_t median = last_use[last_use.size() / 3];
      size_t removed = 0;
      for (auto it = table.begin(); it != table.end();) {
        if (it->last_use < median) {
          delete[] it->key;
          it = table.erase(it);
          removed++;
        } else {
          it++;
        }
      }
      std::cout << "Removed: " << removed << " Cnt" << table.size()
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
  CachedBucket<T> *bucketAlreadyExist(CachedBucket<T> &_cb, size_t hashValue) {
    uint64_t *key = reinterpret_cast<uint64_t *>(hashValue);
    auto it = table.find(Node{key});
    if (it == table.end()) {
      this->m_nbNegativeHit++;
      return 0;
    } else {
      *(size_t*)(&it->last_use) = counter;
      tmp.fc = it->data;

      this->m_nbPositiveHit++;

      delete[] key;

      return &tmp;
    }
    return NULL;
  } // bucketAlreadyExist

  /**
   * Create a bucket and store it in the cache.
   *
   * @param varConnected is the set of variable.
   * @param c is the value we want to store.
   */
  inline void createAndStoreBucket(std::vector<Var> &varConnected, T &c) {
    assert(0);
  } // createBucket

  /**
   * @brief Init the hashTable.
   *
   * @param maxVar is the number of variable.
   */
  void initHashTable(unsigned maxVar) override {
    this->setInfoFormula(maxVar);
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
    assert(0);
  } // removeEntry
};
template <typename T>
std::mt19937_64 CacheListProbLRU<T>::gen = std::mt19937_64(43);

} // namespace d4
  //
