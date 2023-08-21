

#pragma once
#include <absl/container/flat_hash_map.h>

#include <array>
#include <boost/program_options.hpp>
#include <deque>
#include <functional>
#include <string_view>
#include <vector>

#include "BucketManager.hpp"
#include "CacheCleaningManager.hpp"
#include "CacheManager.hpp"
#include "CachedBucket.hpp"
#include "src/caching/cnf/BucketManagerCnf.hpp"
#include "src/hashing/HashString.hpp"
#include "src/specs/SpecManager.hpp"

namespace d4 {
namespace po = boost::program_options;

template <typename H, typename T>
H AbslHashValue(H h, const d4::CachedBucket<T> &c) {
  return H::combine(std::move(h), c.header.info1,
                    std::string_view(c.data, c.header.szData()));
}
template <typename T>
bool operator==(const d4::CachedBucket<T> &lhs,
                const d4::CachedBucket<T> &rhs) {
  return lhs.header == rhs.header && (strcmp(lhs.data, rhs.data) == 0);
}
template <typename T> class SimpleCache {
public:
  virtual bool find(std::vector<Var> &varConnected, CachedBucket<T> &bkt,
                    size_t &hash) = 0;

  virtual bool isActivated(unsigned nbVar) = 0;
  virtual void add(CachedBucket<T> &bkt, size_t hash, T val) = 0;
  virtual void printCacheInformation(std::ostream &out) = 0;
  virtual size_t getNbPositiveHit() = 0;
  virtual size_t getNbNegativeHit() = 0;
  virtual size_t usedMemory() = 0;

  static std::unique_ptr<SimpleCache<T>>
  make(po::variables_map &vm, SpecManager &specs, std::ostream &out);
};
template <typename T> class SimpleCacheConf : public SimpleCache<T> {
  std::unique_ptr<CacheManager<T>> m_impl;

public:
  SimpleCacheConf(po::variables_map &vm, SpecManager &specs,
                  std::ostream &out) {
    m_impl = std::unique_ptr<CacheManager<T>>(CacheManager<T>::makeCacheManager(
        vm, specs.getNbVariable(), &specs, out));
  }

  bool isActivated(unsigned nbVar) final { return m_impl->isActivated(nbVar); }
  bool find(std::vector<Var> &varConnected, CachedBucket<T> &bkt,
            size_t &hash) final {
    TmpEntry<T> entry = m_impl->searchInCache(varConnected);
    bkt = entry.e;
    hash = entry.hashValue;
    return entry.defined;
  }
  void add(CachedBucket<T> &bkt, size_t hash, T val) final {
    TmpEntry<T> e{bkt, (unsigned int)hash, false};
    m_impl->addInCache(e, val);
  }
  void printCacheInformation(std::ostream &out) final {
    m_impl->printCacheInformation(out);
  }
  size_t getNbPositiveHit() final { return m_impl->getNbPositiveHit(); }
  size_t getNbNegativeHit() final { return m_impl->getNbNegativeHit(); }
  size_t usedMemory() final { return m_impl->usedMemory(); }
};

template <typename T> class SimpleCacheFixed : public SimpleCache<T> {
  struct EntryInfo {
    size_t time, hits;
  };
  absl::flat_hash_map<CachedBucket<T>, EntryInfo> m_table;
  size_t max_memory = 8 * (1 << 20); // 16GB;
  BucketManager<T> *m_buckets;
  unsigned m_limitVarCached = 10000;
  size_t m_posHits = 0;
  size_t m_negHits = 0;
  size_t m_time = 0;
  size_t m_last_cull = 0;
  size_t m_max_hits = 0;

public:
  SimpleCacheFixed(po::variables_map &vm, SpecManager &specs,
                   std::ostream &out) {
    m_buckets = BucketManager<T>::makeBucketManager(vm, nullptr, specs, out);
    max_memory = vm["cache-fixed-size"].as<size_t>();
  }
  bool isActivated(unsigned nbVar) final { return nbVar <= m_limitVarCached; }
  size_t getNbPositiveHit() final { return m_posHits; }

  size_t getNbNegativeHit() final { return m_negHits; }

  size_t usedMemory() final { return m_buckets->usedMemory(); }

  bool find(std::vector<Var> &varConnected, CachedBucket<T> &bkt,
            size_t &hash) final {
    if (m_buckets->getComsumedMemory()) {
      m_buckets->reinitComsumedMemory();
      cull();
    }
    bkt = *m_buckets->collectBucket(varConnected);
    auto it = m_table.find(bkt);
    if (it == m_table.end()) {
      m_negHits++;
      return false;
    }
    m_buckets->releaseMemory(bkt.data, bkt.szData());
    bkt = it->first;
    it->second.hits++;
    if (it->second.hits > m_max_hits) {
      m_max_hits = it->second.hits;
    }
    m_posHits++;
    return true;
  } // searchInCache
    //
  void add(CachedBucket<T> &bkt, size_t _hash, T v) final {
    bkt.lockedBucket(v);
    m_table.insert({bkt, {0, m_time}});
    m_time++;
  }
  void printCacheInformation(std::ostream &out) final {
    out << "Cache Info:\n";
    out << "Hits: " << m_posHits << " Misses: " << m_negHits
        << " Mem: " << (double)m_buckets->usedMemory() / double(1 << 20) << "MB"
        << std::endl;
  }
  void decay(size_t &val) { val = val >> 2; }
  void cull() {
    size_t mem = m_buckets->usedMemory();
    double usage = double(mem / 100) / double(max_memory / 100);
    if (usage <= 1.0) {
      return;
    }
    std::cout << "Culling:" << std::endl;

    std::cout << "    usage: " << usage << std::endl;
    double cull_fac = 1 - (1 / (usage + 0.2));
    size_t to_cull = m_table.size() * cull_fac;

    std::cout << "    to_cull: " << to_cull << std::endl;
    std::array<size_t, 64> hist = {};
    size_t scale = std::min<size_t>((m_max_hits + 64 - 1) / 64, 1);
    for (auto &[k, v] : m_table) {
      size_t idx = (v.hits / scale) % 64;
      hist[idx]++;
    }
    for (int i = 1; i < 64; i++) {
      hist[i] += hist[i - 1];
    }
    size_t min_freq = 0;

    std::cout << "    min_freq: " << min_freq << std::endl;
    for (int i = 1; i < 64; i++) {
      min_freq++;
      if (hist[i] > to_cull) {
        break;
      }
    }
    min_freq = min_freq * scale;
    size_t time_diff = m_last_cull - m_time;
    size_t max_time = m_last_cull + (time_diff - time_diff / 10);
    m_last_cull = m_time;
    std::cout << "    max_time: " << min_freq << std::endl;
    std::cout << "    current_time: " << m_time << std::endl;
    for (auto i = m_table.begin(); i != m_table.end();) {
      auto it = i;
      i++;
      if (it->second.time < max_time && it->second.hits < min_freq) {
        m_buckets->releaseMemory(it->first.data, it->first.szData());
        m_table.erase(it);
      } else {
        // decay
        decay(it->second.hits);
      }
    }
    decay(m_max_hits);
  }
};
template <typename T>
std::unique_ptr<SimpleCache<T>> SimpleCache<T>::make(po::variables_map &vm,
                                                     SpecManager &specs,
                                                     std::ostream &out) {

  std::string impl = vm["cache-impl"].as<std::string>();
  if (impl == "conf") {
    return std::make_unique<SimpleCacheConf<T>>(vm, specs, out);
  } else if (impl == "fixed") {
    return std::make_unique<SimpleCacheFixed<T>>(vm, specs, out);
  }
  throw std::runtime_error("unknown caching impl ");
}
} // namespace d4
