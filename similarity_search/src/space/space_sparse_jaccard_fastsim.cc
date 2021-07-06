/**
 * Non-metric Space Library
 *
 * Main developers: Bilegsaikhan Naidan, Leonid Boytsov, Yury Malkov, Ben Frederickson, David Novak
 *
 * For the complete list of contributors and further details see:
 * https://github.com/nmslib/nmslib
 *
 * Copyright (c) 2013-2018
 *
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 */
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <memory>
#include <iomanip>
#include <limits>
#include <random>

#include "object.h"
#include "utils.h"
#include "logging.h"
#include "distcomp.h"
#include "experimentconf.h"
#include "read_data.h"
#include "space/space_sparse_jaccard_fastsim.h"
#include "SpookyV2.h"

namespace similarity {

using namespace std;

template <typename dist_t>
SpaceSparseJaccardFastSim<dist_t>::SpaceSparseJaccardFastSim(const uint32_t &sketchSize)
{
  int p = 1;
  power_of_hash = 0;
  while (p < sketchSize)
  {
    p <<= 1;
    power_of_hash++;
  }

  sketchSize_ = p;

  const int seed = 3;
  static mt19937 randEngine;
  randEngine.seed(seed);

  std::uniform_int_distribution<uint64_t> drawint;

  seeds = vector<uint64_t>(2*sketchSize_);
  for (int i = 0; i < 2*sketchSize_; i++)
  {
    seeds[i] = i;
  }
}

/** End of standard functions to read/write/create objects */

template <typename dist_t>
Object *SpaceSparseJaccardFastSim<dist_t>::CreateObjFromVect(IdType id, LabelType label, const std::vector<int32_t> &InpVect) const
{
  if (label == LABEL_FASTSIM)
    return new Object(id, label, InpVect.size() * sizeof(uint32_t), &InpVect[0]);

  vector<float> sketch(sketchSize_, 2*sketchSize_+2);

  int count = 0;

  for (int l = 0; l < 2*sketchSize_; l++)
  {
    for (int32_t item : InpVect)
    {
      uint64_t hash = SpookyHash::Hash64(&item, sizeof(int32_t), seeds[l]);

      int b;
      if (l < sketchSize_)
      {
        b = (int)(hash & (sketchSize_ - 1));
      }
      else
      {
        b = l - sketchSize_;
      }

      float v = l + ((float) ((hash >> 32) & ((1L << 24) - 1)) / (float) (1L << 24));

      // double v = l + ((hash / nb_hash) & ((1L << 30) - 1)) / ((double)(1L << 30)); // We only keep the 24 most significant bits and then convert to a float between zero and one (see https://docs.oracle.com/javase/7/docs/api/java/util/Random.html#nextFloat())

      // fulfilling of the profile
      if (sketch[b] > 2*sketchSize_ + 1)
      { // first encounter
        count++;
      }
      sketch[b] = min(sketch[b], v);
    }

    if (count == sketchSize_)
    {
      break; // early stop
    }
  }

  return new Object(id, LABEL_FASTSIM, sketch.size() * sizeof(float), &sketch[0]);
};

template class SpaceSparseJaccardFastSim<float>;

}  // namespace similarity
