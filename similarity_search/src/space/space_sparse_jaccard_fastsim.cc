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

namespace similarity {

using namespace std;

template <typename dist_t>
SpaceSparseJaccardFastSim<dist_t>::SpaceSparseJaccardFastSim(const uint32_t &sketchSize)
{
  int p = 1;
  power_of_hash = 1;
  while (p < sketchSize)
  {
    p <<= 1;
    power_of_hash++;
  }

  sketchSize_ = p;

  // INITIALISATION OF THE TABULATION HASHING, SEE https://arxiv.org/pdf/1411.7191.pdf
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<uint32_t> drawint;
  std::uniform_int_distribution<uint64_t> drawlong;

  // We consider 20 bits are sufficient to store items' id
  power_alphabetc = 20 / c;
  alphabet_sizec = 1 << power_alphabetc; // 2 ** power_alphabetc

  power_alphabetd = 32 / d;
  alphabet_sized = 1 << power_alphabetd; // 2 ** (32/d)

  h11 = vector<uint64_t>(sketchSize_ * alphabet_sizec * c);
  h12 = vector<uint32_t>(sketchSize_ * alphabet_sizec * c);
  for (int i = 0; i < sketchSize_ * alphabet_sizec * c; i++)
  {
    h11[i] = drawlong(gen);
    h12[i] = drawint(gen);
  }

  h2 = vector<uint64_t>(alphabet_sized * d);
  for (int i = 0; i < alphabet_sized * d; i++)
  {
    h2[i] = drawlong(gen);
  }

  // COMPUTATION OF THE PROFILES
  maskc = (1 << (power_alphabetc + 1)) - 1;
  maskd = (1 << (power_alphabetd + 1)) - 1;
}

/** End of standard functions to read/write/create objects */

template <typename dist_t>
Object *SpaceSparseJaccardFastSim<dist_t>::CreateObjFromVect(IdType id, LabelType label, const std::vector<int32_t> &InpVect) const
{
  if (label == LABEL_FASTSIM)
    return new Object(id, label, InpVect.size() * sizeof(int32_t), &InpVect[0]);

  vector<int32_t> sketch(sketchSize_, 0);

  int count = 0;

  for (int l = 0; l < 2 * sketchSize_; l++)
  {
    for (uint32_t item : InpVect)
    {
      // hashing
      int val = item;
      long v11 = 0;
      int v12 = 0;
      long v2 = 0;

      // We compute h1
      for (int m = 0; m < c; m++)
      {
        v11 ^= h11[(alphabet_sizec * (l + sketchSize_ * m) + val) & maskc];
        v12 ^= h12[(alphabet_sizec * (l + sketchSize_ * m) + val) & maskc];

        val >>= power_alphabetc;
      }

      // We compute h2
      for (int m = 0; m < d; m++)
      {
        v2 ^= h2[(m * alphabet_sized + v12) & maskd];
        v12 >>= power_alphabetd;
      }

      long hash = v11 ^ v2;

      int b;
      if (l < sketchSize_)
      {
        b = (int)(hash & (sketchSize_ - 1));
      }
      else
      {
        b = l - sketchSize_;
      }

      int v = (int)(hash >> (power_of_hash + 1));
      v |= (l << (31 - power_of_hash));

      // double v = l + ((hash / nb_hash) & ((1L << 30) - 1)) / ((double)(1L << 30)); // We only keep the 24 most significant bits and then convert to a float between zero and one (see https://docs.oracle.com/javase/7/docs/api/java/util/Random.html#nextFloat())

      // fulfilling of the profile
      if (sketch[b] > 2 * sketchSize_)
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

  return new Object(id, LABEL_FASTSIM, sketch.size() * sizeof(int32_t), &sketch[0]);
};

template class SpaceSparseJaccardFastSim<float>;

}  // namespace similarity
