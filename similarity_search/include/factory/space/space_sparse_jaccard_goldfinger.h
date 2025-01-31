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
#ifndef FACTORY_SPACE_SPARSE_JACCARD_GOLDFINGER_H
#define FACTORY_SPACE_SPARSE_JACCARD_GOLDFINGER_H

#include <space/space_sparse_jaccard_goldfinger.h>

namespace similarity {

/*
 * Creating functions.
 */

template <typename dist_t>
Space<dist_t>* CreateSpaceSparseJaccardGoldfinger(const AnyParams& AllParams) {
  AnyParamManager pmgr(AllParams);

  string nbBits;

  pmgr.GetParamRequired("nb_bits", nbBits);
  pmgr.CheckUnused();

  return new SpaceSparseJaccardGoldfinger<dist_t>(stoi(nbBits));
}

/*
 * End of creating functions.
 */

}

#endif
