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
#ifndef FACTORY_SPACE_SPARSE_JACCARD_FASTSIM_H
#define FACTORY_SPACE_SPARSE_JACCARD_FASTSIM_H

#include <space/space_sparse_jaccard_fastsim.h>

namespace similarity {

/*
 * Creating functions.
 */

template <typename dist_t>
Space<dist_t>* CreateSpaceSparseJaccardFastSim(const AnyParams& AllParams) {
  AnyParamManager pmgr(AllParams);

  string sketchSize;

  pmgr.GetParamRequired("sketch_size", sketchSize);
  pmgr.CheckUnused();

  return new SpaceSparseJaccardFastSim<dist_t>(stoi(sketchSize));
}

/*
 * End of creating functions.
 */

}

#endif
