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

#include "object.h"
#include "utils.h"
#include "logging.h"
#include "distcomp.h"
#include "experimentconf.h"
#include "read_data.h"
#include "space/space_sparse_jaccard.h"

namespace similarity {

using namespace std;

/** End of standard functions to read/write/create objects */ 

template <typename dist_t>
Object* SpaceSparseJaccard<dist_t>::CreateObjFromIds(IdType id, LabelType label, const vector<IdType>& InpVect) const {
  return new Object(id, label, InpVect.size() * sizeof(IdType), &InpVect[0]);
};

template class SpaceSparseJaccard<float>;

}  // namespace similarity
