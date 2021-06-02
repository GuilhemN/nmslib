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
#include "space/space_sparse_jaccard_goldfinger.h"

namespace similarity {

using namespace std;

static void ReadIdList(string line, LabelType& label, vector<IdType>& v)
{
  v.clear();

  label = Object::extractLabel(line);

#if 0
  ReplaceSomePunct(line); 
  stringstream str(line);

  str.exceptions(ios::badbit);

  IdType id;


  try {
    while (str >> id) {
      v.push_back(id);
    }
  } catch (const exception &e) {
    LOG(LIB_ERROR) << "Exception: " << e.what();
    LOG(LIB_FATAL) << "Failed to parse the line: '" << line << "'";
  }
#else
  if (!ReadVecDataEfficiently(line, v)) {
    PREPARE_RUNTIME_ERR(err) << "Failed to parse the line: '" << line << "'";
    LOG(LIB_ERROR) << err.stream().str();
    THROW_RUNTIME_ERR(err);
  }
#endif

  sort(v.begin(), v.end());
}

template <typename dist_t>
unique_ptr<Object>
SpaceSparseJaccardGoldfinger<dist_t>::CreateObjFromStr(IdType id, LabelType label, const string &s,
                                                       DataFileInputState *pInpStateBase) const
{
  DataFileInputStateVec*  pInpState = NULL;
  if (pInpStateBase != NULL) {
    pInpState = dynamic_cast<DataFileInputStateVec*>(pInpStateBase);
    if (NULL == pInpState) {
      PREPARE_RUNTIME_ERR(err) << "Bug: unexpected pointer type";
      THROW_RUNTIME_ERR(err);
    }
  }
  vector<IdType>  ids;
  ReadIdList(s, label, ids);
  return unique_ptr<Object>(CreateObjFromIds(id, label, ids));
}

/** End of standard functions to read/write/create objects */

template <typename dist_t>
Object *SpaceSparseJaccardGoldfinger<dist_t>::CreateObjFromIds(IdType id, LabelType label, const vector<IdType> &InpVect) const
{
  if (label == LABEL_GOLDFINGER)
    return new Object(id, label, InpVect.size() * sizeof(IdType), &InpVect[0]);

  int size = 1024;
  vector<uint32_t> sketch((size+31)/32, 0);

  for (IdType id: InpVect) {
    IdType hash = (id*5) % size;
    sketch[hash/32] |= 1 << (hash % 32);
  }

  return new Object(id, LABEL_GOLDFINGER, sketch.size() * sizeof(IdType), &sketch[0]);
};

template class SpaceSparseJaccardGoldfinger<float>;

}  // namespace similarity
