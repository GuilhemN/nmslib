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
#ifndef _SPACE_SPARSE_JACCARD_GOLDFINGER_H_
#define _SPACE_SPARSE_JACCARD_GOLDFINGER_H_

#include <string>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>

#include "global.h"
#include "object.h"
#include "utils.h"
#include "space.h"
#include "space_vector.h"
#include "distcomp.h"
#include "space_bit_jaccard.h"

#define SPACE_SPARSE_JACCARD_GOLDFINGER              "jaccard_sparse_goldfinger"

#define LABEL_GOLDFINGER 1

namespace similarity {


/*
 * This is a Jaccard distance sparse space.
 */

template <typename dist_t>
class SpaceSparseJaccardGoldfinger : public SpaceBitJaccard<dist_t, uint32_t>
{
public:
  explicit SpaceSparseJaccardGoldfinger(const uint32_t& nbBits);
  virtual ~SpaceSparseJaccardGoldfinger() {}

  /* 
   * Space name: It will be used in result files.
   * Consider, including all the parameters where you
   * print the space name:
   */
  virtual std::string StrDesc() const
  {
    return SPACE_SPARSE_JACCARD_GOLDFINGER;
  }

  /** Standard functions to read/write/create objects */ 
  /*
   * Create an object from a (possibly binary) string.
   * If the input state pointer isn't null, we check
   * if the new vector is consistent with previously read vectors.
   * For example, when we start reading vectors,
   * we don't know the number of elements. When, we see the first
   * vector, we memorize dimensionality. If a subsequently read
   * vector has a different dimensionality, an exception will be thrown.
   */
  virtual unique_ptr<Object> CreateObjFromStr(IdType id, LabelType label, const string& s,
                                              DataFileInputState* pInpState) const;
  // Create a string representation of an object
  // The string representation may include external ID.
  // virtual string CreateStrFromObj(const Object *pObj, const string &externId) const;

protected:
  Object *CreateObjFromIds(IdType id, LabelType label, const vector<IdType> &InpVect) const;

  uint32_t nbBits_;

private:
  DISABLE_COPY_AND_ASSIGN(SpaceSparseJaccardGoldfinger);
  
};

}  // namespace similarity

#endif 
