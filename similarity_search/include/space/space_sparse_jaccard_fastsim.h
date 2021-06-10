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
#ifndef _SPACE_SPARSE_JACCARD_FASTSIM_H_
#define _SPACE_SPARSE_JACCARD_FASTSIM_H_

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

#define SPACE_SPARSE_JACCARD_FASTSIM "jaccard_sparse_fastsim"
#define LABEL_FASTSIM 1

namespace similarity {


/*
 * This is a Jaccard distance sparse space.
 */

template <typename dist_t>
class SpaceSparseJaccardFastSim : public VectorSpaceSimpleStorage<dist_t, IdType>
{
public:
  explicit SpaceSparseJaccardFastSim(const uint32_t &sketchSize);
  virtual ~SpaceSparseJaccardFastSim() {}

  /* 
   * Space name: It will be used in result files.
   * Consider, including all the parameters where you
   * print the space name:
   */
  virtual std::string StrDesc() const {
    return SPACE_SPARSE_JACCARD_FASTSIM;
  }

  /*
   * CreateDenseVectFromObj and GetElemQty() are only needed, if
   * one wants to use methods with random projections.
   */
  virtual void CreateDenseVectFromObj(const Object* obj, dist_t* pVect,
                                   size_t nElem) const {
    throw runtime_error("Cannot create vector for the space: " + StrDesc());
  }
  unsigned ComputeOverlap(const Object* pObj1, const Object* pObj2) const {
    const IdType* p1 = reinterpret_cast<const IdType*>(pObj1->data());
    const IdType* p2 = reinterpret_cast<const IdType*>(pObj2->data());

    int c = 0;
    for (int i = 0; i < sketchSize_; i++) {
      c += p1[i] == p2[i];
    }

    return c;
  }
  unsigned ComputeOverlap(const Object* pObj1, const Object* pObj2, const Object* pObj3) const {
    const IdType* p1 = reinterpret_cast<const IdType*>(pObj1->data());
    const IdType* p2 = reinterpret_cast<const IdType*>(pObj2->data());
    const IdType* p3 = reinterpret_cast<const IdType*>(pObj3->data());

    int c = 0;
    for (int i = 0; i < sketchSize_; i++) {
      c += p1[i] == p2[i] && p2[i] == p3[i];
    }

    return c;
  }

  virtual Object *CreateObjFromVect(IdType id, LabelType label, const std::vector<int32_t> &InpVect) const override;

protected:
  int sketchSize_;
  vector<uint64_t> h11;
  vector<uint32_t> h12;
  vector<uint64_t> h2;
  // The number of blocks we use for the hashing
  int c = 4;
  int d = 4;

  int power_alphabetc;
  int alphabet_sizec;
  int power_alphabetd;
  int alphabet_sized;
  int power_of_hash;
  int maskc;
  int maskd;

  /*
  * This function should always be protected.
  * Only children and friends of the Space class
  * should be able to access it.
  */
  virtual dist_t HiddenDistance(const Object* pObj1, const Object* pObj2) const {
    const IdType* p1 = reinterpret_cast<const IdType*>(pObj1->data());
    const IdType* p2 = reinterpret_cast<const IdType*>(pObj2->data());

    int c = 0;
    for (int i = 0; i < sketchSize_; i++) {
      c += p1[i] == p2[i];
    }

    return (dist_t) c / sketchSize_;
  }

private:
  /*
   * One should forbid making copies of the Space object
   */
  DISABLE_COPY_AND_ASSIGN(SpaceSparseJaccardFastSim);
};

}  // namespace similarity

#endif 
