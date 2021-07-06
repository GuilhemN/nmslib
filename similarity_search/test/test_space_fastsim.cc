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
#include <memory>
#include <vector>

#include "space/space_sparse_jaccard_fastsim.h"
#include "bunit.h"
#include "testdataset.h"
#include "object.h"

namespace similarity {

#define FLOAT_TYPE float

TEST(FastSim) {
  vector<float> results;

  for (int i = 0; i < 10000; i++) {
    SpaceSparseJaccardFastSim<FLOAT_TYPE>* space = new SpaceSparseJaccardFastSim<FLOAT_TYPE>(128);

    vector<int32_t> data1 = {1, 2};
    Object *obj1 = space->CreateObjFromVect(0, 0, data1);

    vector<int32_t> data2 = {2, 3};
    Object *obj2 = space->CreateObjFromVect(0, 0, data2);

    results.push_back(space->Distance(obj1, obj2));
  }

  string filename("tmp.json");
  fstream file_out;

  file_out.open(filename, std::ios_base::out);
  if (!file_out.is_open())
  {
    cout << "failed to open " << filename << '\n';
  }
  else
  {

    file_out << "[" << endl;
    bool first = true;
    for (float result: results) {
      if (!first)
        file_out << "," << endl;
      first = false;

      file_out << result;
    }
    file_out << endl << "]";
    cout << "Done Writing!" << endl;
  }
}


}  // namespace similarity

