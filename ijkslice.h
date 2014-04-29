#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>

#include "ijktable.h"
#include "ijkxitIO.h"
#include "ijktimevar_types.h"
using namespace IJK;
using namespace IJKTIMEVAR;
using namespace IJKTABLE;
using namespace std;

namespace IJKSLICE {

void slice_simplices
(const int vertex_dimension, 
 const COORD_TYPE * vertex_coord, const int num_vert, 
 const int simplex_dimension,
 const int * simplex_vert, const int num_simplices, 
 const int slice_axis,  const COORD_TYPE slice_coord,
 const char * isotable_directory,
 int * & intersected_edge_vert, int & num_intersected_edges,
 int * & slice_simplex_vert, int & num_slice_simplices);

}
