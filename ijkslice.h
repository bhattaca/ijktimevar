// procedures for ijkslice, simplex slicer for arbitrary dimension
// slices set of simplices with a hyperplane orthogonal to the x-axis


#include "ijktimevar_types.h"
#include <unordered_map>

using namespace IJK;
using namespace IJKTIMEVAR;
using namespace IJKTABLE;
using namespace std;

namespace IJKSLICE {

  // types
  //typedef float COORD_TYPE;
  //typedef int VERTEX_INDEX;
  //
  //typedef int NUM_TYPE;
  typedef std::pair<VERTEX_INDEX, VERTEX_INDEX> VERTEX_PAIR;

  // **************************************************
  // EDGE HASH TABLE
  // **************************************************

  struct HASH_VERTEX_PAIR {

    std::hash<VERTEX_INDEX> hash_func;

    HASH_VERTEX_PAIR(){};

    size_t operator() (const VERTEX_PAIR & key) const
    {
      return(hash_func(key.first+key.second));
    };
  };

  typedef std::unordered_map<VERTEX_PAIR, IJKTIMEVAR::EDGE_INDEX, HASH_VERTEX_PAIR> 
  EDGE_HASH_TABLE;

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

