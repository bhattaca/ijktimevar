// procedures for ijkslice, simplex slicer for arbitrary dimension
// slices set of simplices with a hyperplane orthogonal to the x-axis

#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "ijkslice.h"
#include "ijktable.h"
#include "ijkxitIO.h"

//using namespace IJK;
//using namespace IJKTABLE;
using namespace std;

namespace IJKSLICE {


  void locate_intersected_edges
  (const int vertex_dimension, 
   const COORD_TYPE * vertex_coord, const int num_vert, 
   const int simplex_dimension,
   const int * simplex_vert, const int num_simplices, 
   const int slice_axis, const COORD_TYPE slice_coord,
   vector<int> & elist, EDGE_HASH_TABLE & edge_hash);

  void construct_sliced_simplices
  (const int vertex_dimension, 
   const COORD_TYPE * vertex_coord, const int num_vert, 
   const int simplex_dimension,
   const int * simplex_vert, const int num_simplices, 
   const int slice_axis, const COORD_TYPE slice_coord,
   vector<int> & elist, EDGE_HASH_TABLE & edge_hash,
   vector<int> & slist);

  bool read_isotable(const int dimension, const char * isotable_directory);
  void insert(vector<int> & elist, const int iv0, const int iv1);
  void insert(vector<int> & elist, EDGE_HASH_TABLE & edge_hash,
              const int v0, const int v1);
  bool find(vector<int> & elist, const int v0, const int v1, int & ie);

  ISOSURFACE_TABLE isotable(2);
  // (initial dimension is ignored)

  void slice_simplices
  (const int vertex_dimension, 
   const COORD_TYPE * vertex_coord, const int num_vert, 
   const int simplex_dimension,
   const int * simplex_vert, const int num_simplices, 
   const int slice_axis,  const COORD_TYPE slice_coord,
   const char * isotable_directory,
   int * & intersected_edge_vert, int & num_intersected_edges,
   int * & slice_simplex_vert, int & num_slice_simplices)
  // vertex_dimension = dimension of the input vertices
  //   Precondition: dimension > 0.
  // vertex_coord = array of vertex coordinates
  //   # coordinates of each vertex = dimension
  //   j'th coord of vertex i = vertex_coord[dimension*i + j]
  // num_vert = number of vertices
  // simplex_dimension = dimension of simplices
  // simplex_vert = array of simplex vertices
  //   # vertices of each simplex = simplex_dimension + 1
  //   j'th vertex of simplex i = simplex_vert[(simplex_dimension+1)*i + j]
  // num_simplices = number of simplices
  // slice_axis = dimension to be sliced. 
  //   In range [0..(vertex_dimension-1)].
  // slice_coord = coordinate of slicing hyperplane 
  //   {x_{slice_axis} = slice_coord}
  // intersected_edge_vert = array of endpoints of original simplex edges
  //   intersected by slice hyperplane {x_{slice_axis} = slice_coord}
  //   allocated and returned by slice_simplices
  //   j'th endpoint of edge i = intersected_edge_vert[2*i+j]  (j = 0 or 1)
  //   slice vertex i lies on intersected edge i
  // num_intersected_edges = number of intersected edges (returned)
  // slice_simplex_vert = array of sliced simplex vertices
  //   allocated and returned by slice_simplices
  //   # vertices of each simplex = simplex_dimension
  //   j'th vertex of simplex i = slice_simplex_vert[simplex_dimension*i + j]
  // num_slice_simplices = number of slice simplices (returned)
  {
    vector<int> intersected_e;
    vector<int> slice_s;
    EDGE_HASH_TABLE edge_hash;

    assert(vertex_dimension >= simplex_dimension && simplex_dimension > 0);
    assert(0 <= slice_axis && slice_axis < vertex_dimension);

    intersected_edge_vert = NULL;
    num_intersected_edges = 0;
    slice_simplex_vert = NULL;
    num_slice_simplices = 0;

    if (num_vert <= 0 || num_simplices <= 0)
      return;

    if (!read_isotable(simplex_dimension, isotable_directory))
      return;

    locate_intersected_edges
      (vertex_dimension, vertex_coord, num_vert, 
       simplex_dimension, simplex_vert, num_simplices, 
       slice_axis, slice_coord, intersected_e, edge_hash);

    construct_sliced_simplices
      (vertex_dimension, vertex_coord, num_vert, 
       simplex_dimension, simplex_vert, num_simplices, 
       slice_axis, slice_coord, intersected_e, edge_hash, slice_s);

    assert(intersected_e.size()%2 == 0);  
    assert(slice_s.size()%simplex_dimension == 0);

    if (intersected_e.size() == 0 || slice_s.size() == 0)
      return;

    // allocate and store intersected_edge_vert
    num_intersected_edges = intersected_e.size()/2;
    intersected_edge_vert = new int[intersected_e.size()];
    for (int i = 0; i < intersected_e.size(); i++)
      intersected_edge_vert[i] = intersected_e[i];


    // allocate and store slice_simplex_vert
    num_slice_simplices = slice_s.size()/simplex_dimension;
    slice_simplex_vert = new int[slice_s.size()];
    for (int i = 0; i < slice_s.size(); i++)
      slice_simplex_vert[i] = slice_s[i];
  }

  void locate_intersected_edges
  (const int vertex_dimension, 
   const COORD_TYPE * vertex_coord, const int num_vert, 
   const int simplex_dimension,
   const int * simplex_vert, const int num_simplices, 
   const int slice_axis, const COORD_TYPE slice_coord,
   vector<int> & elist, EDGE_HASH_TABLE & edge_hash)
  // located intersected edges and store endpoints in elist
  {
    for (int is = 0; is < num_simplices; is++) {
      for (int j = 0; j < simplex_dimension; j++)
        for (int k = j+1; k < simplex_dimension+1; k++) {
          int jv = simplex_vert[(simplex_dimension+1)*is + j];
          int kv = simplex_vert[(simplex_dimension+1)*is + k];
          COORD_TYPE jc = vertex_coord[vertex_dimension*jv+slice_axis];
          COORD_TYPE kc = vertex_coord[vertex_dimension*kv+slice_axis];
          if ((jc < slice_coord && kc >= slice_coord) ||
              (jc >= slice_coord && kc < slice_coord)) {
            insert(elist, edge_hash, jv, kv);
          };
        };
    };
  }

  void insert(vector<int> & elist, EDGE_HASH_TABLE & edge_hash,
              const int v0, const int v1)
  // insert edge (v0,v1) into elist, avoiding duplicates
  // always insert smaller vertex index first
  {
    int u0;
    int u1;

    if (v0 < v1) {
      u0 = v0;
      u1 = v1;
    }
    else {
      u0 = v1;
      u1 = v0;
    };

    VERTEX_PAIR key = std::make_pair(u0, u1);
    EDGE_HASH_TABLE::iterator edge_iter = edge_hash.find(key);
	if (edge_iter == edge_hash.end()) {
		// edge not found; insert edge
		IJKTIMEVAR::EDGE_INDEX ie = elist.size()/2;
		elist.push_back(u0);
		elist.push_back(u1);
		edge_hash.insert(EDGE_HASH_TABLE::value_type(key, ie));
	};
  }

  // find edge (v0,v1) in elist
  // return true if found; return false if not found
  // ie = index of edge found; undefined if edge not found
  bool find(const EDGE_HASH_TABLE & edge_hash, 
            const int v0, const int v1, int & ie)
  {
    bool found_flag = false;

    int u0, u1;

    if (v0 < v1) {
      u0 = v0;
      u1 = v1;
    }
    else {
      u0 = v1;
      u1 = v0;
    };

    VERTEX_PAIR key = std::make_pair(u0, u1);
    EDGE_HASH_TABLE::const_iterator edge_iter = edge_hash.find(key);
    if (edge_iter == edge_hash.end()) {
      // edge not found; insert edge
      return(false);
    }
    else {
      ie = edge_iter->second;
      return(true);
    }

  }

  void bubble_sort(int * vlist, const int numv, bool & orientation)
  {
    bool done;

    do {
      done = true;
      for (int i = 0; i+1 < numv; i++) {
        if (vlist[i] > vlist[i+1]) { 
          std::swap(vlist[i], vlist[i+1]);
          orientation = !orientation;
          done = false;
        }
      }
    } while (!done);

    return;
  }

  void construct_sliced_simplices
  (const int vertex_dimension, 
   const COORD_TYPE * vertex_coord, const int num_vert, 
   const int simplex_dimension,
   const int * simplex_vert, const int num_simplices, 
   const int slice_axis, const COORD_TYPE slice_coord,
   vector<int> & elist, EDGE_HASH_TABLE & edge_hash,
   vector<int> & slist)
  // construct (dimension-1)-simplices in hyperplane 
  //   {x_{slice_axis} = slice_coord}
  {
    assert(vertex_dimension >= simplex_dimension && simplex_dimension > 0);
    assert(0 <= slice_axis && slice_axis < vertex_dimension);

    const int num_simplex_vertices = simplex_dimension+1;
     IJK::ARRAY<int> vlist(num_simplex_vertices);
	//int vlist[num_simplex_vertices];
    const int num_iso_simplex_vertices = num_simplex_vertices-1;
    ARRAY<int> vlist2(num_iso_simplex_vertices);
	//int vlist2[num_iso_simplex_vertices];

    for (int is = 0; is < num_simplices; is++) {
      long it = 0;
      for (int j = 0; j < num_simplex_vertices; j++) {
        vlist[j] = simplex_vert[num_simplex_vertices*is + j];
      };
      bool orientation = true;
      bubble_sort(&(vlist[0]), num_simplex_vertices, orientation);

      for (int j = 0; j < num_simplex_vertices; j++) {
        int v0 = vlist[j];
        if (vertex_coord[vertex_dimension*v0+slice_axis] >= slice_coord) {
          long mask = 1L;
          mask = (mask << j);
          it = (it | mask);
        };
      }

      // store vertex indices of slice simplices in slist
      for (int js = 0; js < isotable.NumSimplices(it); js++) {

        for (int jv = 0; jv < num_iso_simplex_vertices; jv++) {
          int je = isotable.SimplexVertex(it, js, jv);
          int jv0 = isotable.Polyhedron().EdgeEndpoint(je, 0);
          int jv1 = isotable.Polyhedron().EdgeEndpoint(je, 1);
          int v0 = vlist[jv0];
          int v1 = vlist[jv1];
          int ie;
          bool found_flag = find(edge_hash, v0, v1, ie);

          if (!found_flag) {
            cerr << "Programming error." << endl;
            cerr << "  Unable to find edge (" << v0 << "," << v1 
                 << ") in list of intersected edges." << endl;
            cerr << "  Aborting." << endl;
            exit(10);
          };
          assert(0 <= ie && ie < elist.size()/2);

          // ie is index of jv'th vertex in slice simplex
          vlist2[jv] = ie;
        };

        if (num_iso_simplex_vertices > 1) {
          int ilast = num_iso_simplex_vertices-1;
          if (!orientation) { std::swap(vlist2[ilast-1], vlist2[ilast]); }
        }

        for (int i = 0; i < num_iso_simplex_vertices; i++) {
          slist.push_back(vlist2[i]);
        }
      };
    };
  }

  bool read_isotable(const int dimension, const char * isotable_directory)
  // read in isosurface table
  {
    ostringstream isotable_pathname;
    ostringstream isotable_filename;


    isotable_filename.str("");
    isotable_pathname.str("");
    isotable_filename << "iso.simplex." << dimension << "D.xit" << ends;
    isotable_pathname << isotable_directory << "/" 
                      << isotable_filename.str() << ends;

    ifstream isotable_file(isotable_pathname.str().c_str(), ios::in);
    if (!isotable_file) {
      isotable_file.clear();
      isotable_file.open(isotable_filename.str().c_str(), ios::in);
    };

    if (!isotable_file) {
      cerr << "Unable to obtain isosurface table file "
           << isotable_pathname.str() << "." << endl;
      exit(30);
    };

    try {
      IJKXIO::read_xit(isotable_file, isotable);
    }
    catch(...) {
      cerr << "Error reading file: " << isotable_filename.str() << "." << endl;
      throw;
    };

    ERROR error;
    if (!isotable.Check(error)) {
      cerr << "Warning: Data structure inconsistency in isosurface table "
           << isotable_filename.str() << "." << endl;
      error.Print(cerr);
      return(false);
    }

    return(true);
  }

}
