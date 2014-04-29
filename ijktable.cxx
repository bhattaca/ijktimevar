/// \file ijktable.cxx
/// Class containing a table of isosurface patches in a given polyhedron.
/// Version 0.3.1

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012, 2011, 2007, 2006, 2003, 2001 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <ctype.h>
#include <limits>
#include <limits.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <set>
#include <vector>

#include "ijktable.h"

using namespace IJK;
using namespace IJKTABLE;
using namespace std;

#ifndef LONG_BIT

#define LONG_BIT (CHAR_BIT * sizeof(long))

#endif

// local namespace
namespace {
  static const char * standard_encoding_name[] =    
    { "BINARY", "BASE3", "NONSTANDARD" };
}

//**************************************************
// ISOSURFACE_VERTEX
//**************************************************

ISOSURFACE_VERTEX::ISOSURFACE_VERTEX()
  // constructor
{ 
  coord = NULL; 
  num_coord = 0; 
  is_label_set = false; 
}

ISOSURFACE_VERTEX::~ISOSURFACE_VERTEX()
  // destructor
{
  if (coord != NULL) { 
    delete [] coord;
    coord = NULL;
  };

  num_coord = 0;
  is_label_set = false;
}

void ISOSURFACE_VERTEX::SetNumCoord(const int numc)
{
  if (coord != NULL)
    delete [] coord;
  coord = NULL;

  num_coord = numc;
  coord = new COORD_TYPE[numc];
}

//**************************************************
// ISOSURFACE_TABLE
//**************************************************

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::ISOSURFACE_TABLE_ENTRY()
  // constructor
{
  num_simplices = 0;
  simplex_vertex_list = NULL;
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::~ISOSURFACE_TABLE_ENTRY()
  // destructor
{
  FreeAll();
}

bool ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_simplices < 0) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Number of simplices in isosurface table entry must be non-negative.");
    return(false);
  }

  if (num_simplices > 0 && simplex_vertex_list == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for simplex vertex list not allocated.");
    return(false);
  }

  return(true);
}

void ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::FreeAll()
  // free all memory
{
  delete [] simplex_vertex_list;
  simplex_vertex_list = NULL;
  num_simplices = 0;
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE():polyhedron(3)
  // default constructor. dimension = 3
{
  Init(3, 2);
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE(const int d):polyhedron(d)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d, d-1);
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE
(const int dimension, const int simplex_dimension):polyhedron(dimension)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(dimension, simplex_dimension);
}

void ISOSURFACE_TABLE::Init(const int dimension, const int simplex_dimension)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  const char * procname = "ISOSURFACE_TABLE::Init";

  this->simplex_dimension = simplex_dimension;

  max_num_vertices = 20;
  // Note: Even tables for polyhedron of this size are probably impossible 
  //   to compute/store

  num_isosurface_vertices = 0;
  isosurface_vertex = NULL;

  encoding = NONSTANDARD;
  encoding_name = StandardEncodingName(NONSTANDARD);

  num_table_entries = 0;
  entry = NULL;
  is_table_allocated = false;
  if (!CheckDimension())
    throw PROCEDURE_ERROR(procname, "Illegal polyhedron dimension.");
}

ISOSURFACE_TABLE::~ISOSURFACE_TABLE()
  // destructor
{
  FreeAll();
}

std::string ISOSURFACE_TABLE::StandardEncodingName
(const ENCODING encoding)
{
  return(standard_encoding_name[encoding]);
}

void ISOSURFACE_TABLE::SetEncoding(const ENCODING encoding)
{
  this->encoding = encoding;
  encoding_name = StandardEncodingName(encoding);
}

void ISOSURFACE_TABLE::SetNonstandardEncoding(const std::string & name)
{
  this->encoding = NONSTANDARD;
  encoding_name = name;
}

void ISOSURFACE_TABLE::SetNumIsosurfaceVertices(const int num_vertices)
{
  if (isosurface_vertex != NULL)
    delete [] isosurface_vertex;
  isosurface_vertex = NULL;

  num_isosurface_vertices = num_vertices;
  isosurface_vertex = new ISOSURFACE_VERTEX[num_vertices];
}

void ISOSURFACE_TABLE::CheckIsoVerticesAlloc
(const char * procname, const int vstart, const int numv)
  // check allocation of array isosurface_vertices
  // procname = calling procedure name, for error messages
  // vstart = first vertex
  // numv = number of vertices required
{
  if (numv == 0) return;

  if (isosurface_vertex == NULL) {
    throw PROCEDURE_ERROR
      (procname, "Set number of isosurface vertices before storing vertices.");
  }

  if (numv+vstart > NumIsosurfaceVertices()) {
    throw PROCEDURE_ERROR
      (procname, "Illegal isosurface vertex index.");
  }
}

void ISOSURFACE_TABLE::StorePolyVerticesAsIsoVertices(const int vstart)
  // store polyhedron vertices as isosurface vertices
  // store polyhedron vertices starting at isosurface vertex vstart
{
  const int num_polyv = Polyhedron().NumVertices();
  const char * procname = "ISOSURFACE_TABLE::StorePolyVerticesAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polyv);

  for (int iv = 0; iv < num_polyv; iv++) {
    SetIsoVertexType(iv+vstart, ISOSURFACE_VERTEX::VERTEX);
    SetIsoVertexFace(iv+vstart, iv);
  }
}

void ISOSURFACE_TABLE::StorePolyEdgesAsIsoVertices(const int vstart)
  // store polyhedron edges as isosurface vertices
  // store polyhedron edges starting at isosurface vertex vstart
{
  const int num_polye = Polyhedron().NumEdges();
  const char * procname = "ISOSURFACE_TABLE::StorePolyEdgesAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polye);

  for (int ie = 0; ie < num_polye; ie++) {
    SetIsoVertexType(ie+vstart, ISOSURFACE_VERTEX::EDGE);
    SetIsoVertexFace(ie+vstart, ie);
  }
}

void ISOSURFACE_TABLE::StorePolyFacetsAsIsoVertices(const int vstart)
  // store polyhedron facets as isosurface vertices
  // store polyhedron facets starting at isosurface vertex vstart
{
  const int num_polyf = Polyhedron().NumFacets();
  const char * procname = "ISOSURFACE_TABLE::StorePolyFacetsAsIsoVertices";

  CheckIsoVerticesAlloc(procname, vstart, num_polyf);

  for (int jf = 0; jf < num_polyf; jf++) {
    SetIsoVertexType(jf+vstart, ISOSURFACE_VERTEX::FACET);
    SetIsoVertexFace(jf+vstart, jf);
  }
}

void ISOSURFACE_TABLE::SetNumTableEntries(const int num_table_entries)
  // allocate table
{
  const char * procname = "ISOSURFACE_TABLE::SetNumTableEntries";

  if (entry != NULL) delete [] entry;
  entry = NULL;
  this->num_table_entries = 0;
  is_table_allocated = false;

  entry = new ISOSURFACE_TABLE_ENTRY[num_table_entries];
  if (entry == NULL && num_table_entries > 0)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for isosurface table.");

  this->num_table_entries = num_table_entries;
  is_table_allocated = true;
}


void ISOSURFACE_TABLE::SetNumSimplices(const TABLE_INDEX it, const int nums)
  // set number of simplices in table entry it
  // it = table entry
  // nums = number of simplices
{
  const char * procname = "ISOSURFACE_TABLE::SetNumSimplices";

  if (!IsTableAllocated() || entry == NULL) {
    throw PROCEDURE_ERROR
      (procname, "Table must be allocated before entering table entries.");
  };

  if (it < 0 || it >= NumTableEntries())
    throw PROCEDURE_ERROR(procname, "Illegal table index.");
  if (nums < 0)
    throw PROCEDURE_ERROR
      (procname, "Number of simplices must be non-negative.");

  entry[it].num_simplices = 0;
  delete entry[it].simplex_vertex_list;
  entry[it].simplex_vertex_list = NULL;

  if (nums > 0)
    entry[it].simplex_vertex_list = 
      new ISOSURFACE_VERTEX_INDEX[nums*NumVerticesPerSimplex()];

  entry[it].num_simplices = nums;
}


void ISOSURFACE_TABLE::SetSimplexVertex
(const TABLE_INDEX it, const int is, const int k, 
 const ISOSURFACE_VERTEX_INDEX isov)
  // set simplex vertex
  // it = index table entry.  In range [0..NumTableEntries()-1].
  // is = index simplex.  
  // k = k'th simplex vertex.  In range [0..NumVerticesPerSimplex()-1].
  // isov = index of isosurface vertex
{
  entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+k] = isov;
}

bool ISOSURFACE_TABLE::CheckDimension(const int d) const
  // check dimension
{
  if (d < 1)
    return(false);
  else
    return(true);
}

bool ISOSURFACE_TABLE::CheckTable(ERROR & error_msg) const
  // check table
{
  if (polyhedron.NumVertices() >= LONG_BIT) {
    error_msg.ClearAll();
    error_msg.AddMessage("Too many polyhedron vertices");
    return(false);
  }

  if (polyhedron.NumVertices() > MaxNumVertices()) {
    error_msg.ClearAll();
    error_msg.AddMessage("Too many polyhedron vertices");
    return(false);
  }

  if (polyhedron.NumVertices() < 1) {
    error_msg.ClearAll();
    error_msg.AddMessage("Polyhedron must have at least one vertex.");
    return(false);
  }

  if (entry == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for isosurface table not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!entry[it].Check(error_msg)) {
      error_msg.AddMessage("Error detected at isosurface table entry ", 
                           it, ".");
      return(false);
    }

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int is = 0; is < NumSimplices(jt); is++)
      for (int iv = 0; iv < NumVerticesPerSimplex(); iv++) {
        int iso_v = int(SimplexVertex(jt, is, iv));
        if (iso_v < 0 || iso_v >= NumIsosurfaceVertices()) {
          error_msg.ClearAll();
          error_msg.AddMessage
            ("Illegal isosurface vertex ", iso_v, " in isosurface table entry ", jt, "."); 
          return(false);
        }
      };

  return(true);
}

bool ISOSURFACE_TABLE::Check(ERROR & error_msg) const
{
  if (!Polyhedron().Check(error_msg)) return(false);
  if (!CheckTable(error_msg)) return(false);
  return(true);
}

void ISOSURFACE_TABLE::FreeAll()
  // free all memory
{
  if (entry != NULL) {
    for (int i = 0; i < num_table_entries; i++)
      entry[i].FreeAll();
    delete [] entry;
    entry = NULL;
  };
  num_table_entries = 0;
  is_table_allocated = false;

  polyhedron.FreeAll();

  delete [] isosurface_vertex;
  isosurface_vertex = NULL;
  num_isosurface_vertices = 0;
}


//**************************************************
// ISOSURFACE_EDGE_TABLE
//**************************************************

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::ISOSURFACE_EDGE_TABLE_ENTRY()
  // constructor
{
  num_edges = 0;
  edge_endpoint_list = NULL;
}

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::~ISOSURFACE_EDGE_TABLE_ENTRY()
  // destructor
{
  FreeAll();
}

bool ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::Check
(ERROR & error_msg) const
{
  if (num_edges < 0) {
    error_msg.ClearAll();
    error_msg.AddMessage
      ("Number of edges in isosurface table entry must be non-negative.");
    return(false);
  }

  if (num_edges > 0 && edge_endpoint_list == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for edge endpoint list not allocated.");
    return(false);
  }

  return(true);
}

void ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE_ENTRY::FreeAll()
  // free all memory
{
  delete [] edge_endpoint_list;
  edge_endpoint_list = NULL;
  num_edges = 0;
}

ISOSURFACE_EDGE_TABLE::ISOSURFACE_EDGE_TABLE(const int d):ISOSURFACE_TABLE(d)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  Init(d);
}

void ISOSURFACE_EDGE_TABLE::Init(const int d)
  // constructor
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  edge_entry = NULL;
}

ISOSURFACE_EDGE_TABLE::~ISOSURFACE_EDGE_TABLE()
  // destructor
{
  FreeAll();
}

void ISOSURFACE_EDGE_TABLE::SetNumTableEntries(const int num_table_entries)
  // allocate table
{
  const char * procname = "ISOSURFACE_EDGE_TABLE::SetNumTableEntries";

  ISOSURFACE_TABLE::SetNumTableEntries(num_table_entries);

  edge_entry = new ISOSURFACE_EDGE_TABLE_ENTRY[num_table_entries];
  if (edge_entry == NULL)
    throw PROCEDURE_ERROR
      (procname, "Unable to allocate memory for isosurface edge table.");
}

bool ISOSURFACE_EDGE_TABLE::CheckTable(ERROR & error_msg) const
  // check table
{
  if (!ISOSURFACE_TABLE::CheckTable(error_msg)) return(false);

  if (edge_entry == NULL) {
    error_msg.ClearAll();
    error_msg.AddMessage("Memory for isosurface table edge entries not allocated.");
    return(false);
  }

  for (int it = 0; it < NumTableEntries(); it++)
    if (!edge_entry[it].Check(error_msg)) {
      error_msg.AddMessage("Error detected at isosurface table entry ", 
                           it, ".");
      return(false);
    }

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int ie = 0; ie < NumEdges(jt); ie++)
      for (int iend = 0; iend < 2; iend++) {
        int iso_v = int(EdgeEndpoint(jt, ie, iend));
        if (iso_v < 0 || iso_v >= NumIsosurfaceVertices()) {
          error_msg.ClearAll();
          error_msg.AddMessage
            ("Illegal isosurface vertex ", iso_v, " in isosurface table edge entry ", jt, "."); 
          return(false);
        }
      };

  return(true);
}

bool ISOSURFACE_EDGE_TABLE::Check(ERROR & error_msg) const
{
  if (!Polyhedron().Check(error_msg)) return(false);
  if (!CheckTable(error_msg)) return(false);
  return(true);
}

void ISOSURFACE_EDGE_TABLE::FreeAll()
  // free all memory
{
  if (edge_entry != NULL) {
    for (int i = 0; i < num_table_entries; i++) {
      edge_entry[i].FreeAll();
    }
    delete [] edge_entry;
    edge_entry = NULL;  
  };

  ISOSURFACE_TABLE::FreeAll();
}

namespace {
  class EDGE_CONTAINER { 
  public:
    EDGE_INDEX v[2]; 

    EDGE_CONTAINER(){};
    bool operator < (const EDGE_CONTAINER & e) const {
      if (v[0] < e.v[0]) { return(true); }
      else if (v[0] == e.v[0] && v[1] < e.v[1]) { return(true); }
      else {return(false); }
    }
  };

}

void ISOSURFACE_EDGE_TABLE::GenEdgeLists()
  // for each table entry, generate edge lists from simplex lists
{
  PROCEDURE_ERROR error("ISOSURFACE_EDGE_TABLE::GenEdgeLists");
  IJK::ARRAY<EDGE_INDEX> vlist(Dimension());
  typedef set<EDGE_CONTAINER> EDGE_SET;
  EDGE_SET eset;

  if (!IsTableAllocated()) {
    error.AddMessage("Programming error: Isosurface table not allocated.");
    throw error;
  };

  for (int it = 0; it < NumTableEntries(); it++) {
    eset.clear();

    for (int is = 0; is < NumSimplices(it); is++) {

      // get simplex vertices
      for (int iv = 0; iv < Dimension(); iv++) {
        vlist[iv] = SimplexVertex(it, is, iv);
      }
      sort(vlist.Ptr(), vlist.Ptr()+Dimension());

      // store simplex edges
      for (int i0 = 0; i0 < Dimension(); i0++)
        for (int i1 = i0+1; i1 < Dimension(); i1++) {
          EDGE_CONTAINER e;
          e.v[0] = vlist[i0];
          e.v[1] = vlist[i1];

          eset.insert(e);
        }
    }

    edge_entry[it].FreeAll();
    if (eset.size() > 0) {
      edge_entry[it].edge_endpoint_list = new EDGE_INDEX[2*eset.size()];
      if (edge_entry[it].edge_endpoint_list == NULL) {
        error.AddMessage("Unable to allocate memory for edge list in table entry ", it, ".");
        throw error;
      };
      edge_entry[it].num_edges = eset.size();

      int ie = 0;
      for (EDGE_SET::iterator edge_iter = eset.begin(); 
           edge_iter != eset.end(); 
           edge_iter++) {
        edge_entry[it].edge_endpoint_list[2*ie] = edge_iter->v[0];
        edge_entry[it].edge_endpoint_list[2*ie+1] = edge_iter->v[1];
        ie++;
      }
    }
  }

}

//**************************************************
// UTILITY FUNCTIONS
//**************************************************

unsigned long IJKTABLE::calculate_num_entries
(const int num_vert, const int num_colors)
  // calculate num table entries = (num_colors)^(num_vert)
{
  IJK::PROCEDURE_ERROR error("calculate_num_entries");

  unsigned long num_table_entries = 1;
  IJK::int_power(num_colors, num_vert, num_table_entries, error);

  return(num_table_entries);
}

