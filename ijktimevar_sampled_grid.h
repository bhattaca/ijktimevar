
#include "ijktimevar_datastruct.h"
#include <iostream>
#include <vector>
#include "ijktimevarIO.h"
using namespace std;
using namespace IJKTIMEVAR;

//compute the sampled grid, to level 'l', 
bool initialize_sampled_grid ( 
	const int l,
	const int dimension,
	const MC_SCALAR_GRID_BASE & scalar_grid,
	vector <int> &vec_time,
	MC_SCALAR_GRID & sampled_grid); 

void compute_time_steps(
	const int twoK,
	const int time_steps,
	vector<int> &vec_time
	);

void determine_copy_region_size
	(const int dimension,
	const MC_SCALAR_GRID_BASE & full_scalar_grid,
	const int time_axis,
	IJK::ARRAY<int> & copy_region_size);

void populate_grid
	(
	const vector<int> &vec_time,
	const int time_axis,
	const int dimension,
	IJK::ARRAY<int> & copy_region_size,
	const MC_SCALAR_GRID_BASE & full_scalar_grid, 
	MC_SCALAR_GRID & sampled_grid
	);

void rescale_time_axis
	(
	const int k, //level
	const int dimension,
	const vector<int> vec_time,
	MC_ISOSURFACE & mc_isosurface
	);

void copy2_sampled_grid ( 
	const MC_SCALAR_GRID_BASE & scalar_grid, 
	const VERTEX_INDEX vstart_scalar_grid,
	const VERTEX_INDEX vstart_sampled_grid,
	IJK::ARRAY<int> sampled_axis_size,
	//returns
	MC_SCALAR_GRID_BASE & sampled_grid);

void evaluate_multi_slice (
	const int ta, 
	const int tc,
	const SCALAR_TYPE & isovalue,
	POLY_ISOTABLE & isotable,
	const string & isotable_directory,
	const MC_SCALAR_GRID_BASE & scalar_grid, 
	vector<int> rqd_tsteps,
	IO_INFO &io_info
	);
