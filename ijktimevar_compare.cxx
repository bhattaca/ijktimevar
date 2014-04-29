#include"ijktimevar_compare.h"
#include "ijkmcube.h"
#include "ijkIO.txx"
#include "ijktimevarIO.h"
using namespace std;


void slice_timestep
	(
	const int dimension,
	const VERTEX_INDEX numv, 
	const VERTEX_INDEX nums,
	const MC_ISOSURFACE &mc_full_isosurface,
	const float time_step,
	COORD_TYPE * temp_coord, 
	int * temp_slice_simplex_vert,
	VERTEX_INDEX & temp_num_slice_vert,
	VERTEX_INDEX & temp_num_slice_simplices,
	const IO_INFO &io_info
	)
{

}


bool compute_dist 
	(float & dist)
{
  int rand_num = rand();
  if (rand_num%3 == 0)
  {
	  dist = 0.5;
	  return true;
  }
  else
	  return false;
}
void compare_isosurfaces(
	const int dimension,
	const int k,
	const SCALAR_TYPE isovalue,
	const VERTEX_INDEX numv_sampled, 
	const VERTEX_INDEX nums_sampled,
	const MC_ISOSURFACE &mc_sampled_isosurface,

	const VERTEX_INDEX numv_full, 
	const VERTEX_INDEX nums_full,
	const MC_ISOSURFACE &mc_full_isosurface,

	const MC_SCALAR_GRID_BASE & full_scalar_grid, 
	const MC_SCALAR_GRID_BASE & sampled_grid,
	vector<int> &  required_time_steps,
	const IO_INFO &io_info
	)
{
	cout <<"function compare isosurfaces"<< endl;
	VERTEX_INDEX step = 0;
	const VERTEX_INDEX max_time_step = full_scalar_grid.AxisSize(dimension-1);	
	VERTEX_INDEX time_step = pow(2,k)*step + pow(2,k-1);
	while (time_step < max_time_step)
	{
		cout <<"comparing slice with exact "<< time_step << endl;
		cout <<"if required "<< pow(2,k)*step <<", "<< pow(2,k)*step + pow(2,k) 
			<<" , "<< pow(2,k)*step + pow(2,k-1)<< endl;

		const COORD_TYPE slice_coord = float(time_step);
		const int slice_axis = dimension-1;
		int * intersected_edge_vert = NULL;
		int num_intersected_edges = 0;
		int * slice_simplex_vert_sampled = NULL;
		int num_slice_simplices_sampled = 0;
		slice_simplices
			(dimension, &(mc_sampled_isosurface.vertex_coord[0]), numv_sampled, dimension-1, 
			&(mc_sampled_isosurface.simplex_vert[0]), nums_sampled,
			slice_axis, slice_coord, io_info.isotable_directory.c_str(),
			intersected_edge_vert, num_intersected_edges,
			slice_simplex_vert_sampled, num_slice_simplices_sampled);

		int num_slice_vert_sampled = num_intersected_edges;
		COORD_TYPE * coord_sampled = new COORD_TYPE[(dimension-1)*num_slice_vert_sampled];
		int * coord0 = new int[dimension-1];
		int * coord1 = new int[dimension-1];

		for (int iv = 0; iv < num_slice_vert_sampled; iv++) {
			int v0 = intersected_edge_vert[2*iv];
			int v1 = intersected_edge_vert[2*iv+1];
			float c0 = mc_sampled_isosurface.vertex_coord[v0*dimension+slice_axis];
			float c1 = mc_sampled_isosurface.vertex_coord[v1*dimension+slice_axis];
			float s0 = (c1 - slice_coord)/(c1-c0);
			float s1 = (c0 - slice_coord)/(c0-c1);
			int d2 = 0;
			for (int d = 0; d < dimension; d++) {
				if (d != slice_axis) {
					float c = s0*mc_sampled_isosurface.vertex_coord[v0*dimension+d]+
						s1*mc_sampled_isosurface.vertex_coord[v1*dimension+d];
					coord_sampled[iv*(dimension-1) + d2] = c;
					d2++;
				};
			};
		}
		delete [] coord0;
		delete [] coord1;
		delete [] intersected_edge_vert;
		///////////////////////
		num_intersected_edges = 0;
		int * slice_simplex_vert_full = NULL;
		int num_slice_simplices_full = 0;
		slice_simplices
			(dimension, &(mc_full_isosurface.vertex_coord[0]), numv_full, dimension-1, 
			&(mc_full_isosurface.simplex_vert[0]), nums_full,
			slice_axis, slice_coord, io_info.isotable_directory.c_str(),
			intersected_edge_vert, num_intersected_edges,
			slice_simplex_vert_full, num_slice_simplices_full);

		int num_slice_vert_full = num_intersected_edges;
		COORD_TYPE * coord_full = new COORD_TYPE[(dimension-1)*num_slice_vert_full];
		coord0 = new int[dimension-1];
		coord1 = new int[dimension-1];

		for (int iv = 0; iv < num_slice_vert_full; iv++) {
			int v0 = intersected_edge_vert[2*iv];
			int v1 = intersected_edge_vert[2*iv+1];
			float c0 = mc_full_isosurface.vertex_coord[v0*dimension+slice_axis];
			float c1 = mc_full_isosurface.vertex_coord[v1*dimension+slice_axis];
			float s0 = (c1 - slice_coord)/(c1-c0);
			float s1 = (c0 - slice_coord)/(c0-c1);
			int d2 = 0;
			for (int d = 0; d < dimension; d++) {
				if (d != slice_axis) {
					float c = s0*mc_full_isosurface.vertex_coord[v0*dimension+d]+
						s1*mc_full_isosurface.vertex_coord[v1*dimension+d];
					coord_full[iv*(dimension-1) + d2] = c;
					d2++;
				};
			};
		}
		delete [] coord0;
		delete [] coord1;
		delete [] intersected_edge_vert;

		//DEBUG
		/*ostream * out = NULL;
		cout <<"output mesh"<< endl;
		out = new ofstream("sampled.off", ios::out);
		ijkoutOFF(*out, dimension-1, coord_sampled, num_slice_vert_sampled, 
			slice_simplex_vert_sampled, num_slice_simplices_sampled);

		ostream * out1 = NULL;
		cout <<"output mesh"<< endl;
		out1 = new ofstream("full.off", ios::out);
		ijkoutOFF(*out1, dimension-1, coord_full, num_slice_vert_full, 
			slice_simplex_vert_full, num_slice_simplices_full);*/
	

		//distance between the exact and sliced version.
		float dist = 0.0;
		//Returns true if the distance between the original and slice is large.
		bool is_large = compute_dist(dist);
		
		if (is_large)
		{
			required_time_steps.push_back(pow(2,k)*step + pow(2,k));
			required_time_steps.push_back(pow(2,k)*step + pow(2,k-1));
			required_time_steps.push_back(pow(2,k)*step);
		}

		step++;
		time_step = pow(2,k)*step + pow(2,k-1);
	}
}