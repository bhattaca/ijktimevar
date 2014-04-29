#include "ijktimevar_sampled_grid.h"
#include <iostream>
#include <vector>

using namespace std;
void compute_time_steps(
	const int twoK,
	const int time_steps,
	vector<int> &vec_time
	)
{
	int i=0;
	//DEBUG
	cout <<"compute time steps ";
	while (twoK*i < time_steps)
	{
		cout << twoK*i <<" ";
		vec_time.push_back(twoK*i);
		i++;
	}
	cout <<"\n";
}

//Generate sampled grid.
bool  initialize_sampled_grid  (
	const int k , 
	const int dimension,
	const MC_SCALAR_GRID_BASE & full_scalar_grid, 
	//returns
	vector <int> &vec_time,
	MC_SCALAR_GRID & sampled_grid)
{
	vector<int> axis_size; // for the sampled grid.
	for (int i = 0; i < dimension; i++)
	{
		cout <<"i "<< i <<" "<< full_scalar_grid.AxisSize(i) << endl;
		axis_size.push_back(full_scalar_grid.AxisSize(i));
	}
	const int twoK = pow(2,k);
	const int time_axis = dimension-1;
	cout <<"time axis " << dimension-1 << endl;
	const int time_steps = full_scalar_grid.AxisSize(time_axis);

	//Compute time steps.
	compute_time_steps(twoK, time_steps, vec_time);
	cout <<"number of time steps "<< vec_time.size() << endl;

	if (vec_time.size() <= 1)
		return false;

	//Set sampled grid size.
	axis_size[time_axis] = vec_time.size();
	sampled_grid.SetSize(dimension, &(axis_size[0]));
	SCALAR_TYPE zero = 0.0;
	sampled_grid.SetAll(zero);
	cout <<"sampled grid "<< endl;
	for (int d = 0; d < sampled_grid.Dimension(); d++)
	{
		cout <<" "<< sampled_grid.AxisSize(d);
	}
	cout <<"\n";

	return true;
}

void determine_copy_region_size
	(const int dimension,
	const MC_SCALAR_GRID_BASE & full_scalar_grid,
	const int time_axis,
	IJK::ARRAY<int> & copy_region_size)
{
	//set copy_region_size
	for (int d = 0; d < dimension-1; d++)
	{
		copy_region_size[d] = full_scalar_grid.AxisSize(d);
	}
	copy_region_size[time_axis] = 1;
}

void populate_grid
	(
	const vector<int> &vec_time,
	const int time_axis,
	const int dimension,
	IJK::ARRAY<int> & copy_region_size,
	const MC_SCALAR_GRID_BASE & full_scalar_grid, 
	MC_SCALAR_GRID & sampled_grid
	)
{
	for (int t = 0; t < vec_time.size(); t++)
	{
		//cout <<"t "<< t <<" vec_time "<< vec_time[t]<< endl;
		IJK::ARRAY<COORD_TYPE> temp_coord(dimension,0);
		temp_coord[time_axis] = vec_time[t];
		//DEBUG 
		//cout <<"temp_coord "<< temp_coord[0]
		//<<" "<< temp_coord[1]<<" "<< temp_coord[2]<<" "
		//<<temp_coord[3]<<endl;

		VERTEX_INDEX vstart_scalar_grid = full_scalar_grid.ComputeVertexIndex(temp_coord.Ptr());
		//cout <<" vertex index "<< vstart_scalar_grid <<endl;
		temp_coord[time_axis] = t;
		VERTEX_INDEX vstart_sampled_grid = sampled_grid.ComputeVertexIndex(temp_coord.Ptr());
		//DEBUG 
		//cout <<"temp_coord "<< temp_coord[0]
		//<<" "<< temp_coord[1]<<" "<< temp_coord[2]<<" "
		//<<temp_coord[3]<<endl;
		//cout <<" vertex index "<< vstart_sampled_grid <<endl;

		sampled_grid.CopyRegion(full_scalar_grid, vstart_scalar_grid,
			&(copy_region_size[0]), vstart_sampled_grid);
	}
}

void rescale_time_axis
	(
	const int k, //level
	const int dimension,
	const vector<int> vec_time,
	MC_ISOSURFACE & mc_isosurface
	)
{
	const int spacing = pow(2,k);
	const int num_vertices = mc_isosurface.vertex_coord.size() / dimension;
	cout <<"numvertices "<< num_vertices << endl;
	for (int i = 0; i < num_vertices; i++)
	{
		float w = mc_isosurface.vertex_coord[dimension*i + (dimension-1)];
		int w_low = int(w);
		int t = vec_time[w_low];
		float scaled_z_coord = t + (w-w_low)*spacing;
		/*cout <<"w "<< w <<" w_low "<< w_low << " vec "<< vec_time[w_low]
		<< " spacing "<< spacing <<" scaled value "<< scaled_z_coord <<  endl;*/
		mc_isosurface.vertex_coord[dimension*i + (dimension-1)] = scaled_z_coord;
	}
}



void copy2_sampled_grid ( 
	const MC_SCALAR_GRID_BASE & scalar_grid, 
	const VERTEX_INDEX vstart_scalar_grid,
	const VERTEX_INDEX vstart_sampled_grid,
	IJK::ARRAY<int> sampled_axis_size,
	//returns
	MC_SCALAR_GRID_BASE & sampled_grid
	)
{

}
