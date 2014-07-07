#include "ijktimevar_sampled_grid.h"
#include "ijkmcube.h"
#include "ijkslice.h"
#include "ijkIO.txx"
#include "ijktimevarIO.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;
using namespace IJKSLICE;

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

// add a single slice. 
// called by compute_local_grid.
void add_slice(
	const int a, //time axis temp_coord for scalar grid
	const int b, //time axis temp_coord for sampled grid
	const vector<int> &copy_size,
	const MC_SCALAR_GRID_BASE & scalar_grid,
	MC_SCALAR_GRID & sampled_grid)
{
	cout <<"time axis scalar grid "<< a <<", sampled grid "<< b<< endl;
	const int dimension = scalar_grid.Dimension();
	IJK::ARRAY<COORD_TYPE> temp_coord(dimension,0);
	temp_coord[dimension-1] = a;
	cout <<"temp_coord "<< temp_coord[0]<<" "<< temp_coord[1]<<" "<< temp_coord[2]
	<<" "<< temp_coord[3]<< endl;


	VERTEX_INDEX vstart_scalar_grid = scalar_grid.ComputeVertexIndex(temp_coord.Ptr());
	cout <<"vstart scalar grid "<< vstart_scalar_grid << endl;
	temp_coord[dimension-1] = b;
	cout <<"temp_coord "<< temp_coord[0]<<" "<< temp_coord[1]<<" "<< temp_coord[2]
	<<" "<< temp_coord[3]<< endl;
	VERTEX_INDEX vstart_sampled_grid = sampled_grid.ComputeVertexIndex(temp_coord.Ptr());
	cout <<"vstart sample grid "<< vstart_sampled_grid << endl;
	cout <<"copy region "<< endl;
	cout <<"copy region size"<< endl; 

	for (int d = 0; d < dimension; d++)
	{
		cout <<d <<" "<< copy_size[d]<< endl;
	}

	sampled_grid.CopyRegion(scalar_grid, vstart_scalar_grid,
		&(copy_size[0]), vstart_sampled_grid);

}
void compute_local_grid (
	const int t1,
	const int t2, 
	const vector<int> copy_size,
	const MC_SCALAR_GRID_BASE & scalar_grid,
	MC_SCALAR_GRID & sampled_grid)
{

	int b = 0;
	add_slice(t1, b, copy_size, scalar_grid, sampled_grid);

	b = 1;
	add_slice(t2, b, copy_size, scalar_grid, sampled_grid);

}


void time_scaling (
	const int dimension,
	const int numVertices, 
	const int ta,
	vector<COORD_TYPE> & vertex_coord, 
	const int k
	)
{
	cout <<"numvertices "<< numVertices << endl;
	cout <<"range "<<pow(2,k)<<endl;

	for (int i = 0; i < numVertices; i++)
	{
		COORD_TYPE t = vertex_coord[i*dimension + dimension-1];
		int t_low = int(vertex_coord[i*dimension + dimension-1]);
		//cout <<"t "<< t ; 
		//cout <<",t_low "<< t_low;
		vertex_coord[i*dimension + dimension-1] = ta + (t - t_low)* pow(2,k); 
		//cout <<", new coord "<<vertex_coord[i*dimension + dimension-1]<<endl; 
	}
}


//Construct isosurface of the sampled grid.
//return the surface
void construct_isosurface_sampled_gridA
	(const IO_INFO & io_info, const MC_DATA & mc_data,
	MCUBE_TIME & mcube_time, MC_ISOSURFACE & mc_isosurface,
	IO_TIME & io_time)
{
	const int dimension = mc_data.ScalarGrid().Dimension();
	const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
	const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

	io_time.write_time = 0;
	for (unsigned int i = 0; i < io_info.isovalue.size(); i++) {

		const SCALAR_TYPE isovalue = io_info.isovalue[i];

		//MC_ISOSURFACE mc_isosurface;
		MCUBE_INFO mcube_info(dimension);
		mcube_info.grid.num_cubes = num_cubes;
		SNAP_INFO snap_info(dimension);
		snap_info.grid.num_cubes = num_cubes;

		//compute marching cubes,
		//generate the 4D surface.
		marching_cubes(mc_data, isovalue, mc_isosurface, mcube_info);
		mcube_time.Add(mcube_info.time);

		OUTPUT_INFO output_info;
		set_output_info(mc_data.isotable.cube, io_info, i, output_info);

		VERTEX_INDEX nums = 
			mc_isosurface.simplex_vert.size()/numv_per_simplex;

		int grow_factor = 1;
		int shrink_factor = 1;
		if (io_info.flag_subsample) 
		{ grow_factor = io_info.subsample_resolution; }
		if (io_info.flag_supersample) 
		{ shrink_factor = io_info.supersample_resolution; }

		rescale_vertex_coord(grow_factor, shrink_factor, io_info.grid_spacing,
			mc_isosurface.vertex_coord);

	}
}

void evaluate_multi_slice (
	const int ta, 
	const int tc,
	const SCALAR_TYPE & isovalue,
	POLY_ISOTABLE & isotable,
	const string & isotable_directory,
	const MC_SCALAR_GRID_BASE & scalar_grid, 
	vector<int> rqd_tsteps,
	IO_INFO &io_info)
{
	//debug
	cout <<"ta "<< ta <<", tc "<< tc << endl;
	const int dimension = scalar_grid.Dimension();
	vector<int> axis_size; // for the sampled grid.
	vector<int> copy_size; // for the sampled grid.
	for (int i = 0; i < dimension; i++)
	{
		axis_size.push_back(scalar_grid.AxisSize(i));
		copy_size.push_back(scalar_grid.AxisSize(i));
	}
	//Set sampled grid size.
	axis_size[dimension-1] = 2;
	copy_size[dimension-1] = 1;
	int k = 1; 
	const int maxK = 1;
	while(k <= maxK){
		cout <<" Loop k "<< k << endl;
		const int twok = pow(2, k);
		int i = ta/twok;

		while(pow(2,k)*i >= ta && (pow(2,k)*i + pow(2,k)) <= tc){
			//debug 
			//cout <<"i "<< i << endl;
			int t1 = pow(2,k)*i;
			int t3 = pow(2,k)*i + pow(2,k);
			int t2 = pow(2,k)*i + pow(2,k-1);
			//debug 
			cout <<"t1 "<< t1 <<", t2 "<< t2 <<", t3 "<<  t3 << endl;

			MC_SCALAR_GRID local_grid;
			local_grid.SetSize(dimension, &(axis_size[0]));
			SCALAR_TYPE zero = 0.0;
			local_grid.SetAll(zero);
			//debug 
			//compute_local_grid(t1, t3, copy_size, scalar_grid, local_grid);
			//debug
			cout <<"local grid size "<< endl;
			for (int d = 0; d < dimension; d++)
			{
				cout <<d <<" "<< local_grid.AxisSize(d) << endl;
			}
			vector<VERTEX_INDEX> simplex_vert;
			vector<COORD_TYPE> vertex_coord;
			//DEBUG ::  scalar grid instead of local grid/
			cout <<"marching cubes"<< endl;
			/*
			marching_cubes(scalar_grid, isotable.cube, isovalue, simplex_vert, vertex_coord);
			*/
			//
			MCUBE_TIME mcube_time;
			IO_TIME io_time = {0.0, 0.0, 0.0};
			
			MC_DATA mc_data;
			// Note: mc_data.SetScalarGrid must be called before set_mc_data.
			//NOTE this is on the sampled grid
			mc_data.SetScalarGrid
				(scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
				io_info.flag_supersample, io_info.supersample_resolution);
			set_mc_data(io_info, mc_data, mcube_time);
			report_num_cubes(scalar_grid, io_info, mc_data);

			// Note: All flags in mc_data should be set before calling 
			//       read isosurface lookup tables.
			// read isosurface lookup tables
			read_poly_isotable(io_info.isotable_directory, mc_data, io_time);
			MC_ISOSURFACE mc_full_isosurface;

			construct_isosurface_sampled_gridA
				(io_info, mc_data, mcube_time, mc_full_isosurface, io_time);

			const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
			const SCALAR_TYPE isovalue = io_info.isovalue[0];
			MCUBE_INFO mcube_info(dimension);

			OUTPUT_INFO output_info;
			set_output_info(mc_data.isotable.cube, io_info, 0, output_info);

			VERTEX_INDEX nums = 
				mc_full_isosurface.simplex_vert.size()/numv_per_simplex;
			const VERTEX_INDEX numv = mc_full_isosurface.vertex_coord.size()/dimension-1;

			//
			
			//DEBUG PROBLEM
			//const COORD_TYPE  slice_coord = float(t2);
			const COORD_TYPE slice_coord = float (1.4);
			const int slice_axis = dimension-1;
			int * intersected_edge_vert = NULL;
			int num_intersected_edges = 0;
			int * slice_simplex_vert = NULL;
			int num_slice_simplices = 0;
			
			cout <<" numv per simplex "<< numv_per_simplex << endl; 
	/*		VERTEX_INDEX nums = 
				simplex_vert.size()/numv_per_simplex;
			const VERTEX_INDEX numv = vertex_coord.size()/dimension-1;*/

			cout <<"nums "<< nums <<" numv " << numv << endl;
			//time_scaling(dimension, numv, t1, vertex_coord, k);
			cout <<"done time scaling " << endl;
			cout <<"slice axis "<< slice_axis <<" coord "<< slice_coord << endl;
			slice_simplices( dimension, &(mc_full_isosurface.vertex_coord[0]), numv, dimension-1, 
				&(mc_full_isosurface.simplex_vert[0]), nums, slice_axis, slice_coord, isotable_directory.c_str(),
				intersected_edge_vert, num_intersected_edges, slice_simplex_vert, 
				num_slice_simplices);

			cout <<" done slicing "<< endl;
			//print data
			cout <<"print data "<< endl;
			int num_slice_vert_sampled = num_intersected_edges;
			COORD_TYPE * coord_sampled = new COORD_TYPE[(dimension-1)*num_slice_vert_sampled];
			int * coord0 = new int[dimension-1];
			int * coord1 = new int[dimension-1];

			for (int iv = 0; iv < num_slice_vert_sampled; iv++) {
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
						coord_sampled[iv*(dimension-1) + d2] = c;
						d2++;
					};
				};
			}

			cout <<"num slice vert sampled " << num_slice_vert_sampled 
				<<" num_slice_simplices "<< num_slice_simplices << endl;
			ostream * out = NULL;
			cout <<"output mesh"<< endl;
			out = new ofstream("sampled.off", ios::out);
			ijkoutOFF(*out, dimension-1, coord_sampled, num_slice_vert_sampled, 
				slice_simplex_vert, num_slice_simplices);

			cout <<"output 4D mesh"<<endl;
			output_isosurface
				(output_info, mc_data, mc_full_isosurface, mcube_info, io_time);

			// 
			i++;
		}
		k++;
	}
}

