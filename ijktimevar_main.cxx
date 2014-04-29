/// \file ijktimevar.cxx
/// generate isosurface from scalar field
/// Version 0.3.0



#include <iostream>

#include "ijktimevarIO.h"
#include "ijktimevar_sampled_grid.h"
#include "ijktimevar_compare.h"
#include "ijkslice.h"
#include "ijkIO.txx"

using namespace IJK;
using namespace IJKTIMEVAR;
using namespace IJKSLICE;
using namespace std;

// local subroutines
void memory_exhaustion();

void construct_isosurface
	(const IO_INFO & io_info, const MC_DATA & mc_data,
	MCUBE_TIME & mcube_time, IO_TIME & io_time);


void construct_isosurface_sampled_grid
	(const IO_INFO & io_info, const MC_DATA & mc_data,
	MCUBE_TIME & mcube_time, MC_ISOSURFACE & mc_isosurface, IO_TIME & io_time);

void construct_interval_volume
	(const IO_INFO & io_info, const MC_DATA & mc_data,
	MCUBE_TIME & mcube_time, IO_TIME & io_time);


// **************************************************
// MAIN
// **************************************************

int main(int argc, char **argv)
{
	time_t start_time;
	time(&start_time);
	string isotable_filename;

	MCUBE_TIME mcube_time;
	IO_TIME io_time = {0.0, 0.0, 0.0};
	IO_INFO io_info;
	IJK::ERROR error;

	try {

		std::set_new_handler(memory_exhaustion);

		get_isotable_directory(io_info.isotable_directory);
		parse_command_line(argc, argv, io_info);

		MC_SCALAR_GRID full_scalar_grid;
		read_nrrd_file
			(io_info.input_filename, full_scalar_grid,  io_info, io_time);

		if (!check_input(io_info, full_scalar_grid, error)) 
		{ throw(error); };

		cout <<"full scalar grid ";
		for (int d = 0; d < full_scalar_grid.Dimension(); d++)
		{
			cout <<"axis "<< d <<" : "<< full_scalar_grid.AxisSize(d)<< endl;
		}
		cout <<endl;

		const int dimension = full_scalar_grid.Dimension();
		const int time_axis = dimension-1;
		cout <<"dim " << dimension << endl;

		// Generate 4D for the original grid.
		// set mc datastructures and flags
		MC_DATA mc_full_data;
		// Note: mc_data.SetScalarGrid must be called before set_mc_data.
		//NOTE this is on the sampled grid
		mc_full_data.SetScalarGrid
			(full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
			io_info.flag_supersample, io_info.supersample_resolution);
		set_mc_data(io_info, mc_full_data, mcube_time);
		report_num_cubes(full_scalar_grid, io_info, mc_full_data);

		// Note: All flags in mc_data should be set before calling 
		//       read isosurface lookup tables.
		// read isosurface lookup tables
		read_poly_isotable(io_info.isotable_directory, mc_full_data, io_time);
		MC_ISOSURFACE mc_full_isosurface;
		
		construct_isosurface_sampled_grid
			(io_info, mc_full_data, mcube_time, mc_full_isosurface, io_time);
		
		const int numv_per_simplex_full = mc_full_data.isotable.cube.NumVerticesPerSimplex();
		const SCALAR_TYPE isovalue = io_info.isovalue[0];
		MCUBE_INFO mcube_info(dimension);

		OUTPUT_INFO output_info;
		set_output_info(mc_full_data.isotable.cube, io_info, 0, output_info);

		VERTEX_INDEX nums_full = 
			mc_full_isosurface.simplex_vert.size()/numv_per_simplex_full;
		const VERTEX_INDEX numv_full = mc_full_isosurface.vertex_coord.size()/dimension-1;
		/*
		//store the time steps required for qulaity reconstruction.
		vector<int> required_time_steps;
		const int minK = 1;
		const int maxK = 6;
		for (int k = minK; k < maxK; k++)
		{
			MC_SCALAR_GRID sampled_grid;	
			vector <int> vec_time;
			//initialize the sampled grid.
			// returns false if num time slices is zero.
			bool flag_non_zero_timesteps  = initialize_sampled_grid 
				(k, dimension, full_scalar_grid, vec_time, sampled_grid);
	
			if (flag_non_zero_timesteps == false)
			{
				cout <<"break"<<endl;
				break;
			}

			//set copy_region_size.
			IJK::ARRAY<int> copy_region_size(dimension,0);
			determine_copy_region_size
				(dimension, full_scalar_grid, time_axis,
				copy_region_size);
			//populate the sampled grid.
			populate_grid
				(vec_time, time_axis, dimension, copy_region_size, 
				full_scalar_grid, sampled_grid);


			//output_isosurface
			//	(output_info, mc_full_data, mc_full_isosurface, mcube_info, io_time);

			//Generate 4D for the sampled grid.
			MC_DATA mc_sampled_data;
			// Note: mc_data.SetScalarGrid must be called before set_mc_data.
			//NOTE this is on the sampled grid
			mc_sampled_data.SetScalarGrid
				(sampled_grid, io_info.flag_subsample, io_info.subsample_resolution,
				io_info.flag_supersample, io_info.supersample_resolution);
			set_mc_data(io_info, mc_sampled_data, mcube_time);
			report_num_cubes(sampled_grid, io_info, mc_sampled_data);
			read_poly_isotable(io_info.isotable_directory, mc_sampled_data, io_time);
			MC_ISOSURFACE mc_sampled_isosurface;
			construct_isosurface_sampled_grid
				(io_info, mc_sampled_data, mcube_time, mc_sampled_isosurface, io_time);

			//rescale t axis
			rescale_time_axis(k, dimension, vec_time, mc_sampled_isosurface);
			//

			const int numv_per_simplex_sampled = mc_sampled_data.isotable.cube.NumVerticesPerSimplex();
			VERTEX_INDEX nums_sampled = 
				mc_sampled_isosurface.simplex_vert.size()/numv_per_simplex_sampled;
			const VERTEX_INDEX numv_sampled = mc_sampled_isosurface.vertex_coord.size()/dimension-1;
			//
			
			compare_isosurfaces(dimension, k, isovalue, numv_sampled, nums_sampled,
				mc_sampled_isosurface, numv_full, nums_full, mc_full_isosurface,
				full_scalar_grid, sampled_grid, required_time_steps, io_info);
			
		}
		
		MC_SCALAR_GRID final_grid;
		vector<int> axis_size; // for the final grid.
		for (int i = 0; i < dimension; i++)
		{
			cout <<"i "<< i <<" "<< full_scalar_grid.AxisSize(i) << endl;
			axis_size.push_back(full_scalar_grid.AxisSize(i));
		}
		//Set final grid size.
		axis_size[time_axis] = required_time_steps.size();
		final_grid.SetSize(dimension, &(axis_size[0]));
		SCALAR_TYPE zero = 0.0;
		final_grid.SetAll(zero);
		//set copy_region_size.
		IJK::ARRAY<int> copy_region_size(dimension,0);
		determine_copy_region_size
			(dimension, full_scalar_grid, time_axis,
			copy_region_size);
		//populate the final grid.
		populate_grid
			(required_time_steps, time_axis, dimension, copy_region_size, 
			full_scalar_grid, final_grid);

		//
		MC_DATA mc_final_data;
		mc_final_data.SetScalarGrid
			(final_grid, io_info.flag_subsample, io_info.subsample_resolution,
			io_info.flag_supersample, io_info.supersample_resolution);
		set_mc_data(io_info, mc_final_data, mcube_time);
		report_num_cubes(final_grid, io_info, mc_final_data);

		// Note: All flags in mc_data should be set before calling 
		//       read isosurface lookup tables.
		// read isosurface lookup tables
		read_poly_isotable(io_info.isotable_directory, mc_final_data, io_time);
		MC_ISOSURFACE mc_final_isosurface;
		
		construct_isosurface_sampled_grid
			(io_info, mc_final_data, mcube_time, mc_final_isosurface, io_time);
		
		const int numv_per_simplex_final = mc_final_data.isotable.cube.NumVerticesPerSimplex();
		set_output_info(mc_final_data.isotable.cube, io_info, 0, output_info);

		VERTEX_INDEX nums_final = 
			mc_final_isosurface.simplex_vert.size()/numv_per_simplex_final;
		const VERTEX_INDEX numv_final = mc_final_isosurface.vertex_coord.size()/dimension-1;

		cout <<"output 4D mesh"<<endl;
		output_isosurface
			(output_info, mc_final_data, mc_final_isosurface, mcube_info, io_time);
		*/
		//
		//
		cout <<"this is a test"<<endl;
		if (io_info.report_time_flag) {

			time_t end_time;
			time(&end_time);
			double total_elapsed_time = difftime(end_time, start_time);

			cout << endl;
			report_time(io_info, io_time, mcube_time, total_elapsed_time);
		};

	} 
	catch (ERROR & error) {
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		}
		else { error.Print(cerr); }
		cerr << "Exiting." << endl;
		exit(20);
	}
	catch (...) {
		cerr << "Unknown error." << endl;
		exit(50);
	};


}

//Construct isosurface of the sampled grid.
//return the surface
void construct_isosurface_sampled_grid
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



void construct_isosurface
	(const IO_INFO & io_info, const MC_DATA & mc_data,
	MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
	const int dimension = mc_data.ScalarGrid().Dimension();
	const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
	const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

	io_time.write_time = 0;
	for (unsigned int i = 0; i < io_info.isovalue.size(); i++) {

		const SCALAR_TYPE isovalue = io_info.isovalue[i];

		MC_ISOSURFACE mc_isosurface;
		MCUBE_INFO mcube_info(dimension);
		mcube_info.grid.num_cubes = num_cubes;
		SNAP_INFO snap_info(dimension);
		snap_info.grid.num_cubes = num_cubes;

		if (mc_data.Snap()) {

			if (io_info.use_list) {
				std::vector<VERTEX_INDEX> cube_list;

				float preprocessing_time;
				IJKSNAPMC::get_nonempty_snap_cubes
					(mc_data.ScalarGrid(), mc_data.isotable.cube_nep, 
					isovalue, cube_list, preprocessing_time);

				mcube_time.preprocessing += preprocessing_time;
				mcube_time.total += preprocessing_time;

				snapMC(mc_data, isovalue, cube_list, mc_isosurface, snap_info);
			}
			else {
				snapMC(mc_data, isovalue, mc_isosurface, snap_info);
			}

			mcube_time.Add(snap_info.time);
		}
		else {
			marching_cubes(mc_data, isovalue, mc_isosurface, mcube_info);
			mcube_time.Add(mcube_info.time);
		}


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

		if (mc_data.Snap()) {
			output_snapmc_isosurface
				(output_info, mc_data, mc_isosurface, snap_info, io_time);
		}
		else if (mc_data.UseNEP()) {
			output_nep_isosurface
				(output_info, mc_data, mc_isosurface, mcube_info, io_time);
		}
		else if (io_info.flag_color_alternating &&
			mc_isosurface.cube_containing_simplex.size() == nums) {
				output_isosurface_color_alternating
					(output_info, mc_data, mc_isosurface, mcube_info, io_time);
		}
		else {
			output_isosurface
				(output_info, mc_data, mc_isosurface, mcube_info, io_time);
		}
	}
}

void construct_interval_volume
	(const IO_INFO & io_info, const MC_DATA & mc_data,
	MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
	const int dimension = mc_data.ScalarGrid().Dimension();
	const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
	const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

	io_time.write_time = 0;
	for (unsigned int i = 0; i+1 < io_info.isovalue.size(); i++) {

		MC_ISOSURFACE mc_ivol;
		MCUBE_INFO mcube_info(dimension);
		mcube_info.grid.num_cubes = num_cubes;
		SNAP_INFO snap_info(dimension);
		snap_info.grid.num_cubes = num_cubes;

		MCVol(mc_data, io_info.isovalue[i], io_info.isovalue[i+1], 
			mc_ivol, mcube_info);

		mcube_time.Add(mcube_info.time);

		OUTPUT_INFO output_info;
		set_output_info(mc_data.isotable.cube, io_info, i, output_info);

		VERTEX_INDEX nums = 
			mc_ivol.simplex_vert.size()/numv_per_simplex;

		int grow_factor = 1;
		int shrink_factor = 1;
		if (io_info.flag_subsample) 
		{ grow_factor = io_info.subsample_resolution; }
		if (io_info.flag_supersample) 
		{ shrink_factor = io_info.supersample_resolution; }

		rescale_vertex_coord(grow_factor, shrink_factor, io_info.grid_spacing,
			mc_ivol.vertex_coord);

		output_isosurface
			(output_info, mc_data, mc_ivol, mcube_info, io_time);
	}

}


void memory_exhaustion()
{
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}

