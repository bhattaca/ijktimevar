#include <assert.h>

#include <iostream>
#include <vector>

#include "ijktable.h"
#include "ijktimevar_types.h"
#include "ijktimevar_datastruct.h"
#include "ijktimevarIO.h"
#include "ijkIO.txx"
#include "ijkslice.h"

using namespace IJK;
using namespace IJKTABLE;
using namespace IJKTIMEVAR;
using namespace IJKSLICE;
//Compare the exact isosurface with sliced ones.
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
	);