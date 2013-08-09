// --------------------------------------------------------------------
// This file is part of libDistMesh.
//
// libDistMesh is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// libDistMesh is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libDistMesh.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2013 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#ifndef LIBDISTMESH_INCLUDE_DISTMESH_H
#define LIBDISTMESH_INCLUDE_DISTMESH_H

#include "common.h"
#include "dtype.h"
#include "settings.h"
#include "triangulation.h"
#include "utils.h"
#include "distance_functions.h"
#include "edge_length_functions.h"

// namespace distmesh
namespace distmesh {
    // apply the distmesh algorithm
    std::tuple<std::shared_ptr<dtype::array<dtype::real>>,
        std::shared_ptr<dtype::array<dtype::index>>> distmesh(
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> distance_function,
        std::function<dtype::array<dtype::real>(dtype::array<dtype::real>&)> edge_length_function,
        dtype::real initial_edge_length, dtype::array<dtype::real> bounding_box,
        dtype::array<dtype::real> fixed_points=dtype::array<dtype::real>());

    // determine boundary edges of given triangulation
    std::shared_ptr<dtype::array<dtype::index>> boundedges(
        std::shared_ptr<dtype::array<dtype::index>> triangulation);
}

#endif
