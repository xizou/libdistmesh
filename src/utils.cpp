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
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#include <random>
#include <chrono>
#include <set>
#include <array>

#include "distmesh/distmesh.h"
#include "distmesh/settings.h"
#include "distmesh/utils.h"

// calculate factorial recursively
unsigned distmesh::utils::factorial(unsigned const n) {
    if (n <= 1) {
        return 1;
    }
    else {
        return n * factorial(n - 1);
    }
}

// create initial points distribution
Eigen::ArrayXXd distmesh::utils::createInitialPoints(
    Functional const& distanceFunction, double const baseEdgeLength,
    Functional const& edgeLengthFunction, Eigen::Ref<Eigen::ArrayXXd const> const boundingBox,
    Eigen::Ref<Eigen::ArrayXXd const> const fixedPoints) {
    // extract dimension of mesh
    int const dimensions = boundingBox.rows();

    // calculate max number of points per dimension and
    // max total point coun and create initial array
    Eigen::ArrayXi maxPointsPerDimension(dimensions);
    int maxPointCount = 1;
    for (int dim = 0; dim < dimensions; ++dim) {
        maxPointsPerDimension(dim) = ceil((boundingBox(dim, 1) - boundingBox(dim, 0)) /
            (baseEdgeLength * (dim == 0 ? 1.0 : sqrt(3.0) / 2.0)));
        maxPointCount *= maxPointsPerDimension(dim);
    }

    // fill point list with evenly distributed points
    Eigen::ArrayXXd points(maxPointCount, dimensions);
    int sameValueCount = 1;
    for (int dim = 0; dim < dimensions; ++dim) {
        for (int point = 0; point < maxPointCount; ++point) {
            points(point, dim) = boundingBox(dim, 0) +
                baseEdgeLength * (dim == 0 ? 1.0 : sqrt(3.0) / 2.0) *
                ((point / sameValueCount) % maxPointsPerDimension(dim));
            if (dim > 0) {
                points(point, dim - 1) += point / sameValueCount % 2 != 0 ? baseEdgeLength / 2.0 : 0.0;
            }
        }
        sameValueCount *= maxPointsPerDimension(dim);
    }

    // reject points outside of region defined by distance function
    points = selectMaskedArrayElements<double>(points,
        distanceFunction(points) < settings::generalPrecision * baseEdgeLength);

    // clear dublicate points
    Eigen::Array<bool, Eigen::Dynamic, 1> isDublicate =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(points.rows(), true);
    for (int i = 0; i < fixedPoints.rows(); ++i)
    for (int j = 0; j < points.rows(); ++j) {
        isDublicate(j) &= !(fixedPoints.row(i) == points.row(j)).all();
    }
    points = selectMaskedArrayElements<double>(points, isDublicate);

    // add fixed points to final list first
    Eigen::ArrayXXd finalPoints(points.rows() + fixedPoints.rows(), dimensions);
    finalPoints.block(0, 0, fixedPoints.rows(), 2) = fixedPoints;

    // calculate propability to keep points
    Eigen::ArrayXd propability = 1.0 / edgeLengthFunction(points).pow(dimensions);
    propability /= propability.maxCoeff();

    // initialize random number generator
    std::default_random_engine randomGenerator(
        std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> randomDistribution(0.0, 1.0);

    // reject points with wrong propability
    int finalPointCount = 0;
    for (int point = 0; point < propability.rows(); ++point) {
        if (randomDistribution(randomGenerator) < propability(point)) {
            finalPoints.row(finalPointCount + fixedPoints.rows()) = points.row(point);
            finalPointCount++;
        }
    }
    finalPoints.conservativeResize(finalPointCount + fixedPoints.rows(), dimensions);

    return finalPoints;
}

// create array with all unique combinations n over k
Eigen::ArrayXXi distmesh::utils::nOverK(unsigned const n, unsigned const k) {
    // fill an array with all unique combinations n over k,
    // starting with 0, 1, 2, ...
    Eigen::ArrayXXi combinations = Eigen::ArrayXXi::Zero(factorial(n) /
            (factorial(k) * factorial(n - k)), k);
    combinations.row(0).setLinSpaced(k, 0, k - 1);

    for (int combination = 1; combination < combinations.rows(); ++combination) {
        combinations.row(combination) = combinations.row(combination - 1);
        for (int col = k - 1; col >= 0; --col) {
            combinations.block(combination, col, 1, k - col).row(0)
                .setLinSpaced(k - col, combinations(combination, col) + 1,
                    combinations(combination, col) + k - col);
            if (combinations(combination, k - 1) < (int)n) {
                break;
            }
        }
    }

    return combinations;
}

// find unique bars
Eigen::ArrayXXi distmesh::utils::findUniqueBars(Eigen::Ref<Eigen::ArrayXXi const> const triangulation) {
    // find all unique combinations
    auto combinations = nOverK(triangulation.cols(), 2);

    // find unique bars for all combinations
    std::set<std::array<int, 2>> bar_set;
    std::array<int, 2> bar = {{0, 0}};
    for (int combination = 0; combination < combinations.rows(); ++combination)
    for (int triangle = 0; triangle < triangulation.rows(); ++triangle) {
        bar[0] = triangulation(triangle, combinations(combination, 0));
        bar[1] = triangulation(triangle, combinations(combination, 1));

        bar_set.insert(bar);
    }

    // copy set to eigen array
    Eigen::ArrayXXi bar_indices(bar_set.size(), 2);
    Eigen::ArrayXXi::Index bar_index = 0;
    for (const auto& bar : bar_set) {
        bar_indices(bar_index, 0) = bar[0];
        bar_indices(bar_index, 1) = bar[1];
        bar_index++;
    }

    return bar_indices;
}

// project points outside of boundary back to it
void distmesh::utils::projectPointsToFunction(
    Functional const& distanceFunction, double const baseEdgeLength,
    Eigen::Ref<Eigen::ArrayXXd> points) {
    // evaluate distance function at points
    Eigen::ArrayXd distance = distanceFunction(points);

    // check for points outside of boundary
    Eigen::Array<bool, Eigen::Dynamic, 1> outside = distance > 0.0;
    if (outside.any()) {
        // calculate gradient
        Eigen::ArrayXXd gradient(points.rows(), points.cols());
        Eigen::ArrayXXd deltaX(points.rows(), points.cols());
        Eigen::ArrayXXd h;
        deltaX.fill(0.0);
        for (int dim = 0; dim < points.cols(); ++dim) {
            deltaX.col(dim).fill(std::sqrt(std::numeric_limits<double>::epsilon()) * baseEdgeLength);
            h = points + deltaX;
            gradient.col(dim) = (distanceFunction(h) - distance) /
                (std::sqrt(std::numeric_limits<double>::epsilon()) * baseEdgeLength);
            deltaX.col(dim).fill(0.0);
        }

        for (int dim = 0; dim < points.cols(); ++dim) {
            points.col(dim) -= outside.select(
                gradient.col(dim) * distance / gradient.square().rowwise().sum(),
                0.0);
        }
    }
}

// check whether points lies inside or outside of polygon
Eigen::ArrayXXd distmesh::utils::pointsInsidePoly(
    Eigen::Ref<Eigen::ArrayXXd const> const points,
    Eigen::Ref<Eigen::ArrayXXd const> const polygon) {
    Eigen::ArrayXXd inside(points.rows(), 1);
    inside.fill(0.0);

    for (int i = 0, j = polygon.rows() - 1;
        i < polygon.rows(); j = i++) {
        inside = (((points.col(1) < polygon(i, 1)) != (points.col(1) < polygon(j, 1))) &&
            (points.col(0) < (polygon(j, 0) - polygon(i, 0)) * (points.col(1) - polygon(i, 1)) /
            (polygon(j, 1) - polygon(i, 1)) + polygon(i, 0))).select(1.0 - inside, inside);
    }

    return inside;
}
