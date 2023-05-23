// --------------------------------------------------------------------
// This file is an example for using libDistMesh.
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
// along with libDistMesh. If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2023 Xi Zou
// Contact: xi.zou@swansea.ac.uk
// --------------------------------------------------------------------

#include <iostream>
#include <distmesh/distmesh.h>
#include <eigen3/Eigen/Dense>
#include "helper.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    distmesh::helper::HighPrecisionTime time;

    // input parameters
    double polygonEdgeLength = 1.0;
    double angleToYAxis = 75.0 * M_PI / 180.0;
    double edgeSize = 0.05;

    // corner points of polygon
    double YMax = polygonEdgeLength * cos(angleToYAxis);
    double YMin = -1.0 * YMax;
    Eigen::ArrayXXd polygon(4, 2);
    polygon <<
        polygonEdgeLength * (-1.0), 0.0,
        polygonEdgeLength * (-1.0 + sin(angleToYAxis)), YMin,
        polygonEdgeLength * sin(angleToYAxis), YMin,
        0.0, 0.0;

    // edge size function
    // auto sizeFunction =
    //     (edgeSize/2.0 + 0.2 * distmesh::distanceFunction::polygon(polygon).abs())
    //     .min(edgeSize);

    auto sizeFunction = (1.0);

    // create mesh
    Eigen::ArrayXXd points;
    Eigen::ArrayXXi elements;

    std::tie(points, elements) = distmesh::distmesh(
        distmesh::distanceFunction::polygon(polygon),
        edgeSize,
        sizeFunction,
        distmesh::utils::boundingBox(2), polygon);

    // symmetrise the mesh w.r.t X axis
    Eigen::ArrayXXd newPoints = Eigen::ArrayXXd::Zero(2 * points.rows(), points.cols());
    Eigen::ArrayXXi newElements = Eigen::ArrayXXi::Zero(2 * elements.rows(), elements.cols());

    int nNodeOffset = points.rows();
    int newNodeCount = 0;

    Eigen::ArrayXi nodeMap = Eigen::ArrayXi::Zero(nNodeOffset);

    newPoints.topRows(nNodeOffset) = points;

    for (int i = 0; i < points.rows(); i++) {
        if (abs(points(i, 1)) < 1e-10) {
            // if the node is on X axis, do not create a new node
            nodeMap(i) = i;
            newPoints.row(nodeMap(i)) << points(i, 0), 0.0 * points(i, 1);
        } else {
            // create a new node
            nodeMap(i) = newNodeCount + nNodeOffset;
            newNodeCount++;
            newPoints.row(nodeMap(i)) << points(i, 0), -1.0 * points(i, 1);
        }
    }

    for (int i = 0; i < elements.rows(); i++) {
        // check and flip normals
        Eigen::Matrix3d mat;
        mat.col(2) << 1.0, 1.0, 1.0;
        mat.row(0) << points.row(elements(i, 0));
        mat.row(1) << points.row(elements(i, 1));
        mat.row(2) << points.row(elements(i, 2));

        if (mat.determinant() < 0.0) elements.row(i).reverseInPlace();

        for (int j = 0; j < 3; j++) {
            int nodeID = elements(i, j);
            int newNodeID = nodeMap(nodeID);
            newElements(i+elements.rows(), 2-j) = newNodeID;
        }
    }

    newElements.topRows(elements.rows()) = elements;

    points = newPoints.topRows(newNodeCount + nNodeOffset);
    elements = newElements;

    // build a matrix to count edge occurences
    Eigen::ArrayXXi edgeMat = Eigen::ArrayXXi::Zero(points.rows(), points.rows());
    const Eigen::Vector3i endIndex(1, 2, 0);
    for (int i = 0; i < elements.rows(); i++) {
        for (int j = 0; j < 3; j++) {

            int p1 = elements(i, j);
            int p2 = elements(i, endIndex(j));

            if (p1 > p2) {
                p1 = p2;
                p2 = elements(i, j);
            }

            edgeMat(p1, p2)++;
        }
    }

    // build edge lists
    Eigen::ArrayXXi edgeListBoundary(3 * elements.rows(), 2);
    Eigen::ArrayXXi edgeListInterior(3 * elements.rows(), 2);
    int edgeCountBoundary = 0;
    int edgeCountInterior = 0;

    for (int i = 0; i < points.rows(); i++) {
        for (int j = i+1; j < points.rows(); j++) {
            if (edgeMat(i, j) == 1) {
                // edges occurring once must be boundary edges
                edgeListBoundary.row(edgeCountBoundary) << i, j;
                edgeCountBoundary++;
            } else if (edgeMat(i, j) == 2) {
                // edges occurring twice must be interior edges
                edgeListInterior.row(edgeCountInterior) << i, j;
                edgeCountInterior++;
            }
        }
    }

    // collect the edges to output
    Eigen::ArrayXXi edgeBoundary = edgeListBoundary.topRows(edgeCountBoundary);
    Eigen::ArrayXXi edgeInterior = edgeListInterior.topRows(edgeCountInterior);

    // print mesh properties and elapsed time
    std::cout << "Created mesh with " << points.rows() << " points and " << elements.rows() <<
        " elements in " << time.elapsed() * 1e3 << " ms." << std::endl;

    // p = points.data();
    // e = elements.data();

    // save mesh to file
    distmesh::helper::savetxt<double>(points, "points.txt");
    distmesh::helper::savetxt<int>(elements, "triangulation.txt");
    distmesh::helper::savetxt<int>(edgeBoundary, "edgeBoundary.txt");
    distmesh::helper::savetxt<int>(edgeInterior, "edgeInterior.txt");

    // plot mesh using python
    return system("python plot_mesh.py");
}
