# --------------------------------------------------------------------
# This file is part of libDistMesh.
#
# libDistMesh is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# libDistMesh is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libDistMesh. If not, see <http://www.gnu.org/licenses/>.
#
# --------------------------------------------------------------------

from pylab import *


def plot_edges(filename, points, color='b'):
    try:
        edgeList = loadtxt(filename, dtype=int)
        for i in range(edgeList.shape[0]):
            plot(points[edgeList[i, :], 0], points[edgeList[i, :], 1], color)
    except:
        return


def main():
    # load mesh
    points = loadtxt('points.txt')
    triangulation = loadtxt('triangulation.txt', dtype=int)

    # plot mesh
    gca().set_aspect('equal')
    triplot(points[:, 0], points[:, 1], triangulation)

    plot_edges('edgeBoundary.txt', points)
    plot_edges('edgeInterior.txt', points, 'r')

    show()


if __name__ == '__main__':
    main()
