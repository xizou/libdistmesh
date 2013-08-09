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
# along with libDistMesh.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2013 Patrik Gebhardt
# Contact: patrik.gebhardt@rub.de
# --------------------------------------------------------------------

import sconshelper

# create environment
env = Environment()

# use clang++
env.Replace(CXX='clang++')

# create librarty
sconshelper.Library(name='distmesh', env=env, arguments=ARGUMENTS,
    CXXFLAGS=[
        '-std=c++11',
        '-O4',
        ],
    LINKFLAGS=[
        '-O4',
        ],
    CPPPATH=[
        '/usr/include/eigen3/',
        '/usr/local/include/eigen3/',
        ],
    LIBS=[
        'qhull',
        ],
    )
