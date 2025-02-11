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
// along with libDistMesh. If not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2015 Patrik Gebhardt
// Contact: patrik.gebhardt@rub.de
// --------------------------------------------------------------------

#ifndef _f2f87a1d_2511_4053_996f_1156d50b31f1
#define _f2f87a1d_2511_4053_996f_1156d50b31f1

// macro for easier creation of distmesh lambda functions
#define DISTMESH_FUNCTIONAL(function_body) \
    (distmesh::Functional([=](Eigen::Ref<Eigen::ArrayXXd const> const points) -> Eigen::ArrayXd \
        function_body))

namespace distmesh {
    // base class of all function expression for allowing easy function arithmetic
    class Functional {
    public:
        // function type of Functional callable
        typedef std::function<Eigen::ArrayXd(Eigen::Ref<Eigen::ArrayXXd const> const)> function_t;

        // create class from function type
        Functional(function_t const& func) : function_(func) {}
        Functional(double const constant) :
            Functional(DISTMESH_FUNCTIONAL({
                return Eigen::ArrayXd::Constant(points.rows(),constant);
            })) {}

        // copy constructor
        Functional(Functional const& rhs) : function_(rhs.function()) {}
        Functional(Functional&& rhs) : function_(std::move(rhs.function())) {}

        // assignment operator
        Functional& operator=(Functional const& rhs);
        Functional& operator=(Functional&& rhs);

        // evaluate function by call
        Eigen::ArrayXd operator() (Eigen::Ref<Eigen::ArrayXXd const> const points) const;

        // basic arithmetic operations
        Functional operator+() const { return *this; }
        Functional operator-() const;
        Functional& operator+=(Functional const& rhs);
        Functional& operator+=(double const rhs);
        Functional& operator-=(Functional const& rhs);
        Functional& operator-=(double const rhs);
        Functional& operator*=(Functional const& rhs);
        Functional& operator*=(double const rhs);
        Functional& operator/=(Functional const& rhs);
        Functional& operator/=(double const rhs);
        friend Functional operator+(Functional const& lhs, Functional const& rhs);
        friend Functional operator+(Functional const& lhs, double const rhs);
        friend Functional operator+(double const lhs, Functional const& rhs);
        friend Functional operator-(Functional const& lhs, Functional const& rhs);
        friend Functional operator-(Functional const& lhs, double const rhs);
        friend Functional operator-(double const lhs, Functional const& rhs);
        friend Functional operator*(Functional const& lhs, Functional const& rhs);
        friend Functional operator*(Functional const& lhs, double const rhs);
        friend Functional operator*(double const lhs, Functional const& rhs);
        friend Functional operator/(Functional const& lhs, Functional const& rhs);
        friend Functional operator/(Functional const& lhs, double const rhs);
        friend Functional operator/(double const lhs, Functional const& rhs);

        // mathematical methods
        Functional min(Functional const& rhs) const;
        Functional max(Functional const& rhs) const;
        Functional abs() const;

        // geometric transform
        Functional shift(Eigen::Ref<Eigen::ArrayXd const> const offset) const;
        Functional rotate2D(double const angle) const;

        // accessors
        function_t& function() { return this->function_; }
        function_t const& function() const { return this->function_; }

    private:
        // stores std function
        function_t function_;
    };

    Functional operator+(Functional const& lhs, Functional const& rhs);
    Functional operator+(Functional const& lhs, double const rhs);
    Functional operator+(double const lhs, Functional const& rhs);
    Functional operator-(Functional const& lhs, Functional const& rhs);
    Functional operator-(Functional const& lhs, double const rhs);
    Functional operator-(double const lhs, Functional const& rhs);
    Functional operator*(Functional const& lhs, Functional const& rhs);
    Functional operator*(Functional const& lhs, double const rhs);
    Functional operator*(double const lhs, Functional const& rhs);
    Functional operator/(Functional const& lhs, Functional const& rhs);
    Functional operator/(Functional const& lhs, double const rhs);
    Functional operator/(double const lhs, Functional const& rhs);

}

#endif
