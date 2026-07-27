#pragma once
#include "pch.h"

/*
#=====================================
# Type Definition
#=====================================
*/

using Real = double;
using Integer = int;
using Vector2 = Eigen::Vector<Real, 2>;
using Vector3 = Eigen::Vector<Real, 3>;
using Vector4 = Eigen::Vector<Real, 4>;
using Matrix3x3 = Eigen::Matrix<Real, 3, 3, Eigen::RowMajor>;

using DEMMatrix = Matrix3x3;
using EBPMStiffnessMatrix = Eigen::Matrix<Real, 12, 12, Eigen::RowMajor>;
using EBPMForceDisplacementVector = Eigen::Vector<Real, 12>;

// Structual Fields
template<typename T>
using StructField = Eigen::Vector<T, -1>;

template<typename T>
using StructField2 = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

template<typename T>
using StructFieldObj = Eigen::Vector<T*, -1>;

template<typename T>
using StructFieldObj2 = Eigen::Matrix<T*, -1, -1, Eigen::RowMajor>;
