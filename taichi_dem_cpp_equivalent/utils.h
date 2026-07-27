#pragma once
#include "pch.h"

/*
#=====================================
# Utils
#=====================================
*/

Matrix3x3 Zero3x3();

Matrix3x3 quat2RotMatrix(const Vector4& quat);

double sgn(const double& val);