#include "utils.h"

Matrix3x3 Zero3x3()
{
	return Eigen::Matrix<Real, 3, 3, Eigen::RowMajor>::Zero();
}

// Add a math function : quaternion to rotation matrix
// References :
// https ://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
// Lines 511 - 534, https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/master/src/math_extra_liggghts_nonspherical.h
Matrix3x3 quat2RotMatrix(const Vector4& quat)
{
	// w i j k
	// 0 1 2 3
	Real w2 = quat[0] * quat[0];
	Real i2 = quat[1] * quat[1];
	Real j2 = quat[2] * quat[2];
	Real k2 = quat[3] * quat[3];

	Real twoij = 2.0 * quat[1] * quat[2];
	Real twoik = 2.0 * quat[1] * quat[3];
	Real twojk = 2.0 * quat[2] * quat[3];
	Real twoiw = 2.0 * quat[1] * quat[0];
	Real twojw = 2.0 * quat[2] * quat[0];
	Real twokw = 2.0 * quat[3] * quat[0];

	Matrix3x3 result = Zero3x3();
	result(0, 0) = w2 + i2 - j2 - k2;
	result(0, 1) = twoij - twokw;
	result(0, 2) = twojw + twoik;
	result(1, 0) = twoij + twokw;
	result(1, 1) = w2 - i2 + j2 - k2;
	result(1, 2) = twojk - twoiw;
	result(2, 0) = twoik - twojw;
	result(2, 1) = twojk + twoiw;
	result(2, 2) = w2 - i2 - j2 + k2;

	return result;
}

double sgn(const double& val)
{
	if (val < DoublePrecisionTolerance) return 0.0;
	else return val / abs(val);
}