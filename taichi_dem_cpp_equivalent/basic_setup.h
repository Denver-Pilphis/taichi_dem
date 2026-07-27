#pragma once
#include "pch.h"

//======================================================================
// basic setup
//======================================================================

const bool SAVE_FRAMES = false;
const int window_size = 1024; // Number of pixels of the window
const int grid_n = 128;

const double DoublePrecisionTolerance = 1e-12;

// Debug-convenience parameters

// General properties
const double time_increment = 2.56e-6;
const double total_time = 1.28;
const double save_interval_time = 2.56e-3;

const double gravity_x = 0.0;
const double gravity_y = 0.0;
const double gravity_z = -9.81;

// Bond properties
const bool bond = false;

// Wall properties
const double wall_normal_x = 0.0;
const double wall_normal_y = 0.0;
const double wall_notmal_z = -1.0;
const double wall_distance = 0.03125;