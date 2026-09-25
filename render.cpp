#include "render.hpp"

float pix_pos_x(double x) { return float(400 + x * pix::scale); }
float pix_pos_y(double y) { return float(400 - y * pix::scale); }