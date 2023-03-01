#pragma once

#include <prism/core.h>

#include <string>


void export_to_vtu(const prism::mesh::PMesh& pmesh, const std::string& vtu_filename);