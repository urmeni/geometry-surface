//
//  mesh_functions.hpp
//  test1
//
//  Created by Bac Alexandra on 17/02/2022.
//

#ifndef mesh_functions_hpp
#define mesh_functions_hpp

#include <iostream>
#include "starter_defs.h"
#include <CGAL/IO/Color.h>

using namespace std ;

vector<Mesh::Vertex_index> first_vneighbors(Mesh::Vertex_index v, Mesh &m) ;

CGAL::Color color_ramp(double vmin, double vmax, double mean, double stdev, double v, CGAL::Color col1, CGAL::Color col2) ;

#endif /* mesh_functions_hpp */
