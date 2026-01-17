//
//  starter_defs.h
//  test1
//
//  Created by Bac Alexandra on 13/02/2022.
//

#ifndef starter_defs_h
#define starter_defs_h

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

#endif /* starter_defs_h */
