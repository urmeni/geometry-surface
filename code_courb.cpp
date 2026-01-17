//
//  test1.cpp
//  
//
//  Created by Bac Alexandra on 11/02/2022.
//

#include <iostream>
#include <limits>
#include <random>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Kernel/global_functions.h>
#include "starter_defs.h"
#include "mesh_functions.hpp"
#include "distrib.hpp"
#include "courbures.h"

using namespace std ;

int main(int argc, char** argv)
{
	Mesh m;
	string fname ;
	if (argc != 2)
	{
		cout << "Wrong number of arguments: code_courb filename" << endl ;
		return EXIT_FAILURE;
	}
	else
	{
		fname = argv[1] ;
	}

	if(!CGAL::IO::read_polygon_mesh(fname, m))
	{
		std::cerr << "Invalid input file." << std::endl;
		return EXIT_FAILURE;
	}
	
    // Calcul des tenseurs de courbure en chaque point (par ajustement de surface quadratique)
    cout << "----- Starting curvature computation" << endl ;
    MyMesh mym(m) ;
    MyMesh_Norm mymn(mym) ;
    MyMesh_Courb mymc(mymn) ;
    MyMesh_Quad mymq(mymn) ;
    
    // Colors from curvature
    mymq.fit_quads() ;
    std::function<Tcourb(Mesh::Vertex_index)> func = [&] (Mesh::vertex_index v) {return mymq.princip_curv(v) ;} ;
    mymc.set_courbs(func) ;
    mymc.set_K_colors() ;
    mymn.output_mesh("mesh_with_normal.off") ;
    mymc.output_mesh("mesh_with_K.off", "mesh_with_principal_dirs.off") ;


    cout << "----- Writting to off file a quadrtic patch" << endl ;
    // Tirage aléatoire de points
    int nv = m.vertices().size() ;

    std::random_device rd; // obtain a random number from hardw    are
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, nv-1); // define the range

    int ni = distr(gen) ;
    Mesh::Vertex_index vi(ni) ;
    MyQuad *qp = mymq.q[vi] ;
    Mesh &mi(qp->local_mesh(.1, 8)) ;
    CGAL::IO::write_polygon_mesh("sample_quad.off", mi) ;
	
	return 0;
}

