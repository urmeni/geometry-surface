//
//  mesh_functions.cpp
//  test1
//
//  Created by Bac Alexandra on 17/02/2022.
//

#include <algorithm>
#include "mesh_functions.hpp"

vector<Mesh::Vertex_index> first_vneighbors(Mesh::Vertex_index v, Mesh &m)
{
	vector<Mesh::Vertex_index> res ;
	CGAL::Vertex_around_target_circulator<Mesh> vbegin(m.halfedge(v),m), vdone(vbegin) ;
	do {
		res.push_back(*vbegin) ;
		++vbegin ;
	} while (vbegin != vdone) ;
	return res ;
}

CGAL::Color color_ramp(double vmin, double vmax, double mean, double stdev, double v, CGAL::Color col1, CGAL::Color col2)
{
	stdev = stdev/2 ;
	if (vmin < mean-stdev)
		vmin = mean-stdev ;
	if (vmax > mean+stdev)
		vmax = mean+stdev ;
	
	if (v < vmin)
		v = vmin ;
	else if (v > vmax)
		v = vmax ;
	double per ;
	CGAL::Color c ;
	CGAL::Color col_first, col_second ;
	if (v <= mean)
	{
		per = (v-vmin)/(mean-vmin) ;
		col_first = col1 ;
		col_second = CGAL::white() ;
	}
	else
	{
		per = (v-mean)/(vmax-mean) ;
		col_first = CGAL::white() ;
		col_second = col2 ;
	}
	
	c.red() = floor(col_first.red() + per*(col_second.red()-col_first.red())) ;
	c.green() = floor(col_first.green() + per*(col_second.green()-col_first.green())) ;
	c.blue() = floor(col_first.blue() + per*(col_second.blue()-col_first.blue())) ;
	return c ;
}


