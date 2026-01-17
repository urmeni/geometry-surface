//
//  courbures.cpp
//  code_courb
//
//  Created by Bac Alexandra on 07/04/2022.
//

#include <fstream>
#include "courbures.h"
#include "mesh_functions.hpp"
#include "distrib.hpp"


// ////////////////////////////////////////////////////////////////////
// MyMesh_Norm
// ////////////////////////////////////////////////////////////////////

MyMesh_Norm::MyMesh_Norm(MyMesh &m) : MyMesh(m), vnormal(*(new Mesh::Property_map<vertex_descriptor,K::Vector_3>)), vcolor(*(new Mesh::Property_map<vertex_descriptor,CGAL::Color>))
{
	// Allocation de vnormal et initialisation
	bool created_vnorm;
	boost::tie(vnormal, created_vnorm) = _m.add_property_map<vertex_descriptor,K::Vector_3>("v:norm",K::Vector_3(0.,0.,0.));
	assert(created_vnorm);
	normales_locales() ;
    
    // Allocation de vcolor
    bool created_vcol;
    boost::tie(vcolor, created_vcol) = _m.add_property_map<vertex_descriptor,CGAL::Color>("v:color",CGAL::white());
    assert(created_vcol);
    
    // Init local normals
    normales_locales() ;
}

MyMesh_Norm::~MyMesh_Norm() {
    //std::optional<Mesh::Property_map<vertex_descriptor, K::Vector_3>> res_norm = _m.property_map<vertex_descriptor,K::Vector_3>("v:norm");
    // CGAL < 6
	std::pair<Mesh::Property_map<vertex_descriptor,K::Vector_3>, bool> res_norm = _m.property_map<vertex_descriptor,K::Vector_3>("v:norm");
    
    if (res_norm.second)
//    if (res_norm.has_value()) // the property exists
    {
        _m.remove_property_map(vnormal) ;
        delete &vnormal ;
    }
    else
        cerr << "normal property already removed" << endl ;
}

void MyMesh_Norm::normales_locales()
{
	Mesh::Vertex_range vr = _m.vertices() ;
	for (Mesh::Vertex_range::iterator v_it = vr.begin() ; v_it != vr.end(); ++v_it)
	{
        vnormal[*v_it] = compute_normal(*v_it, _m) ;
	}
}

void MyMesh_Norm::set_fixed_colors()
{
    Mesh::Vertex_range vr = _m.vertices() ;

    int i = 0 ;
    for (Mesh::Vertex_range::iterator v_it = vr.begin(); v_it != vr.end(); ++v_it)
    {
        i++ ;
        if (i % 2)
            vcolor[*v_it] = CGAL::green() ;
        else
            vcolor[*v_it] = CGAL::blue() ;
    }
}

void MyMesh_Norm::output_mesh(const string &filename)
{
    CGAL::IO::write_polygon_mesh(filename, _m, //CGAL::parameters::vertex_color_map(vcolor));//
                                 CGAL::parameters::vertex_normal_map(vnormal)) ;
}

// ////////////////////////////////////////////////////////////////////
// MyMesh_Courb
// ////////////////////////////////////////////////////////////////////


MyMesh_Courb::MyMesh_Courb(MyMesh_Norm &mn) : _mn(mn)
{
    // Allocation of d1,d2
    bool created_dir;
    boost::tie(d1, created_dir) = _mn._m.add_property_map<vertex_descriptor,K::Vector_3>("v:d1",K::Vector_3(0,0,0));
    assert(created_dir);
    boost::tie(d2, created_dir) = _mn._m.add_property_map<vertex_descriptor,K::Vector_3>("v:d2",K::Vector_3(0,0,0));
    assert(created_dir);
    
    // Allocation of k1, k2
    bool created_lam;
    boost::tie(k1, created_lam) = _mn._m.add_property_map<vertex_descriptor,double>("v:k1",0.);
    assert(created_lam);
    boost::tie(k2, created_lam) = _mn._m.add_property_map<vertex_descriptor,double>("v:k2",0.);
    assert(created_lam);
}

void MyMesh_Courb::set_K_colors()
{
    distrib<double> dist ;
    
    // Insertion des courbures dans la distribution
    Mesh::Vertex_range vr = _mn._m.vertices() ;
    for (Mesh::Vertex_range::iterator v_it = vr.begin(); v_it != vr.end(); ++v_it)
    {
        dist.add_data(getK(*v_it)) ;
    }
    dist.update_stdev() ;
    
    std::cout << "K min : " << dist.get_min() << " - max : " << dist.get_max() << std::endl ;
    std::cout << "K mean : " << dist.get_mean() << " - sigma : " << dist.get_stdev() << std::endl;

    for (Mesh::Vertex_range::iterator v_it = vr.begin(); v_it != vr.end(); ++v_it)
    {
        _mn.vcolor[*v_it] = color_ramp(dist.get_min(), dist.get_max(), dist.get_mean(), dist.get_stdev(), getK(*v_it), CGAL::blue(), CGAL::red()) ;
    }
}

void MyMesh_Courb::set_vcourb (Mesh::Vertex_index v, std::function<Tcourb(Mesh::Vertex_index)> func)
{
    Tcourb courb (func(v)) ;
    d1[v] = courb.d1 ;
    d2[v] = courb.d2 ;
    k1[v] = courb.k1 ;
    k2[v] = courb.k2 ;
}

void MyMesh_Courb::set_courbs (std::function<Tcourb(Mesh::Vertex_index)> func)
{
    Mesh::Vertex_range vr = _mn._m.vertices() ;
    for (Mesh::Vertex_range::iterator v_it = vr.begin() ; v_it != vr.end(); ++v_it)
    {
        set_vcourb(*v_it, func) ;
    }
}

void MyMesh_Courb::output_mesh(const string &filename, const string &filename_field)
{
    const double scale = .02 ;
    Mesh fields ;
    // Ajout d'une couleur par face
    Mesh::Property_map<face_descriptor,CGAL::Color> fcolor;
    bool created_fcol;
    boost::tie(fcolor, created_fcol) = fields.add_property_map<face_descriptor,CGAL::Color>("f:color",CGAL::white());
    assert(created_fcol);

    Mesh::Vertex_range vr = _mn._m.vertices() ;
    for (Mesh::Vertex_range::iterator v_it = vr.begin(); v_it != vr.end(); ++v_it)
    {
        K::Point_3 p, p1, p2 ;
        p = _mn._m.point(*v_it) ;
        p1 = p+ scale *d1[*v_it] ;
        p2 = p+ scale *d2[*v_it] ;
        arrow_edge(p, p1, _mn.vnormal[*v_it], CGAL::red(), fcolor, fields) ;
        arrow_edge(p, p2, _mn.vnormal[*v_it], CGAL::green(), fcolor, fields) ;
    }
    cout << "writting mesh ..." << endl ;
    CGAL::IO::write_polygon_mesh(filename, _mn._m, CGAL::parameters::vertex_color_map(_mn.vcolor)) ;
    cout << "writting dirs mesh ..." << endl ;
    CGAL::IO::write_polygon_mesh(filename_field, fields, CGAL::parameters::face_color_map(fcolor)) ;
}


void MyMesh_Courb::arrow_edge(K::Point_3 p1, K::Point_3 p2, K::Vector_3 n, CGAL::Color c, Mesh::Property_map<face_descriptor,CGAL::Color> &fcolor, Mesh &m)
{
    const double scale = .005 ;
    K::Point_3 bas(p1-scale*n), haut(p1+scale*n) ;
    Mesh::Vertex_index vb = m.add_vertex(bas), vh = m.add_vertex(haut), v2 = m.add_vertex(p2) ;
    Mesh::Face_index f = m.add_face(vb,vh,v2) ;
    fcolor[f] = c ;
}

// ////////////////////////////////////////////////////////////////////
// MyMesh_Quad
// ////////////////////////////////////////////////////////////////////


MyMesh_Quad::MyMesh_Quad(MyMesh_Norm &mn) : _mn(mn)
{
    // Allocation of vquads
    bool created_vquads;
    boost::tie(q, created_vquads) = _mn._m.add_property_map<vertex_descriptor,MyQuad *>("v:quads",NULL);
    assert(created_vquads);
}

MyMesh_Quad::~MyMesh_Quad ()
{
	Mesh::Vertex_range vr = _mn._m.vertices() ;
	for (Mesh::Vertex_range::iterator v_it = vr.begin(); v_it != vr.end(); ++v_it)
	{
        if (q[*v_it] != NULL)
            delete q[*v_it] ;
	}
}

void MyMesh_Quad::fit_quad(Mesh::Vertex_index v)
{
	std::vector<Mesh::Vertex_index> neigh = get_two_vneighborhood(v, _mn._m) ;
	// Ajout du sommet v
	neigh.push_back(v);

	// Calcul de la matrice de changement de base
	K::Vector_3 n = _mn.vnormal[v] ;
	K::Point_3 p = _mn._m.point(v);

	Eigen::Vector3d ne(n[0], n[1], n[2]), Oz(0,0,1), axis;
	Eigen::Vector3d p_e(p[0], p[1], p[2]), pi_e ;

	axis = ne.cross(Oz) ;
	double sina = axis.norm(), cosa = ne.dot(Oz), angle ;
	if (sina >= 0)
		angle = acos(cosa) ;
	else
		angle = -acos(cosa) ;
	axis = axis.normalized() ;

	Eigen::AngleAxisd r(angle, axis) ;
	Eigen::Translation3d t(-p_e[0], -p_e[1], -p_e[2]) ;
	Eigen::Transform<double, 3, Eigen::Affine> ch_base = r * t ;


	// Calcul de la matrice / vecteur de moindres carrés linéaires (Eigen)
	// A : (x_i^2 x_i.y_i y_i^2 x_i y_i)
	// B : (z_i)

	if (neigh.size() >= 5)
	{
		int n(neigh.size()) ;
		Eigen::MatrixXd A(n,5);
		Eigen::VectorXd B(n);

		for(int i=0; i<neigh.size(); i++)
		{
			K::Point_3 p = _mn._m.point(neigh.at(i)) ;
			pi_e << p[0], p[1], p[2] ;
			pi_e = ch_base * pi_e ; // Application du changement de base
			
            // Remplissage de A et B
            A(i,0) = pi_e(0) * pi_e(0) ; // x^2
            A(i,1) = pi_e(0) * pi_e(1) ; // xy
            A(i,2) = pi_e(1) * pi_e(1) ; // y^2
            A(i,3) = pi_e(0) ; // x
            A(i,4) = pi_e(1) ; // y
            B(i) = pi_e(2) ; // z
		}

		// Résolution aux moindres carrés par méthode SVD
		Eigen::VectorXd coef(5) ;
		coef = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B) ;
		
		// Store the quadratic surface in q
		q[v] = new MyQuad(coef[0], coef[1], coef[2], coef[3], coef[4], r, t) ;
	}
	else
	{
		std::cout << "Quad fitting : not enough neighbors" ;
		throw "Quad fitting : not enough neighbors" ;
	}
}

void MyMesh_Quad::fit_quads()
{
	Mesh::Vertex_range vr = _mn._m.vertices() ;
	for (Mesh::Vertex_range::iterator v_it = vr.begin(); v_it != vr.end(); ++v_it)
	{
		fit_quad(*v_it) ;
	}
}

Tcourb MyMesh_Quad::princip_curv (Mesh::Vertex_index v) const
{
	Tcourb curv ;
	MyQuad  qv(*(q[v])) ;
	std::pair<Eigen::VectorXd, Eigen::MatrixXd> res = qv.get_principal_0 () ;
	Eigen::VectorXd lam(2) ;
	Eigen::MatrixXd dirs(3,2) ;
	lam = res.first ;
	dirs = res.second ;
	curv.k1 = lam(0) ;
	curv.k2 = lam(1) ;
	curv.d1 = Eig2Cgal_vec(dirs.col(0)) ;
	curv.d2 = Eig2Cgal_vec(dirs.col(1)) ;
	return curv ;
}

// ////////////////////////////////////////////////////////////////////
// Tools
// ////////////////////////////////////////////////////////////////////

std::vector<Mesh::Vertex_index> get_two_vneighborhood(Mesh::Vertex_index v, Mesh &m)
{
	// Ajout d'un flag booléen
	Mesh::Property_map<vertex_descriptor,bool> vflag;
	bool created_vflag;
	boost::tie(vflag, created_vflag) = m.add_property_map<vertex_descriptor,bool>("v:flag",false);
	assert(created_vflag);
	
	// Circulateur sur le premier cercle
	std::vector<Mesh::Vertex_index> neigh = first_vneighbors(v, m), neigh2 ;

	// Parcours du premier cercle et ajout du second cercle par circulateurs
	for (int i=0; i<neigh.size(); i++)
	{
		CGAL::Vertex_around_target_circulator<Mesh> v_it(m.halfedge(neigh.at(i)),m), v_end(v_it) ;
		do
		{
			if (!vflag[*v_it]) // sommet non encore rencontré
			neigh2.push_back(*v_it) ; // ajout du point à la liste
			vflag[*v_it] = true ;
			++v_it ;
		} while (v_it != v_end) ;
	}

	// Concaténation des deux cercles
	neigh.insert(neigh.end(), neigh2.begin(), neigh2.end()) ;

	m.remove_property_map(vflag) ;

	return neigh ;
}

K::Vector_3 compute_normal(Mesh::Vertex_index v, Mesh &m)
{
	K::Vector_3 n(0.,0.,0.) ;

    // Visit faces aroung v
    CGAL::Face_around_target_circulator<Mesh> fit(m.halfedge(v), m), fit_end(fit) ;
    do
    {
        Mesh::Face_index fid = *fit ;
        // Visit vertices of the face fid
        CGAL::Vertex_around_face_iterator<Mesh> vit(m.halfedge(fid), m), vit_end(vit) ;
        std::vector<Mesh::Vertex_index> neigh ;
        do
        {
            neigh.push_back(*vit) ;
            ++vit ;
        } while (*vit != *vit_end) ;
        
        // Compute the normal of 3 first points
        if (neigh.size()>=3)
        {
            n += CGAL::normal(m.point(neigh.at(0)), m.point(neigh.at(1)), m.point(neigh.at(2))) ;
        }
        else
        {
            cerr << "face " << fid << "has less than 3 vertices !!!" << endl ;
        }
        ++fit ;
    } while (fit != fit_end) ;
    // Normalisation
    n /= sqrt(n.squared_length()) ;
    return n ;
}
