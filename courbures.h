//
//  courbures.h
//  code_courb
//
//  Created by Bac Alexandra on 07/04/2022.
//

#ifndef courbures_h
#define courbures_h

#include <stdio.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include "starter_defs.h"

// Types
typedef struct zcourb {
	K::Vector_3 d1, d2 ;
	double k1, k2 ; } Tcourb ;

// ////////////////////////////////////////////////////////////////////
// MyQuad
// ////////////////////////////////////////////////////////////////////

class MyQuad
{
	double _coefs[5] ; // a_0 x^2 + a1 xy + a2 y^2 + a3 x + a4 y + a5
	Eigen::AngleAxisd _r ;
	Eigen::Translation3d _t ;
public:
	MyQuad(double *data,
		   const Eigen::AngleAxisd &r = Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)),
		   const Eigen::Translation3d &t = Eigen::Translation3d(0,0,0))
		: _r(r), _t(t)
	{
		for (int i=0; i<5; i++)
			_coefs[i] = data[i] ;
	}
	MyQuad(const Eigen::VectorXd &v,
		   const Eigen::AngleAxisd &r = Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)),
		   const Eigen::Translation3d &t = Eigen::Translation3d(0,0,0))
		: _r(r), _t(t)
	{
		for (int i=0; i<5; i++)
			_coefs[i] = v[i] ;
	}
	MyQuad(double a0=0, double a1=0, double a2=0, double a3=0, double a4=0,
		   const Eigen::AngleAxisd &r = Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)),
		   const Eigen::Translation3d &t = Eigen::Translation3d(0,0,0))
		: _r(r), _t(t)
	{
		_coefs[0] = a0 ;
		_coefs[1] = a1 ;
		_coefs[2] = a2 ;
		_coefs[3] = a3 ;
		_coefs[4] = a4 ;
	}
	MyQuad(const MyQuad &q) : _r(q._r), _t(q._t)
	{
		for (int i=0; i<5; i++)
			_coefs[i] = q._coefs[i] ;
	}

	double & operator[] (int i) {return _coefs[i] ; }

	double quad_fun (double x, double y)
	{
		return _coefs[0] * x*x + _coefs[1] * x*y + _coefs[2] * y*y + _coefs[3] * x + _coefs[4] * y  ;
	}

	double quad_fun (const K::Point_2 &v)
	{

		return quad_fun(v[0], v[1]) ;
	}
	
	std::pair<Eigen::Vector3d,Eigen::Vector3d> quad_tangent_frame_0 ()
	{
		Eigen::Vector3d t1, t2 ;
		// Base du plan tangent : dérivées partielles (1, 0, a3)' (0,1, a4)'
		t1 << 1., 0., _coefs[3] ;
		t2 << 0., 1., _coefs[4] ;
		return std::pair<Eigen::Vector3d,Eigen::Vector3d>(t1,t2) ;
	}
	
	Eigen::Vector3d quad_norm_0 ()
	{
		Eigen::Vector3d n ;
		// Normale en point P (normalisée)
		n << -_coefs[3], -_coefs[4], 1 ;
		return n = n / n.norm() ;
	}
	
	Eigen::MatrixXd quad_first_fond_0 ()
	{
		// Calcul de la première forme fondamentale
		Eigen::MatrixXd A(2,2) ;
		double E, F, G ;

		// Récupération des pentes à l'origine a3 et a4
		double a3 = _coefs[3] ;
		double a4 = _coefs[4] ;

		// Application de la formule : E = 1 + a3^2 ; F = a3*a4 ; G = 1 + a4^2
		E = 1. + a3*a3 ;
        F = a3*a4 ;
        G = 1. + a4*a4 ;
        
		A << E, F,
			 F,G ;
		return A ;
	}
	
	Eigen::MatrixXd quad_sec_fond_0 ()
	{
		// Calcul de la deuxième forme fondamentale
		Eigen::MatrixXd A(2,2) ;
		double g, l, m, n ;

		// Récupération des coefficients quadratiques a0, a1, a2 ... a4
		double a0 = _coefs[0] ;
		double a1 = _coefs[1] ;
		double a2 = _coefs[2] ;
		double a3 = _coefs[3] ;
		double a4 = _coefs[4] ;

		// Calcul des numérateurs (dérivés secondes)
		l = 2.*a0; // dérivée seconde en u
		m = a1; // dérivée croisée uv
		n = 2.*a2; // dérivée seconde en v

		// Caclul du facteur de normalisation (norme du gradient)
		g = std::sqrt(1.0 + a3*a3 + a4*a4);

		A << l, m,
			 m,n ;

		// Projection sur la normale unitaire
		A /= g ;
		return A ;
	}
	
	Eigen::MatrixXd courb_tens_0 ()
	{
		Eigen::MatrixXd I(2,2), II(2,2), A(2,2) ;
		I = quad_first_fond_0 () ;
		II = quad_sec_fond_0 () ;
		A = I.partialPivLu().solve(II) ;
		return A ;
	}

	std::pair<Eigen::VectorXd,Eigen::MatrixXd> get_eigen_0 (Eigen::MatrixXd K0)
	{
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(K0);
		if (eigensolver.info() != Eigen::Success) abort();
		return std::pair<Eigen::VectorXd,Eigen::MatrixXd>(eigensolver.eigenvalues(),eigensolver.eigenvectors());
	}
	
	std::pair<Eigen::VectorXd, Eigen::MatrixXd> get_principal_0 ()
	{
		Eigen::MatrixXd K0(2,2) ;
		K0 << courb_tens_0 () ;
		std::pair<Eigen::VectorXd,Eigen::MatrixXd> eig = get_eigen_0 (K0) ;
		Eigen::VectorXd lambda(2) ;
		lambda << eig.first ;
		Eigen::MatrixXd eig_vects(2,2) ;
		eig_vects << eig.second ;
		
		// Passage du repère du plan -> repère locale de la surface
		// Base de plan tangent
		std::pair<Eigen::Vector3d,Eigen::Vector3d> tmp = quad_tangent_frame_0 () ;
		Eigen::Vector3d t1, t2 ;
		t1 << tmp.first ;
		t2 << tmp.second ;
		// Passage des vecteurs principaux dans cette base
		Eigen::Vector3d d1_loc, d2_loc ;
		d1_loc = eig_vects(0,0) * t1 + eig_vects(1,0) * t2 ;
		d2_loc = eig_vects(0,1) * t1 + eig_vects(1,1) * t2 ;
		// Passage dans le repère global
		Eigen::Transform<double, 3, Eigen::Affine> ch_base = _r * _t ;
		ch_base = ch_base.inverse() ;
		Eigen::Vector3d d1, d2 ;
		d1 = ch_base * d1_loc ;
		d2 = ch_base * d2_loc ;
		d1 = d1/d1.norm() ;
		d2 = d2/d2.norm() ;
		
		// Préparation des retours
		Eigen::MatrixXd D(3,2) ;
		D.col(0) = d1 ;
		D.col(1) = d2 ;
		
		return std::pair<Eigen::VectorXd, Eigen::MatrixXd>(lambda, D) ;
	}
	
	Mesh & local_mesh(double alpha, int npix)
	{
		Mesh &m_loc = *(new Mesh) ;
		double step = alpha/npix ;
		double x,y,z ;
		Mesh::Vertex_index verts[npix+1][npix+1] ;
		Eigen::Transform<double, 3, Eigen::Affine> ch_base = _r * _t ;
		ch_base = ch_base.inverse() ;
		
		// Calcul des coordonnées des points dans le repère local -> passage dans la base globale puis création des sommets correspondants
		for (int i=0; i<=npix; ++i)
		{
			for (int j=0; j<=npix; ++j)
			{
				x = i*step-alpha/2. ;
				y = j*step-alpha/2. ;
				z = quad_fun(x,y) ;
				Eigen::Vector3d p(x,y,z) ;
				// Passage dans la base globale
				p = ch_base * p ;
				// Ajout du sommet
				verts[i][j] = m_loc.add_vertex(K::Point_3(p(0),p(1),p(2))) ;
			}
		}
		// Création du maillage
		for (int i=0; i<npix; ++i)
		{
			for (int j=0; j<npix; ++j)
			{
				// Pour chaque quadrangle donne deux triangles
                m_loc.add_face(verts[i][j],verts[i+1][j], verts[i][j+1]) ;
                m_loc.add_face(verts[i+1][j],verts[i+1][j+1], verts[i][j+1]) ;
			}
		}
		return m_loc ;
	}
};

// ////////////////////////////////////////////////////////////////////
// MyMesh
// ////////////////////////////////////////////////////////////////////

class MyMesh {
public:
	Mesh &_m ;
	MyMesh(Mesh &m) : _m(m) {}
} ;

// WARNING
// Property maps are linked to a given mesh so they cannot be duplicated!

// ////////////////////////////////////////////////////////////////////
// MyMesh_Norm
// ////////////////////////////////////////////////////////////////////

class MyMesh_Norm : public MyMesh
{
public:
	// Data
	Mesh::Property_map<vertex_descriptor,K::Vector_3> & vnormal;
    Mesh::Property_map<vertex_descriptor,CGAL::Color> & vcolor;
	
	// Constructors
	MyMesh_Norm(MyMesh &m) ;
    ~MyMesh_Norm() ;
	
    // Methods
	void normales_locales() ;
    virtual void set_fixed_colors() ;
    virtual void output_mesh(const std::string &filename) ;
} ;

// ////////////////////////////////////////////////////////////////////
// MyMesh_Courb
// ////////////////////////////////////////////////////////////////////

class MyMesh_Courb
{
private:
    MyMesh_Norm & _mn ;
    
    void arrow_edge(K::Point_3 p1, K::Point_3 p2, K::Vector_3 n, CGAL::Color c, Mesh::Property_map<face_descriptor,CGAL::Color> &fcolor, Mesh &m) ;
public:
	// Data
	Mesh::Property_map<vertex_descriptor,K::Vector_3> d1;
	Mesh::Property_map<vertex_descriptor,K::Vector_3> d2;
	Mesh::Property_map<vertex_descriptor,double> k1;
	Mesh::Property_map<vertex_descriptor,double> k2;
	
	// Methods
    MyMesh_Courb(MyMesh_Norm &mn) ;
	double getK (Mesh::Vertex_index v) { return k1[v]*k2[v] ;}
	double getH (Mesh::Vertex_index v) { return (k1[v]+k2[v])/2. ;}
    void set_vcourb (Mesh::Vertex_index v, std::function<Tcourb(Mesh::Vertex_index)> func) ;
    void set_courbs (std::function<Tcourb(Mesh::Vertex_index)> func) ;
    virtual void set_K_colors() ;
    void output_mesh(const std::string &filename, const std::string &filename_field) ;
} ;

// ////////////////////////////////////////////////////////////////////
// MyMesh_Quad
// ////////////////////////////////////////////////////////////////////

class MyMesh_Quad
{
private:
    MyMesh_Norm & _mn ;
public:
	// Data
	Mesh::Property_map<vertex_descriptor,MyQuad *> q ;
	
	// Methods
    MyMesh_Quad(MyMesh_Norm &mn) ;
	~MyMesh_Quad() ;
	
	void fit_quad(Mesh::Vertex_index v) ;
	void fit_quads() ;
	Tcourb princip_curv (Mesh::Vertex_index v) const ;
} ;


// ////////////////////////////////////////////////////////////////////
// Tools
// ////////////////////////////////////////////////////////////////////


inline K::Vector_3 Eig2Cgal_vec(Eigen::VectorXd v) { return K::Vector_3 (v(0),v(1),v(2)) ;}

std::vector<Mesh::Vertex_index> get_two_vneighborhood(Mesh::Vertex_index v, Mesh &m) ;
K::Vector_3 compute_normal(Mesh::Vertex_index v, Mesh &m) ;

#endif /* courbures_hpp */
