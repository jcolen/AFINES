/*
 *  spring.h
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __LINK_H_INCLUDED__
#define __LINK_H_INCLUDED__


//=====================================
// forward declared dependencies
class filament;

//=====================================
//included dependences
#include "globals.h"

//=====================================
//spring class
class spring
{
    public:
        spring();
        
        spring(double len, double stiffness, double max_ext, filament* f, array<int, 2> aindex, array<double, 2> fov, array<int, 2> nq);
        
        virtual ~spring();

        array<double, 2> get_hx();
        
        array<double, 2> get_hy();
        
        double get_kl();
        
        double get_length();
        
        double get_length_sq();
       
        double get_l0();
        
        double get_fene_ext();
        
        double get_stretching_energy();
        
        double get_xcm();
        
        double get_ycm();
        
        string to_string();
        
        string write(string bc, double shear_dist);
        
        void step();
        
        void step(string bc, double shear_dist);
        
        void filament_update();
        
        bool operator==(const spring& that);    
        
        bool is_similar(const spring& that);    

        void update_force(string bc, double shear_dist);
        
        void update_force_fraenkel_fene(string bc, double shear_dist);
        
        double get_stretching_energy_fene(string bc, double shear_dist);
        
        void update_force_marko_siggia(string bc, double shear_dist, double kToverA);

        array<double,2> get_force();

        void set_aindex1(int i);
        
        double get_distance_sq(string bc, double shear_dist, double xp, double yp);

        double get_int_angle(double xp, double yp);
        
//      array<double,2> get_intpoint(string bc, double shear_dist, double xp, double yp);
        array<double,2> get_intpoint();

        void calc_intpoint(string bc, double shear_dist, double xp, double yp);

		bool get_line_intersect(string bc, double delrx, spring * s2);

		double get_r_c(string bc, double delrx, double x, double y);

		array<double, 2> get_point();

        vector<array<int,2> > get_quadrants();
       
        void quad_update(string bc, double shear_dist);
        
        array<double, 2> get_direction();
        
        array<double, 2> get_disp();
        
        array<double, 2> get_neg_disp();

    protected:

        double xcm, ycm, l0, kl, max_ext, eps_ext, llen, llensq;//, force;
       
        array<double,2> fov, hx, hy;
        array<double, 2> disp, force, intpoint, point, direc;

        array<int, 2> nq, half_nq, aindex;
         
        filament *fil;
        
        vector< array<int,2> > quad; //vector of two vectors(x and y quadrants) of integers
};
#endif
