/*------------------------------------------------------------------
 filament_ensemble.cpp : container class for filaments
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "globals.h"
//#include "Link.h"
#include "filament_ensemble.h"
//bead network class

#define SQRT2 1.41421356

filament_ensemble::filament_ensemble(){}

 
filament_ensemble::~filament_ensemble(){ 
    cout<<"DELETING FILAMENT_ENSEMBLE\n";
    
    int s = network.size();
    
    for (int x = 0; x < nq[0]; x++){
        for (int y = 0; y < nq[1]; y++){
            delete springs_per_quad[x]->at(y);
        }
        delete springs_per_quad[x];
        //delete n_springs_per_quad[x];
    }
    
    for (int i = 0; i < s; i++){
        delete network[i];
    }
    
}


filament * filament_ensemble::get_filament(int index)
{
    return network[index];
}


void filament_ensemble::turn_quads_off()
{
    quad_off_flag = true;
}


void filament_ensemble::nlist_init_serial()
{
    for (int x = 0; x < nq[0]; x++){
        springs_per_quad.push_back(new vector< vector<array<int, 2> >* >(nq[1]));
        for (int y = 0; y < nq[1]; y++){
            springs_per_quad[x]->at(y) = new vector<array<int, 2> >();
        }
    }
}

 
void filament_ensemble::quad_update_serial()
{
    int n_quads, net_sz = int(network.size());
    vector<vector<array<int, 2> > > q;
    int x, y;

    //initialize all quadrants to have no springs
    for (x = 0; x < nq[0]; x++){
        for (y = 0; y < nq[1]; y++){
            springs_per_quad[x]->at(y)->clear();
        }
    }
    
    for (int f = 0; f < net_sz; f++){
        q = network[f]->get_quadrants();
        for (int l = 0; l < network[f]->get_nsprings(); l++){
            n_quads = int(q[l].size());
            for (int i = 0; i < n_quads; i++){
                x = q[l][i][0];
                y = q[l][i][1];
                springs_per_quad[x]->at(y)->push_back({{f,l}});
                
            }
        }
    }

}

//given a motor position, and a quadrant
//update the map of {{f, l}} -- > dist
void filament_ensemble::update_dist_map(set<pair<double, array<int,2>>>& t_map, const array<int, 2>& mq, double x, double y){
    
    array<int, 2> fl;
    double dist_sq;
    
    for (int i = 0; i < int(springs_per_quad[mq[0]]->at(mq[1])->size()); i++){

        fl = springs_per_quad[mq[0]]->at(mq[1])->at(i); //fl  = {{filament_index, spring_index}}

        if (fls.find(fl) == fls.end()){
            network[fl[0]]->get_spring(fl[1])->calc_intpoint(network[fl[0]]->get_BC(), delrx, x, y); //calculate the point on the spring closest to (x,y)
            dist_sq = network[fl[0]]->get_spring(fl[1])->get_distance_sq(network[fl[0]]->get_BC(), delrx, x, y); //store the distance to that point
            //cout<<"\nDEBUG : dist = "<<dist;

            t_map.insert(pair<double, array<int, 2> >(dist_sq, fl));
            fls.insert(fl);
        }
    }

}

//given motor head position, return a map between  
//  the INDICES (i.e., {i, j} for the j'th spring of the i'th filament)
//  and their corresponding DISTANCES to the spring at that distance 

set<pair<double, array<int, 2>>> filament_ensemble::get_dist(double x, double y)
{
    fls.clear();
    set<pair<double, array<int, 2>>> t_map;
    
    int mqx = coord2quad(fov[0], nq[0], x);
    int mqy = coord2quad(fov[1], nq[1], y);
    
    update_dist_map(t_map, {{mqx, mqy}}, x, y);
    return t_map;
}


double filament_ensemble::get_llength(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_length();
}


array<double,2> filament_ensemble::get_end(int fil, int spring)
{
    return {{network[fil]->get_spring(spring)->get_hx()[1] , network[fil]->get_spring(spring)->get_hy()[1]}};
}


array<double,2> filament_ensemble::get_force(int fil, int bead)
{
    return network[fil]->get_bead(bead)->get_force();
}


array<double,2> filament_ensemble::get_direction(int fil, int spring)
{
    return network[fil]->get_spring(spring)->get_direction();
}

 
void filament_ensemble::set_y_thresh(double y)
{
    for (unsigned int f = 0; f < network.size(); f++) network[f]->set_y_thresh(y);
}

 
 
double filament_ensemble::get_stretching_energy(){
    return pe_stretch;
}

 
double filament_ensemble::get_bending_energy(){
    return pe_bend;
}

 
vector<int> filament_ensemble::get_broken(){
    return broken_filaments;
}

 
void filament_ensemble::clear_broken(){
    broken_filaments.clear();
}

 
int filament_ensemble::get_nbeads(){
    int tot = 0;
    for (unsigned int f = 0; f < network.size(); f++)
        tot += network[f]->get_nbeads();
    return tot;
}

 
int filament_ensemble::get_nsprings(){
    return this->get_nbeads() - network.size();
}

 
int filament_ensemble::get_nfilaments(){
    return network.size();
}

 
double filament_ensemble::get_delrx(){
    return delrx;
}


void filament_ensemble::update_forces(int f_index, int a_index, double f1, double f2){
    network[f_index]->update_forces(a_index, f1,f2);
}
 

void filament_ensemble::update_filament_stretching(int f){
    vector<filament *> newfilaments = network[f]->update_stretching(t);

    if (newfilaments.size() > 0){ //fracture event occured

        cout<<"\n\tDEBUG: fracturing filament : "<<f;
        filament * broken = network[f];     //store a pointer to the broken filament to delete it with
        network[f] = newfilaments[0];       //replace that pointer with one of the new filaments

        if (newfilaments.size() == 2) network.push_back(newfilaments[1]); //add the second filament to the top of the stack

        broken_filaments.push_back(f);      // record the index, for automatic motor detachment
        delete broken;                      // delete the old filament

    }
}


void filament_ensemble::update_delrx(double drx)
{
    //cout<<"\nDEBUG: SHEARING"; 
    delrx = drx;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_delrx(drx);
    }
}

 
void filament_ensemble::update_d_strain(double g)
{
    //cout<<"\nDEBUG: SHEARING"; 
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_d_strain(g);
    }
}

 
void filament_ensemble::update_shear()
{
    //cout<<"\nDEBUG: SHEARING"; 
    for (unsigned int f = 0; f < network.size(); f++)
    {
        network[f]->update_shear(t);
    }
}

 

/*
Introduce bead-bead interactions. In order to allow formation of a nematic phase,
allow beads to interact via a short-range repulsive force

Loop through quads of the system and find interacting springs
*/
void filament_ensemble::update_filament_interactions()
{
	array <int,2> spring_1, spring_2;
	int f1, f2, l1, l2;	//f = filament number	l = link number
	int x, y, i, j;
	int nsprings_at_quad;
	double par1, par2;

	int dim = this->get_nsprings();
	vector<vector<int>> int_lks (dim, vector<int>(dim, 0));

	//Loop over quadrants of 2D system
	for (x = 0; x < nq[0]; x ++)	{
		for (y = 0; y < nq[1]; y ++)	{
			nsprings_at_quad = int(springs_per_quad[x]->at(y)->size());

			for (i = 0; i < nsprings_at_quad; i ++)	{
				spring_1 = springs_per_quad[x]->at(y)->at(i);

				for (j = i+1; j < nsprings_at_quad; j ++)	{
					spring_2 = springs_per_quad[x]->at(y)->at(j);

					f1 = spring_1[0];
					l1 = spring_1[1];
					f2 = spring_2[0];
					l2 = spring_2[1];

					if (f1 == f2 && abs(l1 - l2) < 2)	continue;

					par1 = f1 * (network[f1]->get_nsprings()) + l1;
					par2 = f2 * (network[f2]->get_nsprings()) + l2;

					if (! int_lks[par1][par2])	{
						int_lks[par1][par2] = 1;
						int_lks[par2][par1] = 1;

						this->update_force_between_filaments(f1, l1, f2, l2);
					}
				}
			}
		}
	}
	int_lks.clear();
}

/*
Calculate forces between actin beads of a pair of filaments
*/
void filament_ensemble::update_force_between_filaments(int f1, int l1, int f2, int l2)	{
	array <spring*, 2> springs;
	array <int, 2> fs, ls;
	array <double, 2> dist, point, hx, hy;
	spring *spring1, *spring2;
	double r=-1., r_c, x1, y1, x2, y2, length, len, r_1, r_2, Fx, Fy;
	int n1, n2, link1, link2, i, j;
	double b = 1. / rmax;

	springs[0] = network[f1]->get_spring(l1);
	springs[1] = network[f2]->get_spring(l2);
	fs[0] = f1; fs[1] = f2;
	ls[0] = l1; ls[1] = l2;

	for (i = 0; i < 1; i ++)	{
		spring1 = springs[i];
		spring2 = springs[(i+1)%2];
		hx = spring1->get_hx();
		hy = spring1->get_hy();
		for (j = 0; j < 1; j ++)	{
			r_c = spring2->get_r_c(BC, delrx, hx[j], hy[j]);
			if (r_c < r || r == -1)	{
				r = r_c;
				length = spring1->get_length();
				point = spring2->get_point();
				x1 = hx[j];
				y1 = hy[j];
				x2 = point[0];
				y2 = point[1];
				len = dist_bc(BC, spring2->get_hx()[0]-x2, spring2->get_hy()[0]-y2, fov[0], fov[1], delrx);
				n1 = fs[(i+1)%2];
				n2 = fs[i];
				link1 = ls[(i+1)%2];
				link2 = ls[i] + j;
			}
		}
	}

	if (r < rmax)	{
		if (springs[0]->get_line_intersect(BC, delrx, springs[1]))	{
			Fx = 2 * kexv / (rmax * SQRT2);
			Fy = 2 * kexv / (rmax * SQRT2);
			network[f1]->update_forces(l1, Fx, Fy);
			network[f1]->update_forces(l1+1, Fx, Fy);
			network[f2]->update_forces(l2, -Fx, -Fy);
			network[f2]->update_forces(l2+1, -Fx, -Fy);
		}
		else	{
			dist = rij_bc(BC, (x2-x1), (y2-y1), fov[0], fov[1], delrx);
			r_2 = len / length;
			r_1 = 1 - r_2;

			Fx = 2 * kexv * dist[0] * b * ((1/r) - b);
			Fy = 2 * kexv * dist[1] * b * ((1/r) - b);

			network[n1]->update_forces(link1, Fx * r_1, Fy * r_1);
			network[n1]->update_forces(link1+1, Fx * r_2, Fy * r_2);
			network[n2]->update_forces(link2, -Fx, -Fy);
		}
		pe_exv += kexv * pow(1 - r * b, 2);
	}

}

/* Overdamped Langevin Dynamics Integrator (Leimkuhler, 2013) */

void filament_ensemble::update()

{      
    int net_sz = network.size();
    // #pragma omp parallel for
  	pe_exv = 0;

	this->update_filament_interactions();

    for (int f = 0; f < net_sz; f++){
      //  if (f==0) cout<<"\nDEBUG: filament updates using "<<omp_get_num_threads()<<" cores"; 
        this->update_filament_stretching(f);
        network[f]->update_bending(t);
        network[f]->update_positions();
    }
    
    this->update_energies();
	this->update_order_parameter();
    
    t += dt;

}


void filament_ensemble::update_energies(){
    pe_stretch = 0;
    pe_bend = 0;
    ke = 0;
    for (unsigned int f = 0; f < network.size(); f++)
    {
        ke += network[f]->get_kinetic_energy();
        pe_bend += network[f]->get_bending_energy();
        pe_stretch += network[f]->get_stretching_energy();
    }
}

void filament_ensemble::update_order_parameter()	{
	/*
		S = \frac{1}{2} < ( 3 cos^2 \theta_i - 1) >
		where theta_i is measured relative to the molecular axis
		
		Here, we will take the molecular axis to be the average orientation 
		in a given quadrant
	*/
	int x, y, i, nsprings_at_quad;
	double cos_t;
	array<double, 2> avg_dir, dir;
	array<int, 2> spring_info;
	
	order_parameter = 0;
	int num_terms = 0;


	for (x = 0; x < nq[0]; x ++)	{
		for (y = 0; y < nq[1]; y ++)	{
			//Establish average direction in this quadrant
			avg_dir = {0, 0};
			nsprings_at_quad = int(springs_per_quad[x]->at(y)->size());
			if (nsprings_at_quad == 0)	continue;
			for (i = 0; i < nsprings_at_quad; i ++)	{
				spring_info = springs_per_quad[x]->at(y)->at(i);
				dir = network[spring_info[0]]->get_spring(spring_info[1])->get_direction();
				//To handle nematic symmetry, fix x > 0
				if (dir[0] < 0)	{
					avg_dir[0] -= dir[0];
					avg_dir[1] -= dir[1];
				} else	{
					avg_dir[0] += dir[0];
					avg_dir[1] += dir[1];
				}
			}
			avg_dir[0] /= nsprings_at_quad;
			avg_dir[1] /= nsprings_at_quad;

			//Add deviations of each spring in this quadrant
			for (i = 0; i < nsprings_at_quad; i ++)	{
				spring_info = springs_per_quad[x]->at(y)->at(i);
				dir = network[spring_info[0]]->get_spring(spring_info[1])->get_direction();
				//Cos theta is dot product of two unit vectors
				cos_t = dir[0] * avg_dir[0] + dir[1] * avg_dir[1];
				order_parameter += 3 * cos_t * cos_t - 1;
			}
			num_terms += nsprings_at_quad;
		}
	}

	order_parameter /= num_terms;
	order_parameter *= 0.5;
}


vector<vector<double> > filament_ensemble::spring_spring_intersections(double len, double prob){

    vector< vector<double> > itrs;
    array<double, 2> r1, r2, s1, s2, direc;
    pair<double, double> mmx1, mmy1, mmx2, mmy2;
    boost::optional<array<double, 2> > inter;
    string bcf1; 
    for (unsigned int f1 = 0; f1 < network.size(); f1++){
        
        for (int l1 = 0; l1 < network[f1]->get_nsprings(); l1++){

            r1 = {{network[f1]->get_spring(l1)->get_hx()[0], network[f1]->get_spring(l1)->get_hy()[0]}};
            r2 = {{network[f1]->get_spring(l1)->get_hx()[1], network[f1]->get_spring(l1)->get_hy()[1]}};
            bcf1 = network[f1]->get_BC();
            for (unsigned int f2 = f1+1; f2 < network.size(); f2++){
                
                for (int l2 = 0; l2 < network[f2]->get_nsprings(); l2++){

                    if (f1 == f2 && fabs(double(l1) - double(l2)) < 2){ //springs should be at least two away to get crosslinked
                        continue;
                    }

                    s1 = {{network[f2]->get_spring(l2)->get_hx()[0], network[f2]->get_spring(l2)->get_hy()[0]}};
                    s2 = {{network[f2]->get_spring(l2)->get_hx()[1], network[f2]->get_spring(l2)->get_hy()[1]}};

                    inter = seg_seg_intersection_bc(bcf1, delrx, fov, r1, r2, s1, s2);
                    
                    if (inter && rng(0,1) <= prob){
                        direc = network[f2]->get_spring(l2)->get_direction();
                        itrs.push_back({inter->at(0), inter->at(1), len*direc[0], len*direc[1], 
                                double(f1), double(f2), double(l1), double(l2)}); 
                    }
                }
            }
        }
    }
    return itrs;
}


void filament_ensemble::write_beads(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_beads(i);
    } 
}

 
void filament_ensemble::write_springs(ofstream& fout)
{
    for (unsigned int i=0; i<network.size(); i++) {
        fout<<network[i]->write_springs(i);
    } 
}

 
void filament_ensemble::write_thermo(ofstream& fout){
    for (unsigned int f = 0; f < network.size(); f++)
        fout<<network[f]->write_thermo(f);
    
}

 
void filament_ensemble::print_network_thermo(){
    cout<<"\nAll Fs\t:\tKE = "<<ke<<"\tPEs = "<<pe_stretch<<"\tPEb = "<<pe_bend<<"\tPEe = "<<pe_exv<<"\tTE = "<<(ke+pe_stretch+pe_bend);
	cout<<"\n\tOrder Parameter = "<<order_parameter;
}

 
////////////////////////////////////////
///SPECIFIC FILAMENT IMPLEMENTATIONS////
////////////////////////////////////////

filament_ensemble::filament_ensemble(int npolymer, int nbeads_min, int nbeads_extra, double nbeads_extra_prob, 
        array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, double spring_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, 
        double frac_force, string bc, double seed, double RMAX, double A) {
    
    fov = myfov;
    view[0] = 1;//(fov[0] - 2*nbeads*len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nbeads*len)/fov[1];
    nq = mynq;
    half_nq = {{nq[0]/2, nq[1]/2}};
    
    double nbeads_mean = nbeads_min + nbeads_extra*nbeads_extra_prob;
    
    visc=vis;
    spring_rest_len = spring_len;
    dt = delta_t;
    temperature = temp;
    shear_stop = 1e10;
    shear_dt = dt;
    t = 0;
    delrx = 0;
	rmax = RMAX;
	kexv = A;
	BC = bc;
    
    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }
    
    
    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Avg number of monomers per filament:"<<nbeads_mean<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<rad<<"\n"; 
   
    int nbeads = 0;
    binomial_distribution<int> distribution(nbeads_extra, nbeads_extra_prob);
    default_random_engine generator(seed+2);

    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament(pos_sets[i], nbeads, fov, nq,
                        visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            
            nbeads = nbeads_min + distribution(generator);
            network.push_back(new filament({{x0,y0,phi0}}, nbeads, fov, nq, visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }
    }
    
    //Neighbor List Initialization
    quad_off_flag = false;
    max_springs_per_quad              = npolymer*(nbeads-1);
    max_springs_per_quad_per_filament = nbeads - 1;
    
    //this->nlist_init();
    this->nlist_init_serial();
    
    pe_stretch = 0;
    pe_bend = 0;
	pe_exv = 0;
    ke = 0;

	this->update_order_parameter();
    
    fls = { };
}

filament_ensemble::filament_ensemble(double density, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double rad, double vis, int nbeads, double spring_len, vector<array<double, 3> > pos_sets, double stretching, double ext, double bending, 
        double frac_force, string bc, double seed, double RMAX, double A) {
    
    fov = myfov;
    view[0] = 1;//(fov[0] - 2*nbeads*len)/fov[0];
    view[1] = 1;//(fov[1] - 2*nbeads*len)/fov[1];
    nq = mynq;
    half_nq = {{nq[0]/2, nq[1]/2}};
   

    visc=vis;
    spring_rest_len =spring_len;
    int npolymer=int(ceil(density*fov[0]*fov[1]) / nbeads);
    dt = delta_t;
    temperature = temp;
    shear_stop = 1e10;
    shear_dt = dt;
    t = 0;
    delrx = 0;
	rmax = RMAX;
	kexv = A;
	BC = bc;
    
    if (seed == -1){
        straight_filaments = true;
    }else{
        srand(seed);
    }
    
    
    cout<<"DEBUG: Number of filament:"<<npolymer<<"\n";
    cout<<"DEBUG: Number of monomers per filament:"<<nbeads<<"\n"; 
    cout<<"DEBUG: Monomer Length:"<<rad<<"\n"; 
   

    int s = pos_sets.size();
    double x0, y0, phi0;
    for (int i=0; i<npolymer; i++) {
        if ( i < s ){
            network.push_back(new filament(pos_sets[i], nbeads, fov, nq,
                        visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }else{
            x0 = rng(-0.5*(view[0]*fov[0]),0.5*(view[0]*fov[0])); 
            y0 = rng(-0.5*(view[1]*fov[1]),0.5*(view[1]*fov[1]));
            phi0 =  rng(0, 2*pi);
            //phi0=atan2(1+x0-y0*y0, -1-x0*x0+y0); // this is just the first example in mathematica's streamplot documentation
            network.push_back(new filament({{x0,y0,phi0}}, nbeads, fov, nq, visc, dt, temp, straight_filaments, rad, spring_rest_len, stretching, ext, bending, frac_force, bc) );
        }
    }
    
    //Neighbor List Initialization
    quad_off_flag = false;
    max_springs_per_quad              = npolymer*(nbeads-1);
    max_springs_per_quad_per_filament = nbeads - 1;
    
    //this->nlist_init();
    this->nlist_init_serial();
    
    pe_stretch = 0;
    pe_bend = 0;
	pe_exv = 0;
    ke = 0;
    
	this->update_order_parameter();

    fls = { };
}

filament_ensemble::filament_ensemble(vector<vector<double> > beads, array<double,2> myfov, array<int,2> mynq, double delta_t, double temp,
        double vis, double spring_len, double stretching, double ext, double bending, double frac_force, string bc, double RMAX, double A) {
    
    fov = myfov;

    visc=vis;
    spring_rest_len = spring_len;
    dt = delta_t;
    temperature = temp;
    t = 0;
    delrx = 0;

	rmax = RMAX;
	kexv = A;
	BC = bc;
	pe_exv = 0;

    view[0] = 1;
    view[1] = 1;

    int s = beads.size(), sa, j;
    int fil_idx = 0;
    vector<bead *> avec;
    
    nq = mynq;
    
    for (int i=0; i < s; i++){
        
        if (beads[i][3] != fil_idx && avec.size() > 0){
            
            network.push_back( new filament( avec, fov, nq, spring_rest_len, stretching, ext, bending, delta_t, temp, frac_force, 0, bc) );
            
            sa = avec.size();
            for (j = 0; j < sa; j++) delete avec[j];
            avec.clear();
            
            fil_idx = beads[i][3];
        }
        avec.push_back(new bead(beads[i][0], beads[i][1], beads[i][2], vis));
    }

    sa = avec.size();
    if (sa > 0)
        network.push_back( new filament( avec, fov, nq, spring_rest_len, stretching, ext, bending, delta_t, temp, frac_force, 0, bc) );
    
    for (j = 0; j < sa; j++) delete avec[j];
    avec.clear();
   
    quad_off_flag = false;
    max_springs_per_quad              = beads.size();
    max_springs_per_quad_per_filament = int(ceil(beads.size() / (fil_idx + 1)))- 1;
    //this->nlist_init();
    this->nlist_init_serial();
    this->update_energies();
	this->update_order_parameter();
    
    fls = { };
} 
