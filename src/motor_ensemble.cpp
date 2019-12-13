/*------------------------------------------------------------------
 motor_ensemble.cpp : container class for motors
 
 Copyright (C) 2016 
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details. 
-------------------------------------------------------------------*/

#include "globals.h"
#include "motor.h"
#include "motor_ensemble.h"
#include "filament_ensemble.h"

//motor_ensemble class

motor_ensemble::motor_ensemble(double mdensity, array<double, 2> myfov, double delta_t, double temp, 
        double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio, 
        double ron, double roff, double rend, 
        double fstall, double rcut,
        double vis, vector<array<double,3> > positions, string BC) {
    
    fov = myfov;
    mld =mlen;
    gamma = 0;
    tMove=0;//10;
    f_network=network;
    v = v0;
    
    ke = 0;
    pe = 0;
	n_attached = 0;

	stress = {0, 0, 0, 0};

    int nm = int(ceil(mdensity*fov[0]*fov[1]));
    cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";

    double alpha = 1, motorx, motory, mang;
    array<double, 3> motor_pos; 
    for (int i=0; i< nm; i++) {
        
        if ((unsigned int)i < positions.size()){
            motorx = positions[i][0];
            motory = positions[i][1];
            mang   = positions[i][2];
        }else{
            motorx = rng(-0.5*(fov[0]*alpha-mld),0.5*(fov[0]*alpha-mld));
            motory = rng(-0.5*(fov[1]*alpha-mld),0.5*(fov[1]*alpha-mld));
            mang   = rng(0,2*pi);
        }
        motor_pos = {{motorx, motory, mang}};

        n_motors.push_back(new motor( motor_pos, mld, f_network,{{0, 0}}, {{-1,-1}}, {{-1,-1}}, fov, delta_t, temp, 
                    v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, BC));
        
    }
}


motor_ensemble::motor_ensemble(vector<vector<double> > motors, array<double, 2> myfov, double delta_t, double temp, 
        double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio, 
        double ron, double roff, double rend, 
        double fstall, double rcut,
        double vis, string BC) {
    
    fov = myfov;
    mld = mlen;
    gamma = 0;
    tMove = 0;
    f_network=network;
    v = v0;

    ke = 0;
    pe = 0;
	n_attached = 0;
   
	stress = {0, 0, 0, 0};

    int nm = motors.size();
    cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";

    array<double, 4> motor_pos;
    array<int, 2> f_index, l_index, state;

    for (int i=0; i< nm; i++) {
        
        motor_pos = {{motors[i][0], motors[i][1], motors[i][2], motors[i][3]}};
        
        f_index = {{int(motors[i][4]), int(motors[i][5])}};
        l_index = {{int(motors[i][6]), int(motors[i][7])}};

        state = {{f_index[0] == -1 && l_index[0] == -1 ? 0 : 1, f_index[1] == -1 && l_index[1] == -1 ? 0 : 1}};  

		if (state[0] == 1 && state[1] == 1)	n_attached += 1;

        n_motors.push_back(new motor( motor_pos, mld, f_network, state, f_index, l_index, fov, delta_t, temp, 
                    v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, BC));
    }

    this->update_energies();
	this->update_stress();
}


motor_ensemble::~motor_ensemble( ){ 
    cout<<"DELETING MOTOR ENSEMBLE\n";
    int s = n_motors.size();
    for (int i = 0; i < s; i++){
        delete n_motors[i];
    }
    n_motors.clear();
};


int motor_ensemble::get_nmotors( ){ 
    return n_motors.size();
}

int motor_ensemble::get_nattached()	{
	return n_attached;	
}

void motor_ensemble::kill_heads(int hd){
    for (unsigned int i = 0; i < n_motors.size(); i++)
        n_motors[i]->kill_head(hd);
}

//check if any motors attached to filaments that no longer exist; 
// if they do, detach them
// Worst case scenario, this is a O(m*n*p) function,
// where m is the number of filaments
//       n is the number of rods per filament
//       p is the number of motors
// However: we don't expect to fracture often, 
// so this loop should rarely if ever be accessed.
    


void motor_ensemble::check_broken_filaments()
{
    vector<int> broken_filaments = f_network->get_broken();
    array<int, 2> f_index;
    
    for (unsigned int i = 0; i < broken_filaments.size(); i++){
        
        for(unsigned int j = 0; j < n_motors.size(); j++){
            
            f_index = n_motors[j]->get_f_index();

            if(f_index[0] == broken_filaments[i]){
                n_motors[j]->detach_head_without_moving(0);
                //cout<<"\nDEBUG: detaching head 0 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }

            if(f_index[1] == broken_filaments[i]){
                n_motors[j]->detach_head_without_moving(1);
                //cout<<"\nDEBUG: detaching head 1 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }
        }
    }


}


void motor_ensemble::motor_walk(double t)
{

    this->check_broken_filaments();
    int nmotors_sz = int(n_motors.size());
    //#pragma omp parallel for
    
	n_attached = 0;

    for (int i=0; i<nmotors_sz; i++) {
       
    //    if(i==0) cout<<"\nDEBUG: motor_walk: using "<<omp_get_num_threads()<<" cores";
        array<int, 2> s = n_motors[i]->get_states();
        
        if (t >= tMove){
            
            //Dynamics
            if (s[0] == 0)
                n_motors[i]->brownian_relax(0);
            if (s[1] == 0)
                n_motors[i]->brownian_relax(1);
            
            n_motors[i]->update_angle();
            n_motors[i]->update_force();
            n_motors[i]->filament_update();
            
            //Attachment or Movement + Detachment
            if (s[0] == 0)
                n_motors[i]->attach(0);
            else
                n_motors[i]->step_onehead(0);

            if (s[1] == 0)
                n_motors[i]->attach(1);
            else 
                n_motors[i]->step_onehead(1);

			n_motors[i]->update_stress();

			s = n_motors[i]->get_states();
			if (s[0] == 1 && s[1] == 1)	n_attached += 1; 
        }

    }
    this->update_energies();
    this->update_stress();
}

/* Used for static, contstantly attached, motors -- ASSUMES both heads are ALWAYS attached */

void motor_ensemble::motor_update()
{

    this->check_broken_filaments();
    int nmotors_sz = int(n_motors.size());
    //#pragma omp parallel for
    
    for (int i=0; i<nmotors_sz; i++) {
       
            n_motors[i]->update_position_attached(0);
            n_motors[i]->update_position_attached(1);
            n_motors[i]->update_angle();
            n_motors[i]->update_force();
			n_motors[i]->update_stress();
            //n_motors[i]->update_force_fraenkel_fene();
            n_motors[i]->filament_update();
    
    }
    this->update_energies();
	this->update_stress();
}

void motor_ensemble::motor_write(ostream& fout)
{
    for (unsigned int i=0; i<n_motors.size(); i++) {
        fout<<n_motors[i]->write();
    } 
}


void motor_ensemble::add_motor(motor * m)
{
    n_motors.push_back(m);
}


void motor_ensemble::set_shear(double g)
{
    for (unsigned int i=0; i<n_motors.size(); i++)
        n_motors[i]->set_shear(g);
    
    gamma = g;
}

 
void motor_ensemble::update_energies()
{
    ke = 0;
    pe = 0;
    for (unsigned int m = 0; m < n_motors.size(); m++)
    {
        ke += n_motors[m]->get_kinetic_energy();
        pe += n_motors[m]->get_stretching_energy();
        //pe += n_motors[m]->get_stretching_energy_fene();
    }
}

void motor_ensemble::update_stress()	{
	unsigned int m, n;
	array<double, 4> m_stress;
	stress = {0, 0, 0, 0};
	for (m = 0; m < n_motors.size(); m++)	{
		m_stress = n_motors[m]->get_stress();
		for (n = 0; n < 4; n ++) {
			stress[n] += m_stress[n];
		}
	}
}
 
double motor_ensemble::get_potential_energy(){
    return pe;
}

 
void motor_ensemble::print_ensemble_thermo(){
    cout<<"\nAll Motors\t:\tKE = "<<ke<<"\tPEs = "<<pe<<"\tPEb = "<<0<<"\tTE = "<<(ke+pe);
	cout<<"\n\tAttached: "<<n_attached<<" / "<<n_motors.size();
	cout<<"\n\tsigma_xx = "<<stress[0]<<" sigma_xy = "<<stress[1];
	cout<<"\n\tsigma_yx = "<<stress[2]<<" sigma_yy = "<<stress[3];
}
