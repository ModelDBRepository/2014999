//% Joachim Behar, 04-05-2017
//% Bioenergetics and Bioelectric Systems Lab
//% Technion, Haifa, Israel
//%
//% When using this model, please cite:
//%   Behar, Joachim, et al. "The Autonomic Nervous System Regulates the Heart
//%   Rate through cAMP-PKA Dependent and Independent Coupled-Clock Pacemaker
//%   Cell Mechanisms." Frontiers in Physiology 7 (2016).
//%
//%
//%    This program is free software; you can redistribute it and/or modify
//%    it under the terms of the GNU General Public License as published by
//%    the Free Software Foundation; either version 2 of the License, or
//%    (at your option) any later version.
//%
//%    This program is distributed in the hope that it will be useful,
//%    but WITHOUT ANY WARRANTY; without even the implied warranty of
//%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//%    GNU General Public License for more details.
//%
//%    You should have received a copy of the GNU General Public License
//%    along with this program; if not, write to the Free Software
//%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//clear all; close all; clc;
//*******************************************************************************
//
// Linked to ICNS model
// Pei-Chi Yang @ Clancy Lab
// May 2020
///******************************************************************************

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <random>
#include <new>

// for mkdir
#include <sys/types.h>
#include <sys/stat.h>

const double pi = 3.141592653589793;

#include "load_constants.h"
#include "update_I.h"
#include "integrate_rk2.h"



using namespace std;

int main( int argc, char *argv[] ) {
    cout.precision(16);
    
    //%% constants
    //con = load_constants;
    cell_con con;
    load_constants( &con );
    
    
    //LINEWITH = 2;
    con.ALL_VAR = 1; // % outputting all the currents? (1:yes, 0:no) flag!!
    
    //% == no caged cAMP
    con.cAB_Conc = 0;
    con.cAB_on_rate = 0;
    con.cAB_off_rate = 0;
    
   
    //
    //%% solver
    double Finit[] = { 0.0001, 0.000223, 0.029, 1.35, 0.00005, 0.22, 0.69, 0.042, 0.089,
        0.032, 0.7499955, 3.4e-6, 1.1e-6, 0.25,
        -65., 0., 1., 1., 0., 0., 1., 0., 1., 0., 1., 1., 0., 0., 1.,
        19.73, 0.23, 0.06, 0.02, 0.06, 1.75e-6,
        0.0004, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.,
        0.042,0.042,0.00065,0.025,0.025,0.00064,0.072,0.0733,0.00065,0.002,0.002,0.0005,0.04,0.04,0.00065,0.0094,0.01,0.00062,0.1,5.0,1.1, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001};
    
    int sizeFinit;
    //if con.ALL_VAR % initial conditions
    if( con.ALL_VAR ) {
        sizeFinit = 73;//
        
        
    } else { // % initial conditions with only the coupled variables
        sizeFinit = 36;
       
    }
    
    
    double t = 0;
    double dt = 1./128;
    
    char *InputFile1 = argv[1];
    FILE * input1 = fopen( InputFile1, "r" );
    
    double t_ISO_s, t_ISO_e, t_ACh_s, t_ACh_e;
    char *InputFile2 = argv[2];
    FILE * input2 = fopen( InputFile2, "r" );
    
    fscanf( input2, "%le", &t_ISO_s );
    fscanf( input2, "%le", &t_ISO_e );
    fscanf( input2, "%le", &t_ACh_s );
    fscanf( input2, "%le", &t_ACh_e );
    
    fclose( input2 );

    t_ISO_s = t_ISO_s*1000;
    t_ISO_e = t_ISO_e*1000;
    t_ACh_s = t_ACh_s*1000;
    t_ACh_e = t_ACh_e*1000;
    
    FILE * output_CL = fopen("cycles_icnsach.txt", "w");
    
    FILE * output1 = fopen( "outputs_icnsach.txt", "w" );
    cout << t << "\t" << Finit[14] << endl;
    
    FILE * output2 = fopen( "alloutputs_icnsach.txt", "w" );
    
    fprintf( output2, "%16.14f\t", t );
    
    for( int i = 0; i < sizeFinit; i++ ) {
        fprintf( output2, "%16.14e\t", Finit[i] );
    }
    fprintf( output2, "\n" );
    
        
    double dvdt1 = 0;
    double dvdt2 = 0;
    double ddvdt21 = 0;
    double ddvdt22 = 0;
    double ct1 = 0;
    double ct2 = 0;
    double CL = 0;    //
    double t1 = -1;
    double t2 = -1;
   
    double ACh;
    double iso;
    double t_read;
    
    
    while( fscanf( input1, "%le", &t_read ) != EOF ) {
                
       
        
        con.ISO = 0; // % 100 [nM] modulates activation of Beta-adrenergic receptor (B-AR)
        con.CCh = 0; // % concentration of carbachol [nM]
        
        fscanf( input1, "%le", &(ACh) );
        fscanf( input1, "%le", &(iso) );

        if ( t >= t_ISO_s && t <= t_ISO_e ) {
            con.ISO = iso * 0.2; //nM concentration
        }
        
        
        if (t >= t_ACh_s && t <= t_ACh_e ) {
            con.CCh = ACh * 0.7; //nM concentration
        } else {
            con.CCh = ACh * 0.1;
        }
                
        
        
        double v1 = Finit[14];
        
        integrate_rk2( update_I, t, sizeFinit, Finit, dt , &con );
        
        t += dt;
        for( int id = 43; id < 61; id++ ) {
            if( Finit[id] < 0 ) {
                Finit[id] = 0;
            }
        }
        
        dvdt2 = (Finit[14] - v1) / dt;
        ddvdt22 = dvdt2 - dvdt1;
        if( dvdt1 > 5. && dvdt2 > 5. && ddvdt22 < 0. && ddvdt21 >= 0. ) {
            ct2 = t;
            
            CL = ct2 - ct1;
            ct1 = ct2;
            
            cout << "t = " << t << " \t ICNS = " << con.ISO << "\t Cycle Length = " << CL << "\t cch = " << con.CCh << endl;
            fprintf (output_CL, "%-16.4e\n", CL);
            
            
        }
        ddvdt21 = ddvdt22;
        dvdt1 = dvdt2;
        
       
        
        if( fmod(t,1) < dt && t > 0e3) {
            
            
            fprintf( output1, "%-16.14f\t%-16.14e\t%-16.14e\t%-16.14e\t%-16.14e\t%-16.14e\t%-16.14e\n", t, con.I_NCX, con.I_f, con.I_CaL, con.I_KACh, con.CCh, con.ISO );
            
            fprintf( output2, "%-16.14e\t", t );
            
            for( int i = 0; i < sizeFinit; i++ ) {
                fprintf( output2, "%-16.14e\t", Finit[i] );
            }
            fprintf( output2, "\n" );
            
            
        }
    }
    
    
    
    fclose( output1);
    fclose( output2);
    fclose( input1);
    
    
    return 0;
}
