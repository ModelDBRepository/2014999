/*
 *  soltis_biophysJ2010_masterCompute.h
 *
 *
 *  Created by Mao-Tsuen Jeng on 12/22/11.
 *  Copyright 2011 __Clancy Lab__. All rights reserved.
 *
 */

/*
 function yfinal = soltis_biophysJ2010_masterCompute
 % This function calls the ode solver and plots results.
 % Dependency: soltis_biophysJ2010_masterODEfile.m
 %
 % Re-implemented by Anthony Soltis <ars7h@virginia.edu> for CaMKII/PKA
 % regulation of ECC model
 %
 % Author: Jeff Saucerman <jsaucerman@virginia.edu>
 % Copyright 2008, University of Virginia, All Rights Reserved
 %
 % Reference: JJ Saucerman and DM Bers, Calmodulin mediates differential
 % sensitivity of CaMKII and calcineurin to local Ca2+ in cardiac myocytes.
 % Biophys J. 2008 Aug 8. [Epub ahead of print]
 % Please cite the above paper when using this model.
 */

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <time.h>
#include <string.h>
#include <random>
#include <new>


#include <sys/types.h>
#include <sys/stat.h>


typedef struct {
    double I_Ca_store, I_to_store[3], I_Na_store, I_K1_store, ibar_store, gates[2];
    double Jserca, IKs_store, Jleak[2], ICFTR, Incx, IKr_store;
    double JCaCyt, JCaSL, JCaDyad;
} pars_rec;

//Structure holding ODE and parameter arrays
typedef struct {
    double y0n[212];
    pars_rec pars1;
    double p[36];
    double v, v_old, ryr_old, I_Total;
    double dvdt, dcaidt, cai_old;
    double ui1_v, uj1_v, ui2_v, uj2_v, dun_v, tu_v; // For computing change of direction.
    double ui1_cai, uj1_cai, ui2_cai, uj2_cai, dun_cai, tu_cai; // For computing change of direction.
    
} Cell;

//Total tissue width and length
const int tw = 1;
const int tl = 1;


typedef struct {
    double t, tt;
    int tStep, counter, beat;
    Cell cellData[tw][tl];
} SimState;


SimState theState;
SimState *S = &theState;

const double pi = 3.141592653589793;

#include "integrate_rk2.h"
#include "soltis_biophysJ2010_masterODEfile.h"


int main( int argc, char *argv[] ) {
    
    
    double freq = 1;                  // [Hz] CHANGE DEPENDING ON FREQUENCY
    double cycleLength = 1.e3/freq;     // [ms]
    
    double Iso_CL;//cycleLength;
    double Ach_CL;//cycleLength;
    
    
    char name1[60];
    char name2[60];
    char name3[60];
    
    Cell * theCell;
    pars_rec * pars1;
    double * y0n;
    double * y0n2;
    S->tStep = 0;
    double * p;
    
    //// Parameters for external modules
    // ECC and CaM modules
    double CaMtotDyad = 418.;           // [uM]
    double BtotDyad = 1.54/8.293e-4;    // [uM]
    double CaMKIItotDyad = 120.;        // [uM]
    double CaNtotDyad = 3.e-3/8.293e-4; // [uM]
    double PP1totDyad = 96.5;           // [uM]
    double CaMtotSL = 5.65;             // [uM]
    double BtotSL = 24.2;               // [uM]
    double CaMKIItotSL = 120.*8.293e-4; // [uM]
    double CaNtotSL = 3.e-3;            // [uM]
    double PP1totSL = 0.57;             // [uM]
    double CaMtotCyt = 5.65;            // [uM]
    double BtotCyt = 24.2;              // [uM]
    double CaMKIItotCyt = 120.*8.293e-4;// [uM]
    double CaNtotCyt = 3.e-3;           // [uM]
    double PP1totCyt = 0.57;            // [uM]
    
    // ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
    char expression[] = "WT";
    double CKIIOE = 0; // Should be zero during 'WT' and 'KO' runs
    
    //if ( expression == "OE" ) {
    if( strcmp(expression, "OE" ) == 0 ) {
        int CKIIOE = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
        CaMKIItotDyad = 120*6;        // [uM]
        CaMKIItotSL = 120*8.293e-4*6; // [uM]
        CaMKIItotCyt = 120*8.293e-4*6;// [uM]
    }	else if ( strcmp( expression, "KO" ) == 0 ) {
        CaMKIItotDyad = 0;          // [uM]
        CaMKIItotSL = 0;            // [uM]
        CaMKIItotCyt = 0;           // [uM]
    }
    // end
    
    // For Recovery from inactivation of LCC
    double recoveryTime = 10;  // initialize to smallest value
    
    // Parameters for CaMKII module
    double LCCtotDyad = 31.4*.9;       // [uM] - Total Dyadic [LCC] - (umol/l dyad)
    double LCCtotSL = 0.0846;          // [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
    double RyRtot = 382.6;             // [uM] - Total RyR (in Dyad)
    double PP1_dyad = 95.7;            // [uM] - Total dyadic [PP1]
    double PP1_SL = 0.57;              // [uM] - Total Subsarcolemmal [PP1]
    double PP2A_dyad = 95.76;          // [uM] - Total dyadic PP2A
    double OA = 0;                     // [uM] - PP1/PP2A inhibitor Okadaic Acid
    double PLBtot = 38;                // [uM] - Total [PLB] in cytosolic units
    
    // Parameters for BAR module
    double Ligtot = 0;                 // [uM] - SET LIGAND CONCENTRATION HERE
    double LCCtotBA = 0.025;           // [uM] - [umol/L cytosol]
    double RyRtotBA = 0.135;           // [uM] - [umol/L cytosol]
    double PLBtotBA = 38;              // [uM] - [umol/L cytosol]
    double TnItotBA = 70;              // [uM] - [umol/L cytosol]
    double IKstotBA = 0.025;           // [uM] - [umol/L cytosol]
    double ICFTRtotBA = 0.025;         // [uM] - [umol/L cytosol]
    double PP1_PLBtot = 0.89;          // [uM] - [umol/L cytosol]
    double PLMtotBA = 48;
    
  
    
    char stateFileName[] = "finalStates.txt"; 

    
    char *InputFile2 = argv[2];
    FILE * input2 = fopen( InputFile2, "r" );
    
    char *InputFile1 = argv[1];
    FILE * input3 = fopen( InputFile1, "r" );
    
    char *InputFile3 = argv[3];
    FILE *input_time = fopen( InputFile3, "r" );
    double simTime;
    fscanf( input_time, "%lf", &simTime );
    fclose( input_time );
    
    double t_ISO_s, t_ISO_e, t_ACh_s, t_ACh_e;
    char *InputFile4 = argv[4];
    FILE * input4 = fopen( InputFile4, "r" );
    
    fscanf( input4, "%le", &t_ISO_s );
    fscanf( input4, "%le", &t_ISO_e );
    fscanf( input4, "%le", &t_ACh_s );
    fscanf( input4, "%le", &t_ACh_e );
    
    fclose( input4 );

    t_ISO_s = t_ISO_s*1000;
    t_ISO_e = t_ISO_e*1000;
    t_ACh_s = t_ACh_s*1000;
    t_ACh_e = t_ACh_e*1000;

    
    int sy0 = 212;
    double dtbase = 1. / 128;
    double dt = dtbase;
    int fold = 1 / dt;
    double t = 0;
    double tb = 0;
     // simulation time [ms]/ dt;
    int iter_max = simTime*1000/dt;
    int tcout = 0;
    time_t t_start, t_end;
    double dif, dv, CL;
    int ww, ll;
    double I_inj;
    double t0 = 0;
    double tt1;
    double t1 = -1;
    double t2 = -1;
    double ACh, cardiac_s;
    
    // Read initial data
    theCell = &(S->cellData[0][0]);
    y0n = theCell->y0n;
    
    
    FILE *fp = fopen( stateFileName, "r");
    
    for ( int idy = 0; idy < sy0; idy++ ) {
        fscanf ( fp, "%lf", &y0n[idy] );
        
    }
    fclose( fp );
    
    theCell->tu_v = 0;
    theCell->tu_cai = 0;
    
    ww = 0;
    ll = 0;
    
    S->cellData[ww][ll].p[0] = cycleLength;
    S->cellData[ww][ll].p[1] = recoveryTime;
    S->cellData[ww][ll].p[2] = CaMtotDyad;
    S->cellData[ww][ll].p[3] = BtotDyad;
    S->cellData[ww][ll].p[4] = CaMKIItotDyad;
    S->cellData[ww][ll].p[5] = CaNtotDyad;
    S->cellData[ww][ll].p[6] = PP1totDyad;
    S->cellData[ww][ll].p[7] = CaMtotSL;
    S->cellData[ww][ll].p[8] = BtotSL;
    S->cellData[ww][ll].p[9] = CaMKIItotSL;
    S->cellData[ww][ll].p[10] = CaNtotSL;
    S->cellData[ww][ll].p[11] = PP1totSL;
    S->cellData[ww][ll].p[12] = CaMtotCyt;
    S->cellData[ww][ll].p[13] = BtotCyt;
    S->cellData[ww][ll].p[14] = CaMKIItotCyt;
    S->cellData[ww][ll].p[15] = CaNtotCyt;
    S->cellData[ww][ll].p[16] = PP1totCyt;
    S->cellData[ww][ll].p[17] = LCCtotDyad;
    S->cellData[ww][ll].p[18] = RyRtot;
    S->cellData[ww][ll].p[19] = PP1_dyad;
    S->cellData[ww][ll].p[20] = PP2A_dyad;
    S->cellData[ww][ll].p[21] = OA;
    S->cellData[ww][ll].p[22] = PLBtot;
    S->cellData[ww][ll].p[23] = LCCtotSL;
    S->cellData[ww][ll].p[24] = PP1_SL;
    S->cellData[ww][ll].p[25] = 0.0;//iso_total;
    S->cellData[ww][ll].p[26] = LCCtotBA;
    S->cellData[ww][ll].p[27] = RyRtotBA;
    S->cellData[ww][ll].p[28] = PLBtotBA;
    S->cellData[ww][ll].p[29] = TnItotBA;
    S->cellData[ww][ll].p[30] = IKstotBA;
    S->cellData[ww][ll].p[31] = ICFTRtotBA;
    S->cellData[ww][ll].p[32] = PP1_PLBtot;
    S->cellData[ww][ll].p[33] = CKIIOE;
    S->cellData[ww][ll].p[34] = PLMtotBA;
    S->cellData[ww][ll].p[35] = 0.0;//ach_total;
    
    sprintf(name1, "ap_icnsach.txt");
    FILE * fpy = fopen( name1, "w" );
    sprintf(name2, "ca_icnsach.txt");
    FILE * fpall = fopen( name2,  "w" );
    
    
    int syid1 = 1;
    int yid1[] = { 39};
    int yid;
    
    time_t startTime;
    time_t previousTime;
    
    const time_t timeSave = 1*60*60;
    startTime = time(NULL);
    previousTime = startTime;
    
    tb = 0;
    double t_read;
    
    for ( int iter = 1; iter <= iter_max; iter ++ ) {
        tcout++;
        
        t = iter * dt;  // total time
        
        t0 = 0;
        
        
        if( tb == 0 ) {
            fscanf( input3, "%le", &CL );
        }
        
        ww = 0;
        ll = 0;
        

        fscanf( input2, "%le", &(t_read) );
        fscanf( input2, "%le", &(ACh));
        fscanf( input2, "%le", &(cardiac_s) );
        
       
            
        S->cellData[ww][ll].p[25] = 0;
            
        if ( S->tt >= t_ISO_s && S->tt <= t_ISO_e ) {
            
            S->cellData[ww][ll].p[25] = cardiac_s * 0.2/1e3;// to uM
            
        }
       

        if ( S->tt >= t_ACh_s && S->tt <= t_ACh_e ) {
       
                S->cellData[ww][ll].p[35] = ACh * 0.7/1e3; //to uM

        } else {
            
            S->cellData[ww][ll].p[35] = ACh * 0.1/1e3; //to uM
        }
        
        
        
        
        if ( tb < 5 && S->tt >= (0e3) && S->tt < 200e3 ) {
            I_inj = -9.5;  //  injection/stimulus current
        }
        else {
            I_inj=0.0;
        }
        
        theCell = &(S->cellData[ww][ll]);
        pars1 = &(theCell->pars1);
        y0n = theCell->y0n;
        p = theCell->p;
        
        theCell->v_old = y0n[38];
        theCell->ryr_old = y0n[14];
        theCell->cai_old = y0n[37];
        
        
        integrate_rk2( soltis_biophysJ2010_masterODEfile, &t, sy0, y0n, dt, p, pars1 ) ;
        theCell->I_Total = (theCell->v_old - y0n[38]) / dt + I_inj;
        theCell->v = y0n[38];
        
        
      
        
        theCell = &(S->cellData[ww][ll]);
        
        dv=dt*(-theCell->I_Total);
        
        theCell->y0n[38] = theCell->v_old + dv;
        
        S->t = tb + dt;
        S->tt = t;
        tb += dt;
        if( tb >= CL ) {
            tb = 0;
        }
        
        y0n = theCell->y0n;
        double voltage1=y0n[38];
        double cai1=y0n[37];
        double srca1=y0n[30];
        double cAMP_cav = y0n[209];
        double cAMP_ecav = y0n[210];
        double cAMP_cyt = y0n[211];
        
        if ( (tcout % 12800 ) == 0 ) {
            
            cout << "Simulation Time: " << t << "\tTime spent: " << (time(NULL) - startTime)/60 << " min "<< (time(NULL) - startTime)%60 << " sec" << endl;
            
            cout << "ICNS: " << theCell->p[25] << endl;
            cout << "ACh: " << theCell->p[35] << endl;
            cout << "CL: " << CL << endl;
           
            
        }
        
        
        double N2, dui, duj;
        
        
        theCell = &(S->cellData[ww][ll]);
        y0n = theCell->y0n;
        double voltage = y0n[38];
        double cai = y0n[37];
        theCell->dvdt = ( voltage - theCell->v_old ) / dt;
        theCell->dcaidt = ( cai - theCell->cai_old ) / dt;
        
        // Computing direction change
        N2 = pow( ( 1. + theCell->dvdt * theCell->dvdt ), 0.5 );
        theCell->ui2_v = 1. / N2;
        theCell->uj2_v = theCell->dvdt / N2;
        dui = theCell->ui1_v - theCell->ui2_v;
        duj = theCell->uj1_v - theCell->uj2_v;
        theCell->dun_v += pow( ( dui * dui + duj * duj ), 0.5 );
        theCell->ui1_v = theCell->ui2_v;
        theCell->uj1_v = theCell->uj2_v;
        
        theCell->dcaidt = 1.e5 * theCell->dcaidt; // To avoid truncation error in N2.
        N2 = pow( ( 1. + theCell->dcaidt * theCell->dcaidt ), 0.5 );
        theCell->ui2_cai = 1. / N2;
        theCell->uj2_cai = theCell->dcaidt / N2;
        dui = theCell->ui1_cai - theCell->ui2_cai;
        duj = theCell->uj1_cai - theCell->uj2_cai;
        theCell->dun_cai += pow( ( dui * dui + duj * duj ), 0.5 );
        theCell->ui1_cai = theCell->ui2_cai;
        theCell->uj1_cai = theCell->uj2_cai;
        
        // Save if there is an significant change of direction or end of a cycle.
        if(  ( ( ( exp( 5. * ( S->tt - theCell->tu_v ) ) * theCell->dun_v ) > 0.05 ) || ( fmod( S->tt, CL ) < dt ) ) ) {
            
            fprintf(fpy, "%-12.6f\t%-12.10e\n", t , voltage);
            theCell->dun_v = 0.;
            theCell->tu_v = S->tt;
            
        }
        
        
        if(  ( ( ( exp( 100. * ( S->tt - theCell->tu_cai ) ) * theCell->dun_cai ) > 0.05 ) || ( fmod( S->tt, CL ) < dt ) ) ) {
            
            fprintf(fpall, "%-12.6f\t%-12.10e\n", t , cai );
            theCell->dun_cai = 0.;
            theCell->tu_cai = S->tt;
            
        }
        
    }
    
    
    fclose (fpy);
    fclose( fpall );
   
    fclose (input2);
    
   
    
    return 0;
}
