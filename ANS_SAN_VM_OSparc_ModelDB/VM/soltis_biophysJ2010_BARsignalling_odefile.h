/*
 *  soltis_biophsJ2010_BARsignalling_odefile.h
 *
 *
 *  Created by Mao-Tsuen Jeng on 12/22/11.
 *  & Pei-Chi Yang on 01/03/12
 *  Copyright 2011 __UCDavis__. All rights reserved.
 *
 */

// ************************************************************
//    Modified from the following matlab function:
//       DAEs are solved using iterative method.
//       ODEs are integrated using RK2 with adaptive step size.
// ************************************************************ 

/*
 % function ydot = soltis_biophysJ2010_BARsignalling_odefile(t,y,pin)
 % Re-implemented by Anthony Soltis <ars7h@virginia.edu> for communication
 % with Shannon (2004) EC coupling model. Final version 07/21/10
 
 % saucerman_circres2004.m
 % coupled beta adrenergic signaling/EC for adult rabbit ventricular
 % myocytes, with extensions for yotiao interactions, KCNQ1-G589D mutation
 %
 % Copyright (2004) The Regents of the University of California
 % All Rights Reserved
 % Permission to use, copy, and modify, any part of this software for
 % academic purposes only, including non-profit  education and research,
 % without fee, and without a written agreement is hereby granted, provided
 % that the above copyright notice, this paragraph and the following three
 % paragraphs appear in all copies.
 % The receiving party agrees not to further distribute nor disclose the
 % source code to third parties without written permission and not to use
 % the software for research under commercial sponsorship or for any other
 % any commercial undertaking.  Those desiring to incorporate this software
 % into commercial products of to use it for commercial purposes should
 % contact the Technology Transfer Office, University of California, San
 % Diego, 9500 Gilman Drive, La Jolla, California, 92093-0910, Ph: (619)
 % 534 5815, referring to Case SDC98008.
 % IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 % FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 % INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 % THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
 % DAMAGE.
 % THE SOFTWARE PROVIDED HEREUNDER IS ON AN AS IS BASIS, AND THE
 % UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
 % SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
 % CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY
 % KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
 % IMPLIED WARRANTIES OF MECHANT ABILITY OR FITNESS FOR A PARTICULAR
 % PURPOSE, OR THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,
 % TRADEMARK OR OTHER RIGHTS.
 %
 % Those using this code should acknowledge Dr. Andrew D. McCulloch and the
 % National Biomedical Computation Resource (NBCR), NIH grant #P41 RR08605.
 %
 % References:
 % -----------
 % Jeffrey J. Saucerman, Sarah N. Healy, Mary E. Belik, Jose L. Puglisi, and
 % Andrew D. McCulloch.  "Proarrhythmic Consequences of a KCNQ1 AKAP-Binding
 % Domain Mutation: Computational Models of Whole Cells and Heterogeneous
 % Tissue", Circ. Res., Vol 95: 1216-1224 (2004).
 %
 % Jeffrey J. Saucerman and Andrew D. McCulloch
 % "Mechanistic systems models of cell signaling networks: a case study of
 % myocyte adrenergic regulation", Prog. Biophys. Mol. Bio., Vol 84: 261-278 (2004).
 %
 % Jeffrey J. Saucerman, Laurence L. Brunton, Anushka P. Michailova, and Andrew D. McCulloch
 % "Modeling beta-adrenergic control of cardiac myocyte contractility in
 % silico", J. Biol. Chem., Vol 278: 47977-48003 (2003).
 %
 % Last modified: 12/20/2004
 % Implemented by: Jeffrey Saucerman <jsaucer@ucsd.edu>
 %
 % Notes:
 % - This code was used for the single cell simulations in the Circ Res ms.
 % - In Matlab 7, you can enable cell mode to jump between modules.
 % - Please email me for any questions or comments.
 */
// ************************************************************************
// ************************************************************************
#ifndef soltis_biophysJ2010_BARsignalling_odefile_H
#define soltis_biophysJ2010_BARsignalling_odefile_H

void soltis_biophysJ2010_BARsignalling_odefile( double tt, double *y, double *pin, double *ydot );

void soltis_biophysJ2010_BARsignalling_odefile( double tt, double *y, double *pin, double *ydot ) {
    
    //// Assign passed in params
    double L_iso = pin[0];
    double LCCtot = pin[1];
    double RyRtot = pin[2];
    double PLBtot = pin[3];
    double TnItot = pin[4];
    double IKstot = pin[5];
    double ICFTRtot = pin[6];
    double PP1_PLBtot = pin[7];
    double PLMtot= pin[8];
    double L_ach = pin[9];
    //cout << L_ach << endl;
    //// Parameters
    //// ----- Signaling model parameters -------
    // b-AR/Gs module
    double p[99];
    p[0] = L_iso;    // Ltotmax   [uM] ** apply agonist concentration here **
    p[1] = 0.028;   // sumb1AR   [uM]
    p[2] = 3.83;    // Gstot     [uM]
    p[3] = 0.285;   // Kl        [uM]
    p[4] = 0.062;   // Kr        [uM]
    p[5] = 33.0;    // Kc        [uM]
    p[6] = 1.1e-3;  // k_barkp   [1/sec]
    p[7] = 2.2e-3;  // k_barkm   [1/sec]
    p[8] = 3.6e-3;  // k_pkap    [1/sec/uM]
    p[9] = 2.2e-3; // k_pkam    [1/sec]
    p[10] = 16.0;   // k_gact    [1/sec]
    p[11] = 0.8;    // k_hyd     [1/sec]
    p[12] = 1.21e3; // k_reassoc [1/sec/uM]
    // cAMP module
    p[13] = 0.047;  // AC_tot    [uM]
    p[14] = 5.0e3;  // ATP       [uM]
    p[15] = 0.036;  // PDE3tot   [uM]   // Changed from .06 to .036
    p[16] = 0.036;  // PDE4tot   [uM]
    p[17] = 0.0;    // IBMXtot   [uM]
    p[18] = 0.0;    // Fsktot    [uM] (10 uM when used)
    p[19] = 0.2;    // k_ac_basal[1/sec]
    p[20] = 8.5;    // k_ac_gsa  [1/sec]
    p[21] = 7.3;    // k_ac_fsk  [1/sec]
    p[22] = 1.03e3; // Km_basal  [uM]
    p[23] = 315.0;  // Km_gsa    [uM]
    p[24] = 860.0;  // Km_fsk    [uM]
    p[25] = 0.4;    // Kgsa      [uM]
    p[26] = 44.0;   // Kfsk      [uM]
    p[27] = 3.5;    // k_pde3    [1/sec]
    p[28] = 0.15;   // Km_pde3   [uM]
    p[29] = 5.0;    // k_pde4    [1/sec]
    p[30] = 1.3;    // Km_pde4   [uM]
    p[31] = 30.0;   // Ki_ibmx   [uM]
    // PKA module
    p[32] = 0.46;   // PKAItot   [uM]
    p[33] = 0.084;  // PKAIItot  [uM]
    p[34] = 0.18;   // PKItot    [uM]
    p[35] = 9.14;   // Ka        [uM]
    p[36] = 1.64;   // Kb        [uM]
    p[37] = 4.375;  // Kd        [uM]
    p[38] = 0.2e-3; // Ki_pki    [uM]
    // PLB module
    p[39] = 10;     // epsilon   [none]
    p[40] = PLBtot; // PLBtot    [uM]
    p[41] = PP1_PLBtot;   // PP1tot    [uM]
    p[42] = 0.3;    // Inhib1tot [uM]
    p[43] = 54;     // k_pka_plb     [1/sec]
    p[44] = 21;     // Km_pka_plb    [uM]
    p[45] = 8.5;    // k_pp1_plb     [1/sec]
    p[46] = 7.0;    // Km_pp1_plb    [uM]
    p[47] = 60;     // k_pka_i1      [1/sec]
    p[48] = 1.0;    // Km_pka_i1     [uM]
    p[49] = 14.0;   // Vmax_pp2a_i1  [uM/sec]
    p[50] = 1.0;    // Km_pp2a_i1    [uM]
    p[51] = 1.0e-3; // Ki_inhib1     [uM]
    // LCC module
    p[52] = LCCtot; // LCCtot        [uM]
    p[53] = 0.025;  // PKAIIlcctot   [uM]
    p[54] = 0.025;  // PP1lcctot     [uM]
    p[55] = 0.025;  // PP2Alcctot    [uM]
    p[56] = 54;     // k_pka_lcc     [1/sec]
    p[57] = 21;     // Km_pka_lcc    [uM]
    p[58] = 8.52;   // k_pp1_lcc     [1/sec]
    p[59] = 3;      // Km_pp1_lcc    [uM]
    p[60] = 10.1;   // k_pp2a_lcc    [1/sec]
    p[61] = 3;      // Km_pp2a_lcc   [uM]
    // RyR module
    p[62] = RyRtot; // RyRtot        [uM]
    p[63] = 0.034;  // PKAIIryrtot   [uM]
    p[64] = 0.034;  // PP1ryr        [uM]
    p[65] = 0.034;  // PP2Aryr       [uM]
    p[66] = 54;     // kcat_pka_ryr  [1/sec]
    p[67] = 21;     // Km_pka_ryr    [uM]
    p[68] = 8.52;   // kcat_pp1_ryr  [1/sec]
    p[69] = 7;      // Km_pp1_ryr    [uM]
    p[70] = 10.1;   // kcat_pp2a_ryr [1/sec]
    p[71] = 4.1;    // Km_pp2a_ryr   [uM]
    // TnI module
    p[72] = TnItot; // TnItot        [uM]
    p[73] = 0.67;   // PP2Atni       [uM]
    p[74] = 54;     // kcat_pka_tni  [1/sec]
    p[75] = 21;     // Km_pka_tni    [uM]
    p[76] = 10.1;   // kcat_pp2a_tni [1/sec]
    p[77] = 4.1;    // Km_pp2a_tni   [uM]
    // Iks module
    p[78] = IKstot; // Iks_tot       [uM]
    p[79] = 0.025;  // Yotiao_tot    [uM]
    p[80] = 0.1e-3; // K_yotiao      [uM] ** apply G589D mutation here **
    p[81] = 0.025;  // PKAII_ikstot  [uM]
    p[82] = 0.025;  // PP1_ikstot    [uM]
    p[83] = 54;     // k_pka_iks     [1/sec]
    p[84] = 21;     // Km_pka_iks    [uM]
    p[85] = 8.52;   // k_pp1_iks     [1/sec]
    p[86] = 7;      // Km_pp1_iks    [uM]
    // Icftr Module - Added 04/30/10 by Anthony Soltis
    p[87] = ICFTRtot;  // CFTR_tot      [uM]
    p[88] = 0.025;  // PKAII_CFTRtot [uM]
    p[89] = 0.025;  // PP1_CFTRtot   [uM]
    p[90] = 54;     // k_pka_CFTR    [1/sec]
    p[91] = 8.5;    // Km_pka_CFTR   [uM]
    p[92] = 8.52;   // k_pp1_CFTR    [1/sec]
    p[93] = 7;      // Km_pp1_CFTR   [uM]
    
    // PLM Module - Added 09/13/12 by Ele Grandi
    p[94] = PLMtot; // PLMtot    [uM]
    p[95] = 54;     // k_pka_plm     [1/sec]
    p[96] = 21;     // Km_pka_plm    [uM]
    p[97] = 8.5;    // k_pp1_plm     [1/sec]
    p[98] = 7.0;    // Km_pp1_plm    [uM]
    
    
    //
    double myeps = 1.e-14;
    //
    int iter = 0;
    
    
    const double volume = 38e-6; //cell volume (L=100um;R=10um) in uL
    const double myop_volume = volume*0.5;
    const double cav_volume = myop_volume*0.01;
    const double ecav_volume = myop_volume*0.02;
    const double myop_membr_volume = myop_volume*0.07;
    
    const double flux1 = 0.9e-14;  //Ecav-Bulk
    const double flux2 = 7.5e-14; //Cav-Bulk
    const double flux3 = 5e-15; //Cav-Ecav
    
    //Rb1_cav
    const double K_H = 0.062; //microM
    const double K_L = 0.567; //microM
    const double K_C = 8.809; //microM
    const double kact1_cav = 0.1; // 1/s
    const double kact2_cav = 5; // 1/s
    const double khydr_cav = 0.8; // 1/s
    const double kreas_cav = 1.21E03 ; // (1/s)*(1/uM)
    
   
    double Rbeta1_tot_Cav;
    
    double RGs_Cav, LRGs_Cav, L_iso_Rbeta1_Cav, Gs_free_Cav;
    
    double beta1_rec_nr_cav = 0.75e5;
    double total_beta1_rec_cav; //concentration Beta 1 receptor - uM
    total_beta1_rec_cav = (beta1_rec_nr_cav/(6.23*cav_volume))*1e-11;
    double total_act_beta1_rec_cav = total_beta1_rec_cav;
    double Rbeta1_free_Cav = total_act_beta1_rec_cav; //uM
    
    double total_gs_nr_cav = 23.5e6;
    double total_gs_cav = 10; //total concentration of Gs - uM
    double total_avail_gs_cav = total_gs_cav;
    
    Gs_free_Cav = total_avail_gs_cav - (y[0] + y[2]);
    
    //R_beta
    RGs_Cav = (Rbeta1_free_Cav * Gs_free_Cav )/(Gs_free_Cav + K_C);
    
    L_iso_Rbeta1_Cav = (L_iso * (Rbeta1_free_Cav - RGs_Cav))/(L_iso + K_L) ;
    
    LRGs_Cav = (L_iso_Rbeta1_Cav*(Gs_free_Cav - RGs_Cav))/((Gs_free_Cav - RGs_Cav)+(K_C*K_H/K_L)) + (L_iso * RGs_Cav)/(L_iso + K_H);
    
    
    Rbeta1_tot_Cav = Rbeta1_free_Cav + L_iso_Rbeta1_Cav + LRGs_Cav + RGs_Cav;
    
    
    //Gs
    ydot[0] = (LRGs_Cav * kact2_cav + RGs_Cav * kact1_cav - y[0] * khydr_cav);
    
    ydot[1] = (LRGs_Cav * kact2_cav + RGs_Cav *kact1_cav - y[2] * y[1] * kreas_cav);
    
    ydot[2] = (y[0] * khydr_cav - y[2] * y[1] * kreas_cav);
    
    
    if (y[0]<0){
        y[0] = 0;
    }
    
    if (y[1]<0){
        y[1] = 0;
    }
    
    if (y[2]<0){
        y[2] = 0;
    }
    
    const double K_H3 = 0.16; //microM
    const double K_L3 = 11; //microM
    const double K_C3 = 30; //microM
    
    const double kact1_cav2 = 0.1; //1/s
    const double kact2_cav2 = 5;
    const double khydr_cav2 = 0.8;
    const double kreas_cav2 = 1.21E03 ; // (1/s)*(1/uM)
    
    //double dGia_GTP, dGi_beta_gamma, dGia_GDP;
    
    double RM2_tot_Cav;
    double RM2_free_Cav; //uM
    double Gi_free_Cav, RGi_Cav, LRGi_Cav, L_ach_RM2_Cav;
    
    double m2_rec_nr_cav = 0.6e5;
    double total_m2_rec_cav;
    double total_act_m2_rec_cav;
    
    total_m2_rec_cav = (m2_rec_nr_cav/(6.23*cav_volume))*1e-11;
    total_act_m2_rec_cav = total_m2_rec_cav;
    RM2_free_Cav = total_act_m2_rec_cav;
    
    double total_gi_nr_cav = 90e6;
    double total_gi_cav; //total concentration of gi - uM
    double total_avail_gi_cav;
    total_gi_cav = 20;
    total_avail_gi_cav = total_gi_cav;
    
    
    
    Gi_free_Cav = total_avail_gi_cav - (y[3] + y[5]);
    //R_M2
    
    RGi_Cav = (RM2_free_Cav * Gi_free_Cav )/(Gi_free_Cav+K_C3);
    
    L_ach_RM2_Cav = (L_ach * (RM2_free_Cav-RGi_Cav))/(L_ach+K_L3) ;
    
    LRGi_Cav = (L_ach_RM2_Cav*(Gi_free_Cav-RGi_Cav))/((Gi_free_Cav-RGi_Cav)+(K_C3*K_H3/K_L3)) + ((L_ach * RGi_Cav)/(L_ach + K_H3));
    
    RM2_tot_Cav = RM2_free_Cav + L_ach_RM2_Cav + LRGi_Cav + RGi_Cav;
    
    
    //Gi
    ydot[3] = (LRGi_Cav * kact2_cav2 + RGi_Cav * kact1_cav2 - y[3] * khydr_cav2);
    
    ydot[4] = (LRGi_Cav * kact2_cav2 + RGi_Cav *kact1_cav2 - y[5] * y[4] * kreas_cav2);
    
    ydot[5] = (y[3] * khydr_cav2 - y[5] * y[4] * kreas_cav2);
    
    if (y[3]<0){
        y[3] = 0;
    }
    
    if (y[4]<0){
        y[4] = 0;
    }
    
    if (y[5]<0){
        y[5] = 0;
    }
    
    
    
    
    //Rb1_ecav
    const double K_H4 = 0.062; //microM
    const double K_L4 = 0.567; //microM
    const double K_C4 = 8.809; //microM
    const double kact1_ecav = 0.1; // 1/s
    const double kact2_ecav = 5; // 1/s
    const double khydr_ecav = 0.8; // 1/s
    const double kreas_ecav = 1.21E03 ; // (1/s)*(1/uM)
    //double dGsa_GTP, dGs_beta_gamma, dGsa_GDP;
    double Rbeta1_tot_Ecav;
    double beta1_rec_nr_ecav = 2.6e5;
    double total_beta1_rec_ecav; //concentration Beta 1 receptor - uM
    total_beta1_rec_ecav = (beta1_rec_nr_ecav/(6.23*ecav_volume))*1e-11;
    double total_act_beta1_rec_ecav = total_beta1_rec_ecav;
    double Rbeta1_free_Ecav = total_act_beta1_rec_ecav;
    
    double total_gs_nr_ecav = 23.5e6;
    double total_gs_ecav; //total concentration of Gs - uM
    total_gs_ecav = 10;
    double total_avail_gs_ecav = total_gs_ecav;
    
    
    double Gs_free_Ecav, RGs_Ecav, LRGs_Ecav, L_iso_Rbeta1_Ecav;
    
    Gs_free_Ecav = total_avail_gs_ecav - (y[6] + y[8]);
    //R_beta
    
    RGs_Ecav = (Rbeta1_free_Ecav * Gs_free_Ecav )/(Gs_free_Ecav + K_C4);
    
    L_iso_Rbeta1_Ecav = (L_iso * (Rbeta1_free_Ecav - RGs_Ecav))/(L_iso + K_L4) ;
    
    LRGs_Ecav = (L_iso_Rbeta1_Ecav*(Gs_free_Ecav - RGs_Ecav))/((Gs_free_Ecav - RGs_Ecav)+(K_C4*K_H4/K_L4)) + (L_iso * RGs_Ecav)/(L_iso + K_H4);
    
    
    Rbeta1_tot_Ecav = Rbeta1_free_Ecav + L_iso_Rbeta1_Ecav + LRGs_Ecav + RGs_Ecav;
    
    
    //Gs
    ydot[6] = (LRGs_Ecav * kact2_ecav + RGs_Ecav * kact1_ecav - y[6] * khydr_ecav);
    
    ydot[7] = (LRGs_Ecav * kact2_ecav + RGs_Ecav *kact1_ecav - y[8] * y[7] * kreas_ecav);
    
    ydot[8] = (y[6] * khydr_ecav - y[8] * y[7] * kreas_ecav);
    
    
    if (y[6]<0){
        y[6] = 0;
    }
    
    if (y[7]<0){
        y[7] = 0;
    }
    
    if (y[8]<0){
        y[8] = 0;
    }
    
    
    
    const double K_H5 = 0.16; //microM
    const double K_L5 = 11; //microM
    const double K_C5 = 30; //microM
    
    const double kact1_ecav2 = 0.1;
    const double kact2_ecav2 = 5;
    const double khydr_ecav2 = 0.8;
    const double kreas_ecav2 = 1.21E03 ; // (1/s)*(1/uM)
    
    //double dGia_GTP, dGi_beta_gamma, dGia_GDP;
    
    double RM2_tot_Ecav;
    double RM2_free_Ecav; // uM
    double Gi_free_Ecav, RGi_Ecav, LRGi_Ecav, L_ach_RM2_Ecav;
    
    double m2_rec_nr_ecav = 1.2e5;
    double total_m2_rec_ecav;
    total_m2_rec_ecav = (m2_rec_nr_ecav/(6.23*ecav_volume))*1e-11;
    double total_act_m2_rec_ecav = total_m2_rec_ecav;
    RM2_free_Ecav = total_act_m2_rec_ecav;
    
    double total_gi_nr_ecav = 90e6;
    
    double total_gi_ecav; //total concentration of gi - uM
    double total_avail_gi_ecav;
    total_gi_ecav = 1;
    total_avail_gi_ecav = total_gi_ecav;
    
    
    
    Gi_free_Ecav = total_avail_gi_ecav - (y[9] + y[11]);
    //R_M2
    
    RGi_Ecav = (RM2_free_Ecav * Gi_free_Ecav )/(Gi_free_Ecav+K_C5);
    
    L_ach_RM2_Ecav = (L_ach * (RM2_free_Ecav-RGi_Ecav))/(L_ach+K_L5) ;
    
    LRGi_Ecav = (L_ach_RM2_Ecav*(Gi_free_Ecav-RGi_Ecav))/((Gi_free_Ecav-RGi_Ecav)+(K_C5*K_H5/K_L5)) + ((L_ach * RGi_Ecav)/(L_ach + K_H5));
    
    RM2_tot_Ecav = RM2_free_Ecav + L_ach_RM2_Ecav + LRGi_Ecav + RGi_Ecav;
    
    
    //Gi
    ydot[9] = (LRGi_Ecav * kact2_ecav2 + RGi_Ecav * kact1_ecav2 - y[9] * khydr_ecav2);
    
    ydot[10] = (LRGi_Ecav * kact2_ecav2 + RGi_Ecav *kact1_ecav2 - y[11] * y[10] * kreas_ecav2);
    
    ydot[11] = (y[9] * khydr_ecav2 - y[11] * y[10] * kreas_ecav2 );
    
    
    if (y[9]<0){
        y[9] = 0;
    }
    
    if (y[10]<0){
        y[10] = 0;
    }
    
    if (y[11]<0){
        y[11] = 0;
    }
    
    
    //Rb1_cyt
    const double K_H6 = 0.062; //microM
    const double K_L6 = 0.567; //microM
    const double K_C6 = 8.809; //microM
    const double kact1_cyt = 0.1; // 1/s
    const double kact2_cyt = 5; // 1/s
    const double khydr_cyt = 0.8; // 1/s
    const double kreas_cyt = 1.21E03 ; // (1/s)*(1/uM)
    
    //double dGsa_GTP, dGs_beta_gamma, dGsa_GDP;
    double Rbeta1_tot_Cyt;
    
    double RGs_Cyt, LRGs_Cyt, L_iso_Rbeta1_Cyt, Gs_free_Cyt;
    
    double beta1_rec_nr_cyt = 5e5;
    double total_beta1_rec_cyt = (beta1_rec_nr_cyt/(6.23*myop_membr_volume))*1e-11;
    double total_act_beta1_rec_cyt = total_beta1_rec_cyt;
    double Rbeta1_free_Cyt = total_act_beta1_rec_cyt;
    
    double total_gs_nr_cyt = 47e6;
    double total_gs_cyt = 10;
    double total_avail_gs_cyt = total_gs_cyt;
    
    Gs_free_Cyt = total_avail_gs_cyt - (y[12] + y[14]);
    
    //R_beta
    RGs_Cyt = (Rbeta1_free_Cyt * Gs_free_Cyt )/(Gs_free_Cyt + K_C6);
    
    L_iso_Rbeta1_Cyt = (L_iso * (Rbeta1_free_Cyt - RGs_Cyt))/(L_iso + K_L6) ;
    
    LRGs_Cyt = (L_iso_Rbeta1_Cyt*(Gs_free_Cyt - RGs_Cyt))/((Gs_free_Cyt - RGs_Cyt)+(K_C6*K_H6/K_L6)) + (L_iso * RGs_Cyt)/(L_iso + K_H6);
    
    
    Rbeta1_tot_Cyt = Rbeta1_free_Cyt + L_iso_Rbeta1_Cyt + LRGs_Cyt + RGs_Cyt;
    
    
    //Gs
    ydot[12] = (LRGs_Cyt * kact2_cyt + RGs_Cyt * kact1_cyt - y[12] * khydr_cyt );
    
    ydot[13] = (LRGs_Cyt * kact2_cyt + RGs_Cyt *kact1_cyt - y[14] * y[13] * kreas_cyt );
    
    ydot[14] = (y[12] * khydr_cyt - y[14] * y[13] * kreas_cyt );
    
    
    if (y[12]<0){
        y[12] = 0;
    }
    
    if (y[13]<0){
        y[13] = 0;
    }
    
    if (y[14]<0){
        y[14] = 0;
    }
    
    
    const double K_H2 = 0.16; //microM
    const double K_L2 = 11; //microM
    const double K_C2 = 30; //microM
    
    const double kact1_cyt2 = 0.1; //1/s
    const double kact2_cyt2 = 5;
    const double khydr_cyt2 = 0.8;
    const double kreas_cyt2 = 1.21E03 ; // (1/s)*(1/uM)
    
    //double dGia_GTP, dGi_beta_gamma, dGia_GDP;
    
    double RM2_tot_Cyt;
    double Gi_tot_Cyt; //uM
    double RM2_free_Cyt; //uM
    double Gi_free_Cyt, RGi_Cyt, LRGi_Cyt, L_ach_RM2_Cyt;
    
    double m2_rec_nr_cyt = 2.5e5;
    double total_m2_rec_cyt;
    total_m2_rec_cyt = (m2_rec_nr_cyt/(6.23*myop_membr_volume))*1e-11;
    double total_act_m2_rec_cyt = total_m2_rec_cyt;
    RM2_free_Cyt = total_act_m2_rec_cyt;
    
    double total_gi_nr_cyt = 90e6;
    double total_gi_cyt; //total concentration of gi - uM
    double total_avail_gi_cyt;
    total_gi_cyt = 10;
    total_avail_gi_cyt = total_gi_cyt;
    
    
    Gi_free_Cyt = total_avail_gi_cyt - (y[30] + y[32]);
    //R_M2
    
    RGi_Cyt = (RM2_free_Cyt * Gi_free_Cyt )/(Gi_free_Cyt+K_C2);
    
    L_ach_RM2_Cyt = (L_ach * (RM2_free_Cyt-RGi_Cyt))/(L_ach+K_L2) ;
    
    LRGi_Cyt = (L_ach_RM2_Cyt*(Gi_free_Cyt-RGi_Cyt))/((Gi_free_Cyt-RGi_Cyt)+(K_C2*K_H2/K_L2)) + ((L_ach * RGi_Cyt)/(L_ach + K_H2));
    
    RM2_tot_Cyt = RM2_free_Cyt + L_ach_RM2_Cyt + LRGi_Cyt + RGi_Cyt;
    
    //Gi
    ydot[30] = (LRGi_Cyt * kact2_cyt2 + RGi_Cyt * kact1_cyt2 - y[30] * khydr_cyt2);
    
    ydot[31] = (LRGi_Cyt * kact2_cyt2 + RGi_Cyt *kact1_cyt2 - y[32]* y[31] * kreas_cyt2);
    
    ydot[32] = (y[30] * khydr_cyt2 - y[32] * y[31] * kreas_cyt2);
    
    if (y[30]<0){
        y[30] = 0;
    }
    
    if (y[31]<0){
        y[31] = 0;
    }
    
    if (y[32]<0){
        y[32] = 0;
    }
    
    
    
    
    //cAMP in Cav
    const double ATP = 5E03; // microM
    const double K_mATP = 315 ; //microM
    const double AF56 = 500;// mg purified protein/mg membrane protein
    const double MW_AC56 = 130; //KDa Molecular weight
    
    double k_AC56, dcAMP_AC56dt;
    double Vmax_cav;
    double basalAC56 = 0;
    
    
    double dcAMP_PDE2_Cavdt, dcAMP_PDE3_Cavdt, dcAMP_PDE4_Cavdt;
    
    
    double PDE2_Cav = 4.5;//uM
    double PDE3_Cav = 5.6;//uM
    double PDE4_Cav = 2.0;//uM
    
    double k_PDE2_cav = 20;// uM
    double K_mPDE2_cav = 50; //uM
    double k_PDE3_cav = 1.25; //uM
    double K_mPDE3_cav = 0.08; //uM
    double k_PDE4_cav = 2.5; //uM
    double K_mPDE4_cav = 2.2; // uM
    
    double dFlux2_cav;
    double dFlux3_cav;
    
    double ac5_6_nr_cav = 3.6e5;
    double AC56_Cav = (ac5_6_nr_cav/(6.23*cav_volume))*1e-11;
    
    
    Vmax_cav = ( (  ( 3.8234 * pow( y[0],0.9787 ) )
                  / ( 0.1986 + pow( y[0],0.9787 ) )
                  + 0.7 )
                * ( 1 + ( (  ( -1.0061 * pow( y[3],0.8356 ) )
                           / ( 0.1918 + pow( y[3],0.8356 ) ) ) / 1.4432 ) ) * 1e-3) ;
    Vmax_cav = ((Vmax_cav + basalAC56) * AF56 * AC56_Cav * MW_AC56)/60;
    
    dcAMP_AC56dt = ((Vmax_cav*ATP)/(K_mATP+ATP));
    
    //Phosphodiesterase
    dcAMP_PDE2_Cavdt = (((k_PDE2_cav * PDE2_Cav) * y[33])/(K_mPDE2_cav + y[33]));
    
    dcAMP_PDE3_Cavdt = (((k_PDE3_cav * PDE3_Cav) * y[33])/(K_mPDE3_cav + y[33]));
    
    dcAMP_PDE4_Cavdt = (((k_PDE4_cav * PDE4_Cav) * y[33])/(K_mPDE4_cav + y[33]));
    
    dFlux2_cav = (1e6*flux2*(y[33]-y[35])/cav_volume);
    dFlux3_cav = (1e6*flux3*(y[33]-y[34])/cav_volume);
    
    
    ydot[33] = (dcAMP_AC56dt - (dcAMP_PDE2_Cavdt + dcAMP_PDE3_Cavdt + dcAMP_PDE4_Cavdt) - dFlux2_cav - dFlux3_cav);
    
    
    // Calculate cAMP concentration in Extra-caveolar
    
    
    const double AF47 = 130;// mg purified protein/mg membrane protein
    const double MW_AC47 = 130; //KDa Molecular weight
    double k_AC47_Ecav;
    double Vmax_ecav;
    double basalAC47 = 0;
    
    //double cAMP_Ecav = 0;
    
    
    double dcAMP_AC47_Ecavdt, dcAMP_PDE2_Ecavdt, dcAMP_PDE4_Ecavdt;
    double k_PDE2_ecav = 20;// uM
    double k_PDE4_ecav = 2.5; //uM
    double K_mPDE2_ecav = 50; //uM
    double K_mPDE4_eacv = 2.2; // uM
    double PDE2_Ecav = 0.002;//uM
    double PDE4_Ecav = 0.01;//uM
    
    double dFlux1_ecav;
    double dFlux3_ecav;
    
    double ac4_7_nr_ecav = 0.475e5;
    double AC47_Ecav = (ac4_7_nr_ecav/(6.23*ecav_volume))*1e-11;
    
    
    Vmax_ecav = (0.063+(2.01*pow(y[6]*1e3,1.0043))/(31.544+pow(y[6]*1e3,1.0043)))*(1+(((49.1*pow(y[10]*1e3,0.8921))/(25.44+pow(y[10]*1e3,0.8921)))/3.01))* 1e-3 ;
    Vmax_ecav = ((Vmax_ecav + basalAC47) * AF47 * AC47_Ecav *MW_AC47)/ 60;
    
    dcAMP_AC47_Ecavdt = ( (Vmax_ecav * ATP) /(K_mATP + ATP));
    
    //Phosphodiesterase
    dcAMP_PDE2_Ecavdt = (((k_PDE2_ecav * PDE2_Ecav) * y[34])/(K_mPDE2_ecav + y[34]));
    dcAMP_PDE4_Ecavdt = (((k_PDE4_ecav * PDE4_Ecav) * y[34])/(K_mPDE4_eacv + y[34]));
    
    dFlux1_ecav = (1e6*flux1*(y[34]-y[35])/ecav_volume);
    dFlux3_ecav = (1e6*flux3*(y[34]-y[33])/ecav_volume);
    
    ydot[34] = (dcAMP_AC47_Ecavdt - (dcAMP_PDE2_Ecavdt + dcAMP_PDE4_Ecavdt) - dFlux1_ecav - dFlux3_ecav);
    
    
    // cAMP concentration Cytoplasmic
    
    
    const double MW = 130; //KDa Molecular weight
    
    double Vmax_AC56;
    double Vmax_AC47;
    
    
    double dcAMP_AC56_Cytdt;
    double dcAMP_AC47_Cytdt;
    double dcAMP_PDE2_Cytdt, dcAMP_PDE3_Cytdt, dcAMP_PDE4_Cytdt;
    
    double k_PDE1 = 2;
    double k_PDE2 = 20;// /s
    double k_PDE3 = 1.25; // /s
    double k_PDE4 = 2.5; // /s
    
    double K_mPDE1 = 0.6;
    double K_mPDE2 = 50; //uM
    double K_mPDE3 = 0.08; //uM
    double K_mPDE4 = 2.2; // uM
    
    double PDE1_Cyt = 0.035;
    double PDE2_Cyt = 0.085;//uM
    double PDE3_Cyt = 0.113;//uM
    double PDE4_Cyt = 0.027;//uM
    
    double dFlux;
    
    double ac5_6_nr = 1.13e5;
    double ac4_7_nr = 8e3;
    double AC56_Cyt, AC47_Cyt;
    
    
    AC56_Cyt = (ac5_6_nr/(6.23*myop_membr_volume))*1e-11;
    AC47_Cyt = (ac4_7_nr/(6.23*myop_membr_volume))*1e-11;
    
    Vmax_AC56 = ( ( 3.8234 * pow( y[12],0.9787))/(0.1986 + pow( y[12],0.9787)) + 0.7 )
    * (1+(((-1.0061 * pow( y[30],0.8356 )) / (0.1918 + pow(y[30],0.8356)))/1.4432))*1e-3;
    
    Vmax_AC56 = ( ( Vmax_AC56 + basalAC56) * AF56 * AC56_Cyt * MW)/60;
    
    dcAMP_AC56_Cytdt = (Vmax_AC56 *ATP)/( K_mATP + ATP);
    
    Vmax_AC47 = (0.063+(2.01 * pow( y[12]*1e3,1.0043)) / (31.544 + pow( y[12]*1e3,1.0043)))
    * (1+ (((49.1 * pow(y[31]*1e3,0.8921)) / (25.44 + pow(y[31]*1e3,0.8921)))
           / 3.01))* 1e-3;
    
    Vmax_AC47 = ( ( Vmax_AC47 + basalAC47 ) * AF47 * AC47_Cyt *MW)/ 60;
    
    dcAMP_AC47_Cytdt = (Vmax_AC47 * ATP )/( K_mATP + ATP);
    
    //	//Phosphodiesterase
    dcAMP_PDE2_Cytdt = (((k_PDE2 * PDE2_Cyt) * y[35])/(K_mPDE2 + y[35]));
    dcAMP_PDE3_Cytdt = (((k_PDE3 * PDE3_Cyt) * y[35])/(K_mPDE3 + y[35]));
    dcAMP_PDE4_Cytdt = (((k_PDE4 * PDE4_Cyt) * y[35])/(K_mPDE4 + y[35]));
    
    dFlux = (((1e6*flux1*(y[34] - y[35]))+
              (1e6*flux2*(y[33] - y[35])))/myop_volume);
    
    
    ydot[35] = ( (dcAMP_AC56_Cytdt + dcAMP_AC47_Cytdt) - ( dcAMP_PDE2_Cytdt + dcAMP_PDE3_Cytdt + dcAMP_PDE4_Cytdt) + dFlux);
    
    //// PKA module
    double PKI = p[34]*p[38]/(p[38]+y[16]+y[17]);
    double A2RC_I = (y[16]/p[37])*y[16]*(1+PKI/p[38]);
    double A2R_I = y[16]*(1+PKI/p[38]);
    double A2RC_II = (y[17]/p[37])*y[17]*(1+PKI/p[38]);
    double A2R_II = y[17]*(1+PKI/p[38]);
    double ARC_I = (p[35]/y[15])*A2RC_I;
    double ARC_II = (p[35]/y[15])*A2RC_II;
    ydot[15] = y[33]-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y[15];
    double PKAtemp = p[35]*p[36]/p[37]+p[35]*y[15]/p[37]+y[15]*y[15]/p[37];
    ydot[16] = 2*p[32]*y[15]*y[15]-y[16]*(1+PKI/p[38])*(PKAtemp*y[16]+y[15]*y[15]);
    ydot[17] = 2*p[33]*y[15]*y[15]-y[17]*(1+PKI/p[38])*(PKAtemp*y[17]+y[15]*y[15]);
    double y15, y16, y17, y15sq;
    double ip37 = 1./p[37];
    double ip38 = 1./p[38];
    
    // Solving DAEs for y[15], y[16], y[17].
    iter = 0;
    double mk;
    while ( iter < 1000 && ( fabs( ydot[15] ) > myeps || fabs( ydot[16] ) > myeps || fabs( ydot[17] ) > myeps ) )
    {
        if ( L_iso == 0) {
            mk = 0.006;
        } else {
            mk = 0.35; // Increase convergent speed when close to fixed point.
        }
        
        PKI = p[34]*p[38]/(p[38]+y[16]+y[17]);
        A2RC_I = (y[16]/p[37])*y[16]*(1+PKI/p[38]);
        A2R_I = y[16]*(1+PKI/p[38]);
        A2RC_II = (y[17]/p[37])*y[17]*(1+PKI/p[38]);
        A2R_II = y[17]*(1+PKI/p[38]);
        ARC_I = (p[35]/y[15])*A2RC_I;
        ARC_II = (p[35]/y[15])*A2RC_II;
        ydot[15] = y[33]-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y[15];
        while ( fabs( mk * ydot[15] / y[15] ) > 0.1 ) {
            mk = 0.1 * mk;
        }
        y[15] += mk * ydot[15];
        y15sq = y[15] * y[15];
        PKAtemp = ( p[35] * ( p[36] + y[15] ) + y15sq )* ip37;
        
        y16 = y[16];
        y[16] = 2*p[32]*y15sq / ( (1+PKI*ip38)*(PKAtemp*y[16]+y15sq) );
        ydot[16] = y[16] - y16;
        y17 = y[17];
        y[17] = 2*p[33]*y15sq / ( (1+PKI*ip38)*(PKAtemp*y[17]+y15sq) );
        ydot[17] = y[17] - y17;
        iter++;
    }
    ydot[15] = 0;
    ydot[16] = 0;
    ydot[17] = 0;
    // End solving DAEs for y[9], y[10], y[11].
    
    // end PKA module
    
    //// PLB module
    double PLB = p[40]-y[18];
    double PLB_PHOSPH = p[43]*y[16]*PLB/(p[44]+PLB);
    double PLB_DEPHOSPH = p[45]*y[21]*y[18]/(p[46]+y[18]);
    ydot[18] = PLB_PHOSPH-PLB_DEPHOSPH;
    
    double Inhib1 = p[42]-y[19];
    double Inhib1p_PP1 = y[20]*y[21]/p[51];
    double Inhib1_PHOSPH = p[47]*y[16]*Inhib1/(p[48]+Inhib1);
    double Inhib1_DEPHOSPH = p[49]*y[19]/(p[50]+y[19]);
    ydot[19] = Inhib1_PHOSPH-Inhib1_DEPHOSPH;
    ydot[20] = y[19]-Inhib1p_PP1-y[20];
    ydot[21] = p[41]-Inhib1p_PP1-y[21];
    
    // Solving DAEs for y[20], y[21].
    double ip51 = 1. / p[51];
    double y20, y21;
    
    y20 = y[20];
    y21 = y[21];
    
    iter = 0;
    while ( iter < 1000 && ( fabs( ydot[20] ) > myeps || fabs( ydot[21] ) > myeps ) )
    {
        y20 = y[19] / ( 1 + y[21] * ip51 );
        y21 = p[41] / ( 1 + y20 * ip51 );
        
        y[20] = y[19] / ( 1 + y21 * ip51 );
        y[21] = p[41] / ( 1 + y[20] * ip51 );
        
        ydot[20] = y[20] - y20;
        ydot[21] = y[21] - y21;
        
        iter++;
    }
    ydot[20] = 0;
    ydot[21] = 0;
    // End solving DAEs for y[20], y[21].
    
    
    // end PLB module
    
    //// LCC module
    double PKAClcc = (p[53]/p[33])*y[17];
    double LCCa = p[52]-y[22];
    double LCCa_PHOSPH = p[39]*p[56]*PKAClcc*LCCa/(p[57] + p[39]*LCCa);
    double LCCa_DEPHOSPH = p[39]*p[60]*p[55]*y[22]/(p[61]+p[39]*y[22]);
    ydot[22] = LCCa_PHOSPH - LCCa_DEPHOSPH;
    
    double LCCb = p[52]-y[23];
    double LCCb_PHOSPH = p[39]*p[56]*PKAClcc*LCCb/(p[57]+p[39]*LCCb);
    double LCCb_DEPHOSPH = p[39]*p[58]*p[54]*y[23]/(p[59]+p[39]*y[23]);
    ydot[23] = LCCb_PHOSPH-LCCb_DEPHOSPH;
    // end LCC module
    
    //// RyR module
    double PKACryr = (p[63]/p[33])*y[17];
    double RyR = p[62]-y[24];
    double RyRPHOSPH = p[39]*p[66]*PKACryr*RyR/(p[67]+p[39]*RyR);
    double RyRDEPHOSPH1 = p[39]*p[68]*p[64]*y[24]/(p[69]+p[39]*y[24]);
    double RyRDEPHOSPH2A = p[39]*p[70]*p[65]*y[24]/(p[71]+p[39]*y[24]);
    ydot[24] = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A;
    // end RyR module
    
    //// TnI module
    double TnI = p[72]-y[25];
    double TnIPHOSPH = p[74]*y[16]*TnI/(p[75]+TnI);
    double TnIDEPHOSPH = p[76]*p[73]*y[25]/(p[77]+y[25]);
    ydot[25] = TnIPHOSPH-TnIDEPHOSPH;
    // end TnI module
    
    //// Iks module
    double IksYot = y[26]*y[27]/p[80];           // [uM]
    ydot[26] = p[78] - IksYot - y[26];    // [uM]
    ydot[27] = p[79] - IksYot - y[27];    // [uM]
    
    // Solving DAEs for y[26], y[27].
    double ip80 = 1. / p[80];
    double y26, y27;
    
    y26 = y[26];
    y27 = y[27];
    
    iter = 0;
    while ( iter < 1000 && ( fabs( ydot[26] ) > myeps || fabs( ydot[27] ) > myeps ) )
    {
        y26 = p[78] / ( 1 + y[27] * ip80 );
        y27 = p[79] / ( 1 + y26 * ip80 );
        
        y[26] = p[78] / ( 1 + y27 * ip80 );
        y[27] = p[79] / ( 1 + y[26] * ip80 );
        
        ydot[26] = y[26] - y26;
        ydot[27] = y[27] - y27;
        
        iter++;
    }
    ydot[26] = 0;
    ydot[27] = 0;
    // End solving DAEs for y[26], y[27].
    
    double PKACiks = (IksYot/p[78])*(p[81]/p[33])*y[17];
    double PP1iks = (IksYot/p[78])*p[82];
    double Iks = p[78]-y[28];
    double IKS_PHOSPH = p[39]*p[83]*PKACiks*Iks/(p[84]+p[39]*Iks);
    double IKS_DEPHOSPH = p[39]*p[85]*PP1iks*y[28]/(p[86]+p[39]*y[28]);
    ydot[28] = IKS_PHOSPH-IKS_DEPHOSPH;
    
    //// CFTR module (included 04/30/10)
    double CFTRn = p[87] - y[29];  // Non-phos = tot - phos
    double PKAC_CFTR = (p[88]/p[33])*y[17];    // (PKACFTRtot/PKAIItot)*PKAIIact
    double CFTRphos = p[39]*CFTRn*PKAC_CFTR*p[90]/(p[91]+p[39]*CFTRn);
    double CFTRdephos = p[89]*p[92]*p[39]*y[29]/(p[93] + p[39]*y[29]);
    ydot[29] = CFTRphos - CFTRdephos;
    
    // PLM module (included 09/16/12)
    double PLM = p[94]-y[36];
    double PLM_PHOSPH = p[95]*y[16]*PLM/(p[96]+PLM);
    double PLM_DEPHOSPH = p[97]*y[21]*y[36]/(p[98]+y[36]);
    ydot[36] = PLM_PHOSPH-PLM_DEPHOSPH;
    
    //// Gather odes
    // Need to convert all ydot terms that are ODEs (not DAEs) to miliseconds
    // odes = [4,5,6,7,8,9,13,14,15,19,20,23,24,25,26,29,30];
    // ydot(odes) = ydot(odes).*1e-3;
    for ( int i = 0; i < 37; i++ ) {
        ydot[i] = ydot[i] * 0.001;
    }
}
#endif
