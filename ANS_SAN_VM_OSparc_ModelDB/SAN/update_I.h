//function dF = update_I(t,F,con)
//    % Joachim Behar, 04-05-2017
//    % Bioenergetics and Bioelectric Systems Lab
//    % Technion, Haifa, Israel
//    %
//    %
//    %    This program is free software; you can redistribute it and/or modify
//    %    it under the terms of the GNU General Public License as published by
//    %    the Free Software Foundation; either version 2 of the License, or
//    %    (at your option) any later version.
//    %
//    %    This program is distributed in the hope that it will be useful,
//    %    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    %    GNU General Public License for more details.
//    %
//    %    You should have received a copy of the GNU General Public License
//    %    along with this program; if not, write to the Free Software
//    %    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
struct cell_c {
    double Cai, Ca_sub, Ca_JSR, Ca_NSR, Cam;
};

struct cell_ac {
    double cAMP, PLB_p, PKA, ATP, k_oCa;
};
#include <iostream>
using namespace std;

#include "update_PKA.h"

#include "update_I_CaL.h"
#include "update_I_CaT.h"
#include "update_I_Kr.h"
#include "update_I_Ks.h"
#include "update_I_4AP.h"
#include "update_I_f.h"
#include "update_I_st.h"
#include "update_I_bNa.h"
#include "update_I_NaK.h"
#include "update_I_bCa.h"
#include "update_I_NCX.h"
#include "update_I_KACh.h"


//void update_I( double * dF, double t, double * F, cell_con * con );
//
//void update_I( double * dF, double t, double * F, cell_con * con ){
void update_I( double t, double * F, cell_con * con, double * dF);

void update_I( double t, double * F, cell_con * con, double * dF){
    
    
    //
    //    if rem(t,1000)==0 && t>1000
    //        disp('-->');
    //    end
    //
    //    % == adding a timer
    //    global elapsedtime
    //
    //    if isempty(elapsedtime)
    //        elapsedtime = tic;
    //    end
    //
    //%     if toc(elapsedtime) > 5000
    //%         error('Stopped. Taking too long.');
    //%     end
    //
    //    %disp(num2str(t));
    //
    //    % = concentration variables
    //    c.Cai = F(1); c.Ca_sub = F(2); c.Ca_JSR = F(3); c.Ca_NSR = F(4); c.Cam = F(5); % concentrations variables
    cell_c c;
    c.Cai = F[0]; c.Ca_sub = F[1]; c.Ca_JSR = F[2]; c.Ca_NSR = F[3]; c.Cam = F[4]; //% concentrations variables
    //    f_TMC = F(6); f_TMM = F(7); f_CMi = F(8); f_CMs = F(9);
    //    f_CQ = F(10); R = F(11); OO = F(12); S = F(13); RI = F(14); % buffers variables
    //    V_m = F(15); % membrane potential
    //    d_L = F(16); f_L = F(17); f_Ca = F(18); p_aF = F(19);
    //    p_aS = F(20); p_i = F(21); n = F(22); y = F(23); d_T = F(24);
    //    f_T = F(25); q = F(26); r = F(27); q_a = F(28); q_i = F(29); % gating variables
    //    ac.cAMP = F(30); ac.PLB_p = F(31); % PKA signaling variables
    //    A = F(32); TT = F(33); U = F(34); SL = F(35); % force variables - % FIXME SL -> it is a constant in the code
    //    w = F(36);
    double f_TMC = F[5]; double f_TMM = F[6]; double f_CMi = F[7]; double f_CMs = F[8];
    double f_CQ = F[9]; double R = F[10]; double OO = F[11]; double S = F[12]; double RI = F[13]; // % buffers variables
    double V_m = F[14]; // % membrane potential
    double d_L = F[15]; double f_L = F[16]; double f_Ca = F[17]; double p_aF = F[18];
    double p_aS = F[19]; double p_i = F[20]; double n = F[21]; double y = F[22]; double d_T = F[23];
    double f_T = F[24]; double q = F[25]; double r = F[26]; double q_a = F[27]; double q_i = F[28]; // % gating variables
    cell_ac ac;
    ac.cAMP = F[29]; ac.PLB_p = F[30]; // % PKA signaling variables
    double A = F[31]; double TT = F[32]; double U = F[33]; double SL = F[34]; // % force variables - % FIXME SL -> it is a constant in the code
    double w = F[35];
    //
    //    if con->ALL_VAR; f_cAB = F(43); end
    double f_cAB;
    if( con->ALL_VAR ) {
        f_cAB = F[42];
        // end
    }
    
    double df_cAB;
    //
    //%% ==
    //%     if t>900*1000 & t<905*1000
    //%         cAB_Conc = con->cAB_Conc;
    //%         df_cAB = con->cAB_on_rate*(ac.cAMP/600)*(1-f_cAB)-con->cAB_off_rate*f_cAB; % cAMP Buffer
    //%     else
    //%         cAB_Conc = con->cAB_Conc;
    //%         df_cAB = con->cAB_on_rate*(ac.cAMP/600)*(1-f_cAB); % cAMP Buffer
    //%     end
    //%% ==
    //
    //    % = AC-cAMP-PKA cycling
    //    ac.PKA = update_PKA(ac.cAMP);
    ac.PKA = update_PKA( ac.cAMP );
    //    ac.ATP = con->ATP_max*(con->kATP*(ac.cAMP*100/con->cAMPb)^con->n_ATP/...
    //    (con->k_ATP05+(ac.cAMP*100/con->cAMPb)^con->n_ATP)-con->K_ATPmin)/100; % ATP outputed in mM
    ac.ATP = con->ATP_max * ( con->kATP * pow( ( ac.cAMP * 100. / con->cAMPb ) , con->n_ATP ) /
                             ( con->k_ATP05 + pow( ( ac.cAMP * 100. / con->cAMPb ) , con->n_ATP ) ) - con->K_ATPmin ) / 100. ; // % ATP outputed in mM
    //    if con->ISO==0; k_iso = 0; else k_iso = 0.0070+0.1181.*con->ISO.^0.8664./(48.1212.^0.8664+con->ISO.^0.8664); end;
    double k_iso;
    if( con->ISO==0 ) {
        k_iso = 0;
    } else {
        k_iso = 0.0070 + 0.1181 * pow( con->ISO, 0.8664 ) / ( pow( 48.1212 , 0.8664 ) + pow( con->ISO , 0.8664 ) );
        // end;
    }
    //    k_CCh = 0.0146.*con->CCh.^1.4402./(51.7331.^1.4402+con->CCh.^1.4402);
    double k_CCh = 0.0146 * pow( con->CCh, 1.4402 ) / ( pow( 51.7331, 1.4402 ) + pow( con->CCh, 1.4402 ) );
    //    k_1 = con->K_ACI+con->K_AC/(1+exp((con->K_Ca-con->k_bCM*f_CMi/(con->k_fCM*(1-f_CMi)))/con->K_AC_Ca)); % transformation of ATP into cAMP
    double k_1 = con->K_ACI + con->K_AC / ( 1. + exp( ( con->K_Ca - con->k_bCM * f_CMi / ( con->k_fCM * ( 1. - f_CMi ) ) )
                                                     / con->K_AC_Ca ) ); // % transformation of ATP into cAMP
    //    k_2 = 1.1*237.9851.*ac.cAMP.^5.1010./(20.1077.^6.1010+ac.cAMP.^6.1010);
    double k_2 = 1.1 * 237.9851 * pow( ac.cAMP , 5.1010 ) / ( pow( 20.1077, 6.1010 ) + pow( ac.cAMP, 6.1010 ) );
    //    k_3 = con->k_PKA*(ac.cAMP^(con->n_PKA-1))/(con->k_PKA_cAMP^(con->n_PKA)+ac.cAMP^con->n_PKA); % cAMP transformation into PKA
    double k_3 = con->k_PKA * ( pow( ac.cAMP , (con->n_PKA-1) ) )
    / ( pow( con->k_PKA_cAMP , (con->n_PKA) ) + pow( ac.cAMP , con->n_PKA ) ); // % cAMP transformation into PKA
    //    k_4 = (con->k_PLBp*ac.PKA^con->n_PLB)/(con->k_PKA_PLB^con->n_PLB+ac.PKA^(con->n_PLB)); % PLB phosphorylation
    double k_4 = ( ( con->k_PLBp * pow( ac.PKA , con->n_PLB ) )
                  / ( pow( con->k_PKA_PLB , con->n_PLB ) + pow( ac.PKA , (con->n_PLB) ) ) ); // % PLB phosphorylation
    //    k_5 = con->k_PP1*con->PP1*ac.PLB_p/(con->k_pp1_PLB+ac.PLB_p);
    double k_5 = con->k_PP1 * con->PP1 * ac.PLB_p / ( con->k_pp1_PLB + ac.PLB_p );
    //
    //    dcAMP = k_iso*(ac.ATP*0.6*1e3)+k_1*(ac.ATP*0.6*1e3)-k_2*ac.cAMP-k_3*ac.cAMP-k_CCh*(ac.ATP*0.6*1e3);
    double dcAMP = ( k_iso * ( ac.ATP * 0.6 * 1e3 ) + k_1 * ( ac.ATP * 0.6 * 1e3 )
                    - k_2 * ac.cAMP - k_3 * ac.cAMP - k_CCh * ( ac.ATP * 0.6 * 1e3 ) );
    //    %dcAMP = dcAMP/(60e3)-(cAB_Conc*0.6*1e3)*df_cAB; % Transfer 1/min to 1/msec + cAMP caging
    //    dcAMP = dcAMP/(60e3);
    dcAMP = dcAMP / ( 60e3 );
    //    dPLB_p = k_4-k_5;
    double dPLB_p = k_4 - k_5;
    //    dPLB_p = dPLB_p/(60e3); % Transfer 1/min to 1/msec
    dPLB_p = dPLB_p / ( 60e3 ); // % Transfer 1/min to 1/msec
    //
    //    % = compute currents
    //    [I_CaL,a_l,b_l,a_fl,b_fl,a_fCa,b_fCa,ac] = update_I_CaL(con,V_m,c,ac,d_L,f_L,f_Ca); % L-type Ca2+ current (I_CaL)
    //    [I_CaT,a_dT,b_dT,a_fT,b_fT] = update_I_CaT(con,V_m,d_T,f_T); % T-type Ca2+ current (ICaT)
    //    [I_Kr,a_paF,b_paF,a_paS,b_paS,a_pi,b_pi] = update_I_Kr(con,V_m,p_aF,p_aS,p_i); % rapidly activating delayed rectifier K+ current (I_Kr)
    //    [I_Ks,a_n,b_n] = update_I_Ks(con,V_m,ac,n); % Slowly activating delayed rectifier K+ current
    //    [I_to,I_sus,a_q,b_q,a_r,b_r] = update_I_4AP(con,V_m,q,r); % 4-aminopyridine-sensitive currents (I_4AP = I_to + I_sus)
    //    [I_f,~,~,a_y,b_y]= update_I_f(con,V_m,ac,y); % Hyperpolarization-activated, funny current (I_f)
    //    [I_st,a_qa,b_qa,a_qi,b_qi] = update_I_st(con,V_m,q_a,q_i); % Sustained inward current (I_st)
    //    I_bNa = update_I_bNa(con,V_m); % Na + -dependent background current (I_bNa)
    //    I_NaK = update_I_NaK(con,V_m,ac); % Na + -K + pump current (I_NaK)
    //    I_bCa = update_I_bCa(con,V_m); % Ca 2+ - background current (I_bCa)
    //    [I_NCX] = update_I_NCX(con,V_m,c); % Na + -Ca 2+ exchanger current (I_NCX)
    //    [I_KACh,a_w,b_w] = update_I_KACh(con,V_m,w); % Acetylcholine-activated K+ current
    
    double I_CaL , I_CaT , I_bCa , I_f , I_st , I_Kr;
    double I_Ks , I_NaK , I_NCX , I_bNa , I_to , I_sus , I_KACh;
    double a_l, b_l, a_fl, b_fl, a_fCa, b_fCa;
    double a_dT, b_dT, a_fT, b_fT ;
    double a_paF, b_paF, a_paS, b_paS, a_pi, b_pi;
    double a_n, b_n;
    double a_q, b_q, a_r, b_r;
    double a_y, b_y;
    double a_qa, b_qa, a_qi, b_qi;
    double a_w, b_w;
    
    
    update_I_CaL( con, V_m, &c, &ac, d_L, f_L, f_Ca, &I_CaL, &a_l, &b_l, &a_fl, &b_fl, &a_fCa, &b_fCa ); //% L-type Ca2+ current (I_CaL)
    update_I_CaT( con, V_m, d_T, f_T, &I_CaT, &a_dT, &b_dT, &a_fT, &b_fT ); // % T-type Ca2+ current (ICaT)
    update_I_Kr( con, V_m, p_aF, p_aS, p_i, &I_Kr, &a_paF, &b_paF, &a_paS, &b_paS, &a_pi, &b_pi ); // % rapidly activating delayed rectifier K+ current (I_Kr)
    update_I_Ks( con, V_m, &ac, n, &I_Ks, &a_n, &b_n ); // % Slowly activating delayed rectifier K+ current
    
    update_I_4AP( con, V_m, q, r, &I_to, &I_sus, &a_q, &b_q, &a_r, &b_r ); // % 4-aminopyridine-sensitive currents (I_4AP = I_to + I_sus)
    update_I_f( con, V_m, &ac, y, &I_f, &a_y, &b_y); // % Hyperpolarization-activated, funny current (I_f)
    update_I_st( con, V_m, q_a, q_i, &I_st, &a_qa, &b_qa, &a_qi, &b_qi ); // % Sustained inward current (I_st)
    update_I_bNa( con, V_m, &I_bNa ); // % Na + -dependent background current (I_bNa)
    update_I_NaK( con, V_m, &ac, &I_NaK ); // % Na + -K + pump current (I_NaK)
    update_I_bCa( con, V_m, &I_bCa ); // % Ca 2+ - background current (I_bCa)
    update_I_NCX( con, V_m, &c, &I_NCX ); // % Na + -Ca 2+ exchanger current (I_NCX)
    update_I_KACh( con, V_m, w, &I_KACh, &a_w, &b_w ); // % Acetylcholine-activated K+ current
    
    
    //
    //    % = RyR phosphorylation
    //    k_CaSR = con->MaxSR-(con->MaxSR-con->MinSR)/(1+(con->EC_50SR/c.Ca_JSR)^con->HSR);
    double k_CaSR = con->MaxSR - ( con->MaxSR - con->MinSR ) / ( 1. + pow( ( con->EC_50SR / c.Ca_JSR ) , con->HSR ) );
    //    ac.k_oCa = con->k_oCa_max*(con->RyR_min-con->RyR_max*ac.PKA^con->n_RyR/...
    ac.k_oCa = con->k_oCa_max * ( con->RyR_min - con->RyR_max * pow( ac.PKA , con->n_RyR ) /
                                 ( pow( con->k_05Ry , con->n_RyR ) + pow( ac.PKA , con->n_RyR ) ) + 1. );
    //    k_oSRCa = ac.k_oCa/k_CaSR; % `k_oCa' is modulated by PKA activity
    double k_oSRCa = ac.k_oCa / k_CaSR; // % `k_oCa' is modulated by PKA activity
    //    k_iSRCa = con->k_iCa*k_CaSR;
    double k_iSRCa = con->k_iCa * k_CaSR;
    //    dR = (con->k_im*RI-k_iSRCa*c.Ca_sub*R)-(k_oSRCa*c.Ca_sub^2*R-con->k_om*OO);
    double dR = ( con->k_im * RI - k_iSRCa * c.Ca_sub * R ) - ( k_oSRCa * c.Ca_sub * c.Ca_sub * R - con->k_om * OO );
    //    dOO = (k_oSRCa*c.Ca_sub^2*R-con->k_om*OO)-(k_iSRCa*c.Ca_sub*OO-con->k_im*S);
    double dOO = ( k_oSRCa * c.Ca_sub * c.Ca_sub * R - con->k_om * OO ) - ( k_iSRCa * c.Ca_sub * OO - con->k_im * S );
    //    dS = (k_iSRCa*c.Ca_sub*OO-con->k_im*S)-(con->k_om*S-k_oSRCa*c.Ca_sub^2*RI);
    double dS = ( k_iSRCa * c.Ca_sub * OO - con->k_im * S ) - ( con->k_om * S - k_oSRCa * c.Ca_sub * c.Ca_sub * RI );
    //    dRI = (con->k_om*S-k_oSRCa*c.Ca_sub^2*RI)-(con->k_im*RI-k_iSRCa*c.Ca_sub*R);
    double dRI = ( con->k_om * S - k_oSRCa * c.Ca_sub * c.Ca_sub * RI ) - ( con->k_im * RI - k_iSRCa * c.Ca_sub * R );
    //
    //    % = Force
    double Ve = 0;
    double dSL = -Ve;
    //    NXB = (SL-con->SL_lo)/2*1e3*(TT+U)*con->Nc;
    double NXB = ( SL - con->SL_lo ) / 2 * 1e3 * ( TT + U ) * con->Nc;
    //    K_Ca = con->FK0+con->FK1*(NXB^con->FN)/(con->FK05^con->FN+NXB^con->FN);
    double K_Ca = con->FK0 + con->FK1 * pow( NXB , con->FN ) / ( pow( con->FK05 , con->FN ) + pow( NXB , con->FN ) );
    double k_l = con->Fkl / K_Ca;
    //    dA = con->Fkl*c.Cai*(1-A-TT-U)-(con->Ff+k_l)*A+(con->Fg0+con->Fg1*Ve)*TT;
    double dA = con->Fkl * c.Cai * ( 1. - A - TT - U ) - ( con->Ff + k_l ) * A + ( con->Fg0 + con->Fg1 * Ve ) * TT;
    //    dTT = con->Ff*A-(con->Fg0+con->Fg1*Ve+k_l)*TT+con->Fkl*c.Cai*U;
    double dTT = con->Ff * A - ( con->Fg0 + con->Fg1 * Ve + k_l ) * TT + con->Fkl * c.Cai * U;
    //    dU = k_l*TT-(con->Fg0+con->Fg1*Ve+con->Fkl*c.Cai)*U;
    double dU = k_l * TT - ( con->Fg0 + con->Fg1 * Ve + con->Fkl * c.Cai ) * U;
    //
    
    
    con->I_CaL = I_CaL;
    con->I_NCX = I_NCX;
    con->I_KACh = I_KACh;
    con->I_f = I_f;
    
    
    //    % = overall current
    double I = I_CaL + I_CaT + I_bCa + I_f + I_st + I_Kr
    + I_Ks + I_NaK + I_NCX + I_bNa + I_to + I_sus + I_KACh;
    double dV = -( 1. / con->C ) * I;
    //
    //    % = intracellular Ca2+ flux
    //    modfac = 1.6980.*ac.PLB_p.^13.5842./(0.2240.^13.5842+ac.PLB_p.^13.5842); % for ACh
    double modfac = 1.6980 * pow( ac.PLB_p , 13.5842 ) / ( pow( 0.2240 , 13.5842 ) + pow( ac.PLB_p , 13.5842 ) ); // % for ACh
    //    if ac.PLB_p>0.23
    //        % for ISO
    //        modfac = 3.3931.*ac.PLB_p.^4.0695./(0.2805.^4.0695+ac.PLB_p.^4.0695)-.0952;
    //    end
    if( ac.PLB_p > 0.23 ){
        //        % for ISO
        modfac = 3.3931 * pow( ac.PLB_p , 4.0695 ) / ( pow( 0.2805 , 4.0695 ) + pow( ac.PLB_p , 4.0695 ) ) - .0952;
        //    end
    }
    //    %modfac = 1;
    //
    //    j_SRCarel = con->k_s*OO*(c.Ca_JSR-c.Ca_sub); % Ca2+ fluxes in the SR
    double  j_SRCarel = con->k_s * OO * ( c.Ca_JSR - c.Ca_sub ); // % Ca2+ fluxes in the SR
    //    j_Ca_dif = (c.Ca_sub-c.Cai)/con->tho_difCa;
    double j_Ca_dif = ( c.Ca_sub - c.Cai ) / con->tho_difCa;
    //    j_up = 0.9*con->P_up*modfac/(1+con->K_up/c.Cai);
    double j_up = 0.9 * con->P_up * modfac / ( 1. + con->K_up / c.Cai );
    //    j_tr = (c.Ca_NSR-c.Ca_JSR)/con->tho_tr;
    double j_tr = ( c.Ca_NSR - c.Ca_JSR ) / con->tho_tr;
    //    j_uni = con->P_Ca*2*con->phi_m/con->E_T*(con->alpha_m*c.Cam...
    double j_uni = con->P_Ca * 2. * con->phi_m / con->E_T * ( con->alpha_m * c.Cam
                                                             * exp( -2. * con->phi_m / con->E_T ) - con->alpha_e * c.Cai ) / ( exp( -2. * con->phi_m / con->E_T ) - 1. );
    //    j_NaCam = con->Q_mo*(c.Cam/(con->K_Cam+c.Cam));
    double j_NaCam = con->Q_mo * ( c.Cam / ( con->K_Cam + c.Cam ) );
    //
    //    % = Ca2+ buffering
    //    df_TC = con->Fkl*c.Cai*(1-A-TT)-k_l*(A+TT);
    //    df_TMC = con->k_fTMC*c.Cai*(1-f_TMC-f_TMM)-con->k_bTMC*f_TMC;
    //    df_TMM = con->k_fTMM*con->Mgi*(1-f_TMC-f_TMM)-con->k_bTMM*f_TMM;
    //    df_CMi = con->k_fCM*c.Cai*(1-f_CMi)-con->k_bCM*f_CMi;
    //    df_CMs = con->k_fCM*c.Ca_sub*(1-f_CMs)-con->k_bCM*f_CMs;
    //    df_CQ = con->k_fCQ*c.Ca_JSR*(1-f_CQ)-con->k_bCQ*f_CQ;
    double df_TC = con->Fkl * c.Cai * ( 1. - A - TT ) - k_l * ( A + TT );
    double df_TMC = con->k_fTMC * c.Cai * ( 1. - f_TMC - f_TMM ) - con->k_bTMC * f_TMC;
    double df_TMM = con->k_fTMM * con->Mgi * ( 1. - f_TMC - f_TMM ) - con->k_bTMM * f_TMM;
    double df_CMi = con->k_fCM * c.Cai * ( 1. - f_CMi ) - con->k_bCM * f_CMi;
    double df_CMs = con->k_fCM * c.Ca_sub * ( 1. - f_CMs ) - con->k_bCM * f_CMs;
    double df_CQ = con->k_fCQ * c.Ca_JSR * ( 1. - f_CQ ) - con->k_bCQ * f_CQ;
    //
    //    % = dynamics of Ca2+ concentrations in cell compartments
    //    dCai = ((j_Ca_dif*con->V_sub-j_up*con->V_nSR)/con->V_i)-...
    //        (con->CM_tot*df_CMi+con->TC_tot*df_TC)...
    //        -(j_uni-j_NaCam)*con->V_myto/con->V_i;
    //    dCa_sub = j_SRCarel*con->V_jSR/con->V_sub-(I_CaL+I_CaT+I_bCa-2*I_NCX)/...
    //        (2*con->F*con->V_sub)-(j_Ca_dif+con->CM_tot*df_CMs);
    //    dCa_JSR = j_tr-j_SRCarel-con->CQ_tot*df_CQ;
    //    dCa_NSR = j_up-j_tr*con->V_jSR/con->V_nSR;
    //    dCam = (j_uni-j_NaCam)*(con->V_myto/con->V_i);
    double dCai = ( ( j_Ca_dif * con->V_sub - j_up * con->V_nSR ) / con->V_i ) -
    ( con->CM_tot * df_CMi + con->TC_tot * df_TC )
    - ( j_uni - j_NaCam ) * con->V_myto / con->V_i;
    double dCa_sub = j_SRCarel * con->V_jSR / con->V_sub - ( I_CaL + I_CaT + I_bCa - 2. * I_NCX ) /
    ( 2. * con->F * con->V_sub ) - ( j_Ca_dif + con->CM_tot * df_CMs );
    double dCa_JSR = j_tr - j_SRCarel - con->CQ_tot * df_CQ;
    double dCa_NSR = j_up - j_tr * con->V_jSR / con->V_nSR;
    double dCam = ( j_uni - j_NaCam ) * ( con->V_myto / con->V_i );
    //
    //    % = derivative of gating variables
    //    dd_L  = a_l*(1-d_L)-b_l*d_L;
    //    df_L  = a_fl*(1-f_L)-b_fl*f_L;
    //    df_Ca = a_fCa*(1-f_Ca)-b_fCa*f_Ca;
    //    dd_T = a_dT*(1-d_T)-b_dT*d_T;
    //    df_T = a_fT*(1-f_T)-b_fT*f_T;
    //    dp_aF = a_paF*(1-p_aF)-b_paF*p_aF;
    //    dp_aS = a_paS*(1-p_aS)-b_paS*p_aS;
    //    dp_i = a_pi*(1-p_i)-b_pi*p_i;
    //    dn = a_n*(1-n)-b_n*n;
    //    dq = a_q*(1-q)-b_q*q;
    //    dr = a_r*(1-r)-b_r*r;
    //    dy = a_y*(1-y)-b_y*y;
    //    dq_a = a_qa*(1-q_a)-b_qa*q_a;
    //    dq_i = a_qi*(1-q_i)-b_qi*q_i;
    //    dw = a_w*(1-w)-b_w*w;
    double dd_L  = a_l * ( 1. - d_L ) - b_l * d_L;
    double df_L  = a_fl * ( 1. - f_L ) - b_fl * f_L;
    double df_Ca = a_fCa * ( 1. - f_Ca ) - b_fCa * f_Ca;
    double dd_T = a_dT * ( 1. - d_T ) - b_dT * d_T;
    double df_T = a_fT * ( 1. - f_T ) - b_fT * f_T;
    double dp_aF = a_paF * ( 1. - p_aF ) - b_paF * p_aF;
    double dp_aS = a_paS * ( 1. - p_aS ) - b_paS * p_aS;
    double dp_i = a_pi * ( 1. - p_i ) - b_pi * p_i;
    double dn = a_n * ( 1. - n ) - b_n * n;
    double dq = a_q * ( 1. - q ) - b_q * q;
    double dr = a_r * ( 1. - r ) - b_r * r;
    double dy = a_y * ( 1. - y ) - b_y * y;
    double dq_a = a_qa * ( 1. - q_a ) - b_qa * q_a;
    double dq_i = a_qi * ( 1. - q_i ) - b_qi * q_i;
    double dw = a_w * ( 1. - w ) - b_w * w;
    //
    //    % additional outputs that I want (not in coupled equations)
    //    if con->ALL_VAR
    //        dj_SRCarel = j_SRCarel-F(37);
    //        dj_Ca_dif = j_Ca_dif-F(38);
    //        dj_up = j_up-F(39);
    //        dj_tr = j_tr-F(40);
    //        dj_uni = j_uni-F(41);
    //        dj_NaCam = j_NaCam-F(42);
    //
    //        % = return outputs
    //        dF = [dCai; dCa_sub; dCa_JSR; dCa_NSR; dCam
    //            df_TMC; df_TMM; df_CMi; df_CMs; ...
    //            df_CQ; dR; dOO; dS; dRI; dV; dd_L; df_L; ...
    //            df_Ca; dp_aF; dp_aS; dp_i; dn; dy; dd_T; df_T; dq; dr; dq_a; dq_i; ...
    //            dcAMP; dPLB_p; dA; dTT; dU; dSL; dw; ...  % adding force parameters - 36 coupled equations
    //            dj_SRCarel; ... % additional outputs that I want to follow up (not in coupled equations)
    //            dj_Ca_dif; dj_up; dj_tr; dj_uni; dj_NaCam; df_cAB];
    //    else
    //
    //        dF = [dCai; dCa_sub; dCa_JSR; dCa_NSR; dCam
    //            df_TMC; df_TMM; df_CMi; df_CMs; ...
    //            df_CQ; dR; dOO; dS; dRI; dV; dd_L; df_L; ...
    //            df_Ca; dp_aF; dp_aS; dp_i; dn; dy; dd_T; df_T; dq; dr; dq_a; dq_i; ...
    //            dcAMP; dPLB_p; dA; dTT; dU; dSL; ...
    //            dw];
    //    end
    dF[0] = dCai;
    dF[1] = dCa_sub;
    dF[2] = dCa_JSR;
    dF[3] = dCa_NSR;
    dF[4] = dCam;
    dF[5] = df_TMC;
    dF[6] = df_TMM;
    dF[7] = df_CMi;
    dF[8] = df_CMs;
    dF[9] = df_CQ;
    dF[10] = dR;
    dF[11] = dOO;
    dF[12] = dS;
    dF[13] = dRI;
    dF[14] = dV;
    dF[15] = dd_L;
    dF[16] = df_L;
    dF[17] = df_Ca;
    dF[18] = dp_aF;
    dF[19] = dp_aS;
    dF[20] = dp_i;
    dF[21] = dn;
    dF[22] = dy;
    dF[23] = dd_T;
    dF[24] = df_T;
    dF[25] = dq;
    dF[26] = dr;
    dF[27] = dq_a;
    dF[28] = dq_i;
    dF[29] = dcAMP;
    dF[30] = dPLB_p;
    dF[31] = dA;
    dF[32] = dTT;
    dF[33] = dU;
    dF[34] = dSL;
    dF[35] = dw; //...  % adding force parameters - 36 coupled equations
    
    if( con->ALL_VAR ){
        double dj_SRCarel = j_SRCarel - F[36];
        double dj_Ca_dif = j_Ca_dif - F[37];
        double dj_up = j_up - F[38];
        double dj_tr = j_tr - F[39];
        double dj_uni = j_uni - F[40];
        double dj_NaCam = j_NaCam - F[41];
        
        //        % = return outputs
        dF[36] = 0.001 * dj_SRCarel; // ... % additional outputs that I want to follow up (not in coupled equations)
        dF[37] = 0.001 * dj_Ca_dif;
        dF[38] = 0.001 * dj_up;
        dF[39] = 0.001 * dj_tr;
        dF[40] = 0.001 * dj_uni;
        dF[41] = 0.001 * dj_NaCam;
        dF[42] = df_cAB;
        //    else
        //
        //        dF = [dCai; dCa_sub; dCa_JSR; dCa_NSR; dCam
        //            df_TMC; df_TMM; df_CMi; df_CMs; ...
        //            df_CQ; dR; dOO; dS; dRI; dV; dd_L; df_L; ...
        //            df_Ca; dp_aF; dp_aS; dp_i; dn; dy; dd_T; df_T; dq; dr; dq_a; dq_i; ...
        //            dcAMP; dPLB_p; dA; dTT; dU; dSL; ...
        //            dw];
        //    end
    }
    //end
    //
}
