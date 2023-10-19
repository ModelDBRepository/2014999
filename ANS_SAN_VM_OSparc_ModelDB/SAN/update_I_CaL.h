
void update_I_CaL( cell_con * con, double V_m, cell_c * c, cell_ac * ac, double d_L, double f_L, double f_Ca, double *p_I_CaL, double *p_a_l, double *p_b_l, double *p_a_fl, double *p_b_fl, double *p_a_fCa, double *p_b_fCa );

void update_I_CaL( cell_con * con, double V_m, cell_c * c, cell_ac * ac, double d_L, double f_L, double f_Ca, double *p_I_CaL, double *p_a_l, double *p_b_l, double *p_a_fl, double *p_b_fl, double *p_a_fCa, double *p_b_fCa ){

    double I_CaL, a_l, b_l, a_fl, b_fl, a_fCa, b_fCa;
//% -20% decrease
    double b_CaL = -0.2152 + 1.6913 * pow( ac->PKA , 10.0808 ) / ( pow( 0.8836 , 10.0808 ) + pow( ac->PKA , 10.0808 ) );

    I_CaL = con->C * con->g_CaL * ( b_CaL + 1. ) * ( V_m - con->E_CaL ) * d_L * f_L * f_Ca;

    double alpha_dL = ( -0.02839 * ( V_m + 35. ) / ( exp( -( V_m + 35. ) / 2.5 ) - 1. )
                       - 0.0849 * V_m / ( exp( -V_m / 4.8 ) - 1. ) );
    
    double beta_dL;
    // if( V_m != 5. ) {
    if( fabs( V_m - 5. ) > 1E-8 ) {
        beta_dL = 0.01143 * ( V_m - 5. ) / ( exp( ( V_m - 5. ) / 2.5 ) - 1. );
    } else {
        beta_dL = 0.02858;
    }
    double d_Linf = 1. / ( 1. + exp( -( V_m + 13.5) / 6. ) );
    double f_Linf = 1. / ( 1. + exp( ( V_m + 35. ) / 7.3 ) );
   

    
    double f_Cainf = con->K_mfCa / ( con->K_mfCa + c->Ca_sub );
    double tho_dL = 1. / ( alpha_dL + beta_dL );
    double tho_fL = 257.1 * exp( - ( ( V_m + 32.5 ) / 13.9 ) * ( ( V_m + 32.5 ) / 13.9 ) ) + 44.3;
    double tho_fCa = f_Cainf / con->alpha_fCa;

    a_l = d_Linf / tho_dL;
    b_l = ( 1. - d_Linf ) / tho_dL;
    a_fl = f_Linf / tho_fL;
    b_fl = ( 1. - f_Linf ) / tho_fL;
    a_fCa = f_Cainf / tho_fCa;
    b_fCa = (1. - f_Cainf ) / tho_fCa;


    *p_I_CaL = I_CaL;
    *p_a_l = a_l;
    *p_b_l = b_l;
    *p_a_fl = a_fl;
    *p_b_fl = b_fl;
    *p_a_fCa = a_fCa;
    *p_b_fCa = b_fCa;

}










