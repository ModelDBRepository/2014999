void update_I_CaT( cell_con * con, double V_m, double d_T, double f_T, double *p_I_CaT, double *p_a_dT, double *p_b_dT, double *p_a_fT, double *p_b_fT );

void update_I_CaT( cell_con * con, double V_m, double d_T, double f_T, double *p_I_CaT, double *p_a_dT, double *p_b_dT, double *p_a_fT, double *p_b_fT ){

    double I_CaT, a_dT, b_dT, a_fT, b_fT;
    
    I_CaT = con->C * con->g_CaT * ( V_m - con->E_CaT ) * d_T * f_T;

    double d_Tinf = 1. / ( 1. + exp( -( V_m + 26.3 ) / 6. ) );
    double f_Tinf = 1. / ( 1. + exp( ( V_m + 61.7 ) / 5.6 ) );
    double tho_dT = 1. / ( 1.068 * exp( ( V_m + 26.3 ) / 30. ) + 1.068 * exp( -( V_m + 26.3 ) / 30. ) );
    double tho_fT = 1. / ( 0.0153 * exp( -( V_m + 61.7 ) / 83.3 ) + 0.015 * exp( ( V_m + 61.7 ) / 15.38 ) );

    a_dT = d_Tinf / tho_dT;
    b_dT = ( 1. - d_Tinf ) / tho_dT;
    a_fT = f_Tinf / tho_fT;
    b_fT = ( 1. - f_Tinf ) / tho_fT;
    
    *p_I_CaT = I_CaT;
    *p_a_dT = a_dT;
    *p_b_dT = b_dT;
    *p_a_fT = a_fT;
    *p_b_fT = b_fT;

}
