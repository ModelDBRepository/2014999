void update_I_4AP( cell_con * con, double V_m, double q, double r, double *p_I_to, double *p_I_sus, double *p_a_q, double *p_b_q, double *p_a_r, double *p_b_r );

void update_I_4AP( cell_con * con, double V_m, double q, double r, double *p_I_to, double *p_I_sus, double *p_a_q, double *p_b_q, double *p_a_r, double *p_b_r ){

    double I_to = con->C * con->g_to * ( V_m - con->E_K ) * q * r;
    double I_sus = con->C * con->g_sus * ( V_m - con->E_K ) * r;

    double q_inf = 1. / ( 1. + exp( ( V_m + 49. ) / 13. ) );
    double r_inf = 1. / ( 1. + exp( -( V_m - 19.3 ) / 15. ) );
    double tho_q = 39.102 / ( 0.57 * exp( -0.08 * ( V_m + 44. ) )
                             + 0.065 * exp( 0.1 * ( V_m + 45.93 ) ) ) + 6.06;
    double tho_r = 14.40516 / ( 1.037 * exp( 0.09 * ( V_m + 30.61 ) )
                               + 0.369 * exp( -0.12 * ( V_m + 23.84 ) ) ) + 2.75352;

    double a_q = q_inf / tho_q;
    double b_q = ( 1. - q_inf ) / tho_q;
    double a_r = r_inf / tho_r;
    double b_r = ( 1. - r_inf ) / tho_r;

    *p_I_to = I_to;
    *p_I_sus = I_sus;
    *p_a_q = a_q;
    *p_b_q = b_q;
    *p_a_r = a_r;
    *p_b_r = b_r;
    
}
