void update_I_Kr( cell_con * con, double V_m, double p_aF, double p_aS, double p_i, double *p_I_Kr, double *p_a_paF, double *p_b_paF, double *p_a_paS, double *p_b_paS, double *p_a_pi, double *p_b_pi );

void update_I_Kr( cell_con * con, double V_m, double p_aF, double p_aS, double p_i, double *p_I_Kr, double *p_a_paF, double *p_b_paF, double *p_a_paS, double *p_b_paS, double *p_a_pi, double *p_b_pi ){

    
    double I_Kr = con->C * con->g_Kr * ( V_m - con->E_K ) * ( 0.6 * p_aF + 0.4 * p_aS ) * p_i;

    double pa_inf = 1. / ( 1. + exp( -( V_m + 23.2 ) / 10.6 ) );
    double pi_inf = 1. / ( 1. + exp( ( V_m + 28.6 ) / 17.1 ) );

    double tho_paF = 0.84655354 / ( 0.0372 * exp( V_m / 15.9 ) + 0.00096 * exp( -V_m / 22.5 ) );
    double tho_paS = 0.84655354 / ( 0.0042 * exp( V_m / 17.0 ) + 0.00015 * exp( -V_m / 21.6 ) );
    double tho_pi = 1. / ( 0.1 * exp( -V_m / 54.645 ) + 0.656 * exp( V_m / 106.157 ) );

    double a_paF = pa_inf / tho_paF;
    double b_paF = ( 1. - pa_inf ) / tho_paF;
    double a_paS = pa_inf / tho_paS;
    double b_paS = ( 1. - pa_inf ) / tho_paS;
    double a_pi = pi_inf / tho_pi;
    double b_pi = ( 1. - pi_inf ) / tho_pi;

    *p_I_Kr = I_Kr;
    *p_a_paF = a_paF;
    *p_b_paF = b_paF;
    *p_a_paS = a_paS;
    *p_b_paS = b_paS;
    *p_a_pi = a_pi;
    *p_b_pi = b_pi;
    
}
