void update_I_Ks( cell_con * con, double V_m, cell_ac * ac, double n, double *p_I_Ks, double *p_a_n, double *p_b_n );

void update_I_Ks( cell_con * con, double V_m, cell_ac * ac, double n, double *p_I_Ks, double *p_a_n, double *p_b_n ){

    double I_Ks = con->C * con->g_Ks * ( V_m - con->E_Ks ) * n * n;
    double V_shift = 0.; // % FIXME: for now I put that to zero

    double alpha_n = 0.014 / ( 1. + exp( -( V_m + V_shift - 40. ) / 9. ) );
    double beta_n = 0.001 * exp( -( V_m + V_shift ) / 45. );
    double n_inf = alpha_n / ( alpha_n + beta_n );
    double tho_n = 1. / ( alpha_n + beta_n );

    double a_n = n_inf / tho_n;
    double b_n = ( 1. - n_inf ) / tho_n;
    
    *p_I_Ks = I_Ks;
    *p_a_n = a_n;
    *p_b_n = b_n;
    
}
