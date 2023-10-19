void update_I_KACh( cell_con * con, double V_m, double w, double *p_I_KACh, double *p_a_w, double *p_b_w );

void update_I_KACh( cell_con * con, double V_m, double w, double *p_I_KACh, double *p_a_w, double *p_b_w ){

//function [I_KACh,a_w,b_w] = update_I_KACh(con,V_m,w)

    double I_KACh = con->C * w * con->g_KACh_max * ( V_m - con->E_K );
    double beta_w = 0.001 * 12.32 / ( 1. + 0.0042 / ( con->CCh * 1E-6 ) ); // % CCh in mM
    double alpha_w = 0.001 * 17. * exp( 0.0133 * ( V_m + 40. ) );
    double w_inf = beta_w / ( alpha_w + beta_w );
    double tho_w = 1. / ( alpha_w + beta_w );

    double a_w = w_inf / tho_w;
    double b_w = ( 1. - w_inf ) / tho_w;

    *p_I_KACh = I_KACh;
    *p_a_w = a_w;
    *p_b_w = b_w;
    
}
