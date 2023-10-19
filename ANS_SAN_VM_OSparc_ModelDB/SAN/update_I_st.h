void update_I_st( cell_con * con, double V_m, double qa, double qi, double *p_I_st, double *p_a_qa, double *p_b_qa, double *p_a_qi, double *p_b_qi );

void update_I_st( cell_con * con, double V_m, double qa, double qi, double *p_I_st, double *p_a_qa, double *p_b_qa, double *p_a_qi, double *p_b_qi ){

    double I_st = con->C * con->g_st * ( V_m - con->E_st ) * qa * qi;

    double qa_inf = 1. / ( 1. + exp( -( V_m + 57. ) / 5. ) );
    double alpha_qa = 1. / ( 0.15 * exp( -V_m / 11. ) + 0.2 * exp( -V_m / 700. ) );
    double beta_qa = 1. / ( 16. * exp( V_m / 8. ) + 15. * exp( V_m / 50. ) );
    double tho_qa = 1. / ( alpha_qa + beta_qa );

    double alpha_qi = 1. / ( 3100. * exp( V_m / 13. ) + 700. * exp( V_m / 70. ) );
    double beta_qi = 1. / ( 95. * exp( -V_m / 10. ) + 50. * exp( -V_m / 700. ) ) + 0.000229 / ( 1. + exp( -V_m / 5. ) );
    double tho_qi = 6.65 / ( alpha_qi + beta_qi );
    double qi_inf = alpha_qi / ( alpha_qi + beta_qi );

    double a_qa = qa_inf / tho_qa;
    double b_qa = ( 1. - qa_inf ) / tho_qa;
    double a_qi = qi_inf / tho_qi;
    double b_qi = ( 1. - qi_inf ) / tho_qi;

    *p_I_st = I_st;
    *p_a_qa = a_qa;
    *p_b_qa = b_qa;
    *p_a_qi = a_qi;
    *p_b_qi = b_qi;
    
}
