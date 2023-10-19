void update_I_f( cell_con * con, double V_m, cell_ac * ac, double y, double *p_I_f, double *p_a_y, double *p_b_y );

void update_I_f( cell_con * con, double V_m, cell_ac * ac, double y, double *p_I_f, double *p_a_y, double *p_b_y ){

    double V_shift = ( con->K_if
                      * pow( ac->cAMP , con->n_if )
                      / ( pow( con->K_05if , con->n_if ) + pow( ac->cAMP , con->n_if ) )
                      - 18.76 );

    double V_If12 = V_shift + con->V_If12;

    double y_inf = 1. / ( 1. + exp( ( V_m - V_If12 ) / 13.5 ) );
    double tho_y = 0.7166529 / ( exp( -( V_m + 386.9 ) / 45.302 ) + exp( ( V_m - 73.08 ) / 19.231 ) );

    double I_fNa = con->C * 0.3833 * con->g_If * ( V_m - con->E_Na ) * y * y;
    double I_fK = con->C * 0.6167 * con->g_If * ( V_m - con->E_K ) * y * y;

    double a_y = y_inf / tho_y;
    double b_y = ( 1. - y_inf ) / tho_y;

    double I_f = I_fNa + I_fK;
 
    *p_I_f = I_f;
    *p_a_y = a_y;
    *p_b_y = b_y;
    
}
