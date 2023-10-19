void update_I_NaK( cell_con * con, double V_m, cell_ac * ac, double *p_I_NaK );

void update_I_NaK( cell_con * con, double V_m, cell_ac * ac, double *p_I_NaK ){

//%INAK_ATP = 1.005/(1+con.k_NaK_ATP/ac.ATP*(1+(con.CATPi-ac.ATP)/con.k_NaK_ADP));
    double INAK_ATP = 1.;
    *p_I_NaK = (con->C
                * con->I_NaKmax
                * INAK_ATP
                / ( ( 1. + pow( ( con->K_mKp / con->Ko ) , 1.2 ) )
                   * ( 1. + pow( ( con->K_mNap / con->Nai ) , 1.3 ) )
                   * ( 1. + exp( -( V_m - con->E_Na + 120. ) / 30. ) ) ) );

}
