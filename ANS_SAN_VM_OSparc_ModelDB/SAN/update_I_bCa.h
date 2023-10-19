void update_I_bCa( cell_con * con, double V_m, double *p_I_bCa );

void update_I_bCa( cell_con * con, double V_m, double *p_I_bCa ){

 *p_I_bCa = con->C * con->g_bCa * ( V_m - con->E_CaL );

}
