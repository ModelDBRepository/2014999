void update_I_bNa( cell_con * con, double V_m, double *p_I_bNa );

void update_I_bNa( cell_con * con, double V_m, double *p_I_bNa ){

    *p_I_bNa = con->C * con->g_bNa * ( V_m - con->E_Na );
}
