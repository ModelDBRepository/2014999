void update_I_NCX( cell_con * con, double V_m, cell_c * c, double *p_I_NCX );

void update_I_NCX( cell_con * con, double V_m, cell_c * c, double *p_I_NCX ){

    double do1 = ( 1.
                 + ( con->Cao / con->K_co ) * ( 1. + exp( con->Q_co * V_m / con->E_T ) )
                 + ( con->Nao / con->K_1no ) * ( 1. + ( con->Nao / con->K_2no )
                                                * ( 1. + con->Nao / con->K_3no ) ) );

    double k_43 = con->Nai / ( con->K_3ni + con->Nai );
    
    double k_41 = exp( -con->Q_n * V_m / ( 2. * con->E_T ) );
    
    double k_34 = con->Nao / ( con->K_3no + con->Nao );
    
    double k_21 = ( con->Cao / con->K_co ) * exp( con->Q_co * V_m / con->E_T ) / do1;
    
    double k_23 = ( ( con->Nao / con->K_1no )
                   * ( con->Nao / con->K_2no )
                   * ( 1. + con->Nao / con->K_3no )
                   * exp( -con->Q_n * V_m / ( 2. * con->E_T ) )
                   / do1 );
    
    double k_32 = exp( con->Q_n * V_m / ( 2. * con->E_T ) );

    double x_1 = k_34 * k_41 * ( k_23 + k_21 ) + k_21 * k_32 * ( k_43 + k_41 );

    double di = ( 1.
                 + ( c->Ca_sub / con->K_ci )
                 * ( 1. + exp( -con->Q_ci * V_m / con->E_T ) + con->Nai / con->K_cni )
                 + ( con->Nai / con->K_1ni )
                 * ( 1. + ( con->Nai / con->K_2ni ) * ( 1. + con->Nai / con->K_3ni ) ) );

    double k_12 = ( c->Ca_sub / con->K_ci ) * exp( -con->Q_ci * V_m / con->E_T ) / di;
    
    double k_14 = ( ( con->Nai / con->K_1ni )
                   * ( con->Nai / con->K_2ni )
                   * ( 1. + con->Nai / con->K_3ni )
                   * exp( con->Q_n* V_m / ( 2. * con->E_T ) )
                   / di );

    double x_2 = k_43 * k_32 * ( k_14 + k_12 ) + k_41 * k_12 * ( k_34 + k_32 );
    double x_3 = k_43 * k_14 * ( k_23 + k_21 ) + k_12 * k_23 * ( k_43 + k_41 );
    double x_4 = k_34 * k_23 * ( k_14 + k_12 ) + k_21 * k_14 * ( k_34 + k_32 );

    *p_I_NCX = con->C * con->k_NCX * ( k_21 * x_2 - k_12 * x_1 ) / ( x_1 + x_2 + x_3 + x_4 );

    
}
