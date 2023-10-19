//    function y = driving_function_stochastic_input(v, n, s, w, p, W, s_out, x_PN, I_inj,ref, A, A_1, N)
//        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        % This function sets up the driving function for the 2-neuron model
//        
//        % Inputs:
//        % v,n,s,u - state variables
//        % p - data structure containing all of the parameters in the system
//        % A - Adjacency matrix that determines which neurons are coupled
//        % Iapp - applied current
//        
//        % Outputs:
//        % y - vector valued output of driving function for system
//        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
void driving_function_stochastic_input( double *v, double *n, double *s, double *w, structp &p, structW &W, double *s_out, double *x_PN, double *I_inj, double *ref, structA &A, structM3PN &A_1, structN &N, double *y);

void driving_function_stochastic_input( double *v, double *n, double *s, double *w, structp &p, structW &W, double *s_out, double *x_PN, double *I_inj, double *ref, structA &A, structM3PN &A_1, structN &N, double *y) {
    using namespace std;
    
    int M = N.CNS + N.STELLATE + N.PN;
    double s_PN_E[ M ], s_PN_I[ M ], s1_out[ M ];
    double sum1;
    
//          s_PN_E = [zeros(N.CNS + N.STELLATE,1); A.PN_to_SICNS_EE' * x_PN(1:N.PN) + A.PN_to_SICNS_EI' * x_PN(1:N.PN)];
    for( int id1 = 0; id1 < (N.CNS+N.STELLATE); id1++ ) {
        s_PN_E[ id1 ] = 0;
    }
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        sum1 = 0;
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            sum1 += ( A.PN_to_SICNS_EE[id2][id1] + A.PN_to_SICNS_EI[id2][id1] ) * x_PN[ id2 ] ;
        }
        s_PN_E[ id1 + N.CNS+N.STELLATE ] = sum1;
    }
//          s_PN_I = [zeros(N.CNS + N.STELLATE,1); A.PN_to_SICNS_IE' * x_PN(1:N.PN) + A.PN_to_SICNS_II' * x_PN(1:N.PN)];
    for( int id1 = 0; id1 < (N.CNS+N.STELLATE); id1++ ) {
        s_PN_I[ id1 ] = 0;
    }
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        sum1 = 0;
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            sum1 += ( A.PN_to_SICNS_IE[id2][id1] + A.PN_to_SICNS_II[id2][id1] ) * x_PN[ id2 ] ;
            //cout << A.PN_to_SICNS_IE[id2][id1] << "  ";
        }
        //cout << endl;
        
        s_PN_I[ id1 + N.CNS+N.STELLATE ] = sum1;
    }

//          s_out = [zeros(N.CNS,1);  A_1' * s_out(1:N.CNS); s_out(N.CNS + N.STELLATE + 1 : end)];
    for( int id1 = 0; id1 < N.CNS; id1++ ) {
        sum1 = 0;
        for( int id2 = 0; id2 < N.CNS; id2++ ) {
            sum1 += A_1.M[id2][id1] * s_out[id2];
        }
        s1_out[ id1 + N.CNS ] = sum1;
    }

    for( int id1 = 0; id1 < N.CNS; id1++ ) {
        s1_out[ id1 ] = 0;
    }
    for( int id1 = (N.CNS+N.STELLATE); id1 < M; id1++ ) {
        s1_out[ id1 ] = s_out[ id1 ];
    }
//
//        y = (-p.gL * (v - p.EL) ...
//            - p.gK * (n.^2) .* (v - p.EK)...
//           - (p.gM_S .* w) .* (v - p.EK) ...
//           - p.gSyn_E * ref .* ((W.EE' + W.EI') * s + s_PN_E) .* (v - p.ESyn_E) ...
//           - p.gSyn_I * ref .* ((W.IE' + W.II') * s  + s_PN_I) .* (v - p.ESyn_I) ...
//           - ref .* (p.g_out .* s_out) .* (v - p.E_out)...
//           + ref .* I_inj)/p.C;
//    end
    double WEEWEIs[M], WIEWIIs[M];
    for( int id1 = 0; id1 < M; id1++ ) {
        WEEWEIs[id1] = 0;
        WIEWIIs[id1] = 0;
        
        for( int id2 = 0; id2 < M; id2++ ) {
            WEEWEIs[id1] += ( W.EE[id2][id1] + W.EI[id2][id1] ) * s[id2];
            WIEWIIs[id1] += ( W.IE[id2][id1] + W.II[id2][id1] ) * s[id2];
        }
    }
    //cout << "s_PN_I  " ;
    for( int id1 = 0; id1 < M; id1++ ) {
        y[ id1 ] = ( -p.gL * ( v[ id1 ] - p.EL )
                    - p.gK * n[id1] * n[id1] * ( v[id1] - p.EK )
                    - (p.gM_S[id1] * w[id1] ) * ( v[id1] - p.EK )
                    - p.gSyn_E * ref[id1] * ( WEEWEIs[id1] + s_PN_E[id1] ) * ( v[id1] - p.ESyn_E )
                    - p.gSyn_I * ref[id1] * ( WIEWIIs[id1] + s_PN_I[id1] ) * ( v[id1] - p.ESyn_I )
                    - ref[id1] * ( p.g_out * s1_out[ id1 ] ) * ( v[id1] - p.E_out )
                    + ref[id1] * I_inj[id1] 
                    ) / p.C;
        //cout << s_PN_I[id1] << "  ";
    }
    //cout << endl;
}
