//function y = driving_function_PN_efferent(v, n, p, gVagal, s, w, s_S, s_out, ref, A, B, N)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% This function sets up the driving function for the 2-neuron model
//
//% Inputs:
//% v,n,s,u - state variables
//% p - data structure containing all of the parameters in the system
//% A - Adjacency matrix that determines which neurons are coupled
//% Iapp - applied current
//
//% Outputs:
//% y - vector valued output of driving function for system
//% G - conductances for neuron. (except leakage and "out" (external
//% input)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//

void driving_function_PN_efferent(double *v, double *n,  structp &p, double *gVagal, double *s, double *w, double *s_S, double *s_out, double *ref, structA &A, structB &B, structN &N, double *y);

void driving_function_PN_efferent(double *v, double *n,  structp &p, double *gVagal, double *s, double *w, double *s_S, double *s_out, double *ref, structA &A, structB &B, structN &N, double *y){
    double M;
    M = N.CNS + N.STELLATE + N.CNS;
    using namespace std;
    // cout << s_out[1] << "\t";
    
    double As1, As2, As3, As4, As5, As6, As7, As8;
    double Bv;
    
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        As1 = 0;
        As2 = 0;
        As3 = 0;
        As4 = 0;
        As5 = 0;
        As6 = 0;
        As7 = 0;
        As8 = 0;
        Bv = 0;
        
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            As1 += A.PN_EE[id2][id1] * s[id2];
            As2 += A.PN_EI[id2][id1] * s[id2];
            As3 += A.SICNS_to_PN_EE[id2][id1] * s_S[N.CNS + N.STELLATE + id2];
            As4 += A.SICNS_to_PN_EI[id2][id1] * s_S[N.CNS + N.STELLATE + id2];
            As5 += A.PN_IE[id2][id1] * s[id2];
            As6 += A.PN_II[id2][id1] * s[id2];
            As7 += A.SICNS_to_PN_IE[id2][id1] * s_S[N.CNS + N.STELLATE + id2];
            As8 += A.SICNS_to_PN_II[id2][id1] * s_S[N.CNS + N.STELLATE + id2];
            Bv += B.PN[id2][id1] * ( v[id2] - p.ESyn_E );
        }
        
        //y = (-p.gL * (v - p.EL) ... % leakage
        y[ id1 ] = ( -p.gL * ( v[ id1 ] - p.EL ) // leakage
                    //    - p.gK * (n.^2) .* (v - p.EK) ... % delayed rectifier
                    - p.gK * (n[id1] * n[id1]) * (v[id1] - p.EK) // delayed rectifier
                    //    - p.gSyn_E * ref .*  (A.PN_EE' * s + A.PN_EI' * s + A.SICNS_to_PN_EE' * s_S(N.CNS + N.STELLATE + 1: M) + A.SICNS_to_PN_EI' * s_S(N.CNS + N.STELLATE + 1: M)) .*  (v - p.ESyn_E)... % excitatory connections
                    - p.gSyn_E * ref[id1] *  ( As1 + As2 + As3 + As4) *  (v[id1] - p.ESyn_E)  //% excitatory connections
                    //    - p.gSyn_I * ref .* (A.PN_IE' * s + A.PN_II' * s + A.SICNS_to_PN_IE' * s_S(N.CNS + N.STELLATE + 1: M) + A.SICNS_to_PN_II' * s_S(N.CNS + N.STELLATE + 1: M)) .* (v - p.ESyn_I) ... % inhibitory connections
                    - p.gSyn_I * ref[id1] * (As5 + As6 + As7 + As8) * (v[id1] - p.ESyn_I) //% inhibitory connections
                    //    - gVagal .* ref .* B.PN' *(v - p.ESyn_E)... % efferent
                    - gVagal[id1] * ref[id1] * Bv // efferent
                    //    - p.g_out .* ref .* s_out .* (v - p.ESyn_E)... % cardiac
                    - p.g_out * ref[id1] * s_out[id1] * (v[id1] - p.ESyn_E) // cardiac
                    //    - (p.gM_PN .* w) .* (v - p.EK))/p.C; % M-current
                    - (p.gM_PN[id1] * w[id1]) * (v[id1] - p.EK) ) / p.C ; // M-current
    }
}
