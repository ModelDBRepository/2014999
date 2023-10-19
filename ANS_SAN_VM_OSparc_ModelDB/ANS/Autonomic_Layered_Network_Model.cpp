#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <random>
#include <new>


#include <sys/types.h>
#include <sys/stat.h>

const int size_PN = 50;

typedef struct {
    double ICNS; //: 0.3333
    double STELLATE; //: 0.3333
    double CNS; //: 0.3333
    double PN; //: 0.3333
    double efferent_STELLATE; //: 0.2500
    double efferent_ICNS; //: 0.2500
    double afferent_STELLATE;
    double S_to_P; //: 0.2000
    double P_to_S; //: 0.2000
    double v_reset; //: -68
    double v_T; //: -60
    double EL; //: -70
    double EK; //: -90
    double ESyn_E; //: 0
    double ESyn_I; //: -80
    double E_out; //: 0
    double nss; //: 0
    double tau_n; //: 75
    double tau_w; //: 300
    double dn; //: 0.2500
    double gL; //: 0.1000
    double gK; //: 0.1000
    double g_out; //: 5.0000e-04
    double C; //: 0.7727
    double gM_S[ size_PN * 3 ]; //: [60×1 double]
    double gM_PN[ size_PN ]; //: [20×1 double]
    double gSyn_E; //: 5.0000e-04
    double gSyn_I; //: 5.0000e-04
    double tau_w_reset;
} structp;

typedef struct{
    double ds; //: 0.1500
    double lambda1; //: 0.0250
    double lambda2; //: 0.0100
    double dISO;
    double dACh;
} structparam;

typedef struct{
    int ICNS; //: 20
    int STELLATE; //: 20
    int CNS; //: 20
    int PN; //: 20
} structN;

typedef struct{
    double PN[size_PN][size_PN]; //: [20×20 double]
    double PN_EE[size_PN][size_PN]; //: [20×20 double]
    double PN_II[size_PN][size_PN]; //: [20×20 double]
    double PN_EI[size_PN][size_PN]; //: [20×20 double]
    double PN_IE[size_PN][size_PN]; //: [20×20 double]
    
    double SICNS_to_PN[size_PN][size_PN]; //: [20×20 double]
    double SICNS_to_PN_EE[size_PN][size_PN]; //: [20×20 double]
    double SICNS_to_PN_II[size_PN][size_PN]; //: [20×20 double]
    double SICNS_to_PN_EI[size_PN][size_PN]; //: [20×20 double]
    double SICNS_to_PN_IE[size_PN][size_PN]; //: [20×20 double]
    
    double PN_to_SICNS[size_PN][size_PN]; //: [20×20 double]
    double PN_to_SICNS_EE[size_PN][size_PN]; //: [20×20 double]
    double PN_to_SICNS_II[size_PN][size_PN]; //: [20×20 double]
    double PN_to_SICNS_EI[size_PN][size_PN]; //: [20×20 double]
    double PN_to_SICNS_IE[size_PN][size_PN]; //: [20×20 double]
    
    double SICNS_to_PN_E[size_PN][size_PN]; //: [20×20 double]
    double SICNS_to_PN_I[size_PN][size_PN]; //: [20×20 double]
} structA;

typedef struct{
    double M[size_PN*3][size_PN*3];
} structM3PN;

typedef struct{
    double PN[size_PN][size_PN]; //: [20×20 double]
} structB;

typedef struct {
    double EE; //: 1
    double EI; //: 5
    double IE; //: 300
    double II; //: 10
} structdistribution;

typedef struct{
    double EE[size_PN*3][size_PN*3]; //: [60×60 double]
    double EI[size_PN*3][size_PN*3]; //: [60×60 double]
    double IE[size_PN*3][size_PN*3]; //: [60×60 double]
    double II[size_PN*3][size_PN*3]; //: [60×60 double]
} structW;

typedef struct{
    double strong; //: 5
    double weak; //: 0.5000
    double PN; //: 1
} structweights;

typedef struct{
    double SICNS; // = 0.1;
    double STELLATE; // = 0.5;
    double PICNS; // = 0;
} structtone;

typedef struct {
    double max; // = 0.8;
    double min; // = 0.2;
} structP_ICNS;

typedef struct {
    double max; // = 0.8;
    double min; // = 0.4;
} structP_STELLATE;

// w_ss = @(v) 1 ./(1 + exp(-(v + 45)/2.4));
double w_ss( double v ) {
    return ( 1. / ( 1. + exp( - ( v + 45. ) / 2.4 ) ) );
}

double intensity_function( double t, double p1, double p2, double p_drop ) {
    double t1 = fmod( t, 1000. );
    double intensity = ( ( p2 / 50. * t1 + p1 ) * (double)( t1 <= 50 )
                        + ( ( p_drop - 1. ) * ( p1 + p2 ) / 200. * ( t1 - 50 ) + p1 + p2 )
                        * (double)( 50. < t1 && t1 <= 250 )
                        + ( ( p1 - p_drop * ( p1 + p2 ) ) / 250. * ( t1 - 500 ) + p1 )
                        * (double)( 250 < t1 && t1 <= 500 )
                        + p1 * (double)( 500 < t1 )
                        );
    return intensity;
}

void intensity_function_vector( double t, structtone &tone, structP_ICNS &P_ICNS, structP_STELLATE &P_STELLATE, double alpha, structN &N, int K, double * y );

#include "Erdos_Renyi.h"
#include "sympathetic_neural_network.h"
#include "couple_two_networks.h"

#include "s_update_stochastic_input.h"
#include "u_update_stochastic_input.h"
#include "driving_function_PN_efferent.h"
#include "driving_function_stochastic_input.h"

int main( int argc, char *argv[] ) {
    using namespace std;
    cout.precision(6);
    
    std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> urand(0.0,1.0);
    
    std::gamma_distribution<double> gammarand1(0.1,1.0);
    std::gamma_distribution<double> gammarand2(0.2,1.0);
    
    std::default_random_engine generatorg1(floor( urand(generator) * 1E10 ));
    std::default_random_engine generatorg2(floor( urand(generator) * 1E10 ));
    
    char *outputFolder = argv[2];
    mkdir( outputFolder, 0777 );

    char filename[90];
    sprintf( filename, "%s/%s", outputFolder, "outputfile.txt" );
    
    FILE *output;
    output = fopen( filename, "w" );
    
    char *InputFile1 = argv[1];
    FILE *fp1;
    fp1 = fopen( InputFile1, "r" );
    
    if ( fp1 == NULL ) {
        
        printf( "%s%s%s", "\nFile \'", InputFile1 , "\' NOT found. \n \n" );
        
        exit(1);
    }
    
    cout << "stimulus Parameters: " << endl;
    
    double input1[1];
    for ( int idy = 0; idy < 1; idy++ ) {
        fscanf ( fp1, "%lf", &input1[idy] );
        cout << idy << "\t" << input1[idy] << endl;
    }
    fclose( fp1 );
    
    // Simulation time
    double T = input1[0];  // end time in sec
    T = T * 1000; // end time in ms
    double dt = 1. / 128;
    // Option to set up new network, E and I markers,synaptic weights, external input receivers.
    // This is set up to allow the rerunning of the simulation with the same
    // network, inputs, E/I distribution, and weights.  If a new network is made
    // with different size, you must re-generate all of the rest of the
    // parameters/inputs.
    
    char new_neural_connections[] = "yes";
    char set_up_new_network[5], set_up_new_E_I_dist[5], set_up_new_weights[5];
    char set_up_new_M_current_dist[5];
    
    if( strcmp(new_neural_connections, "yes") == 0 ) {
        // generate a new configuration model random network?
        strcpy( set_up_new_network, "yes" );
        
        // generate new index of excit. and inhib. neurons? (must say yes if new network is generated)
        strcpy( set_up_new_E_I_dist , "yes" );
        
        // set up new weights for network? (must say yes if new network is generated)
        strcpy( set_up_new_weights , "yes" );
        
        // set up a new distribution for the M-current?
        strcpy( set_up_new_M_current_dist , "yes" );
        
    } else {
        // generate a new configuration model random network?
        strcpy( set_up_new_network , "no" );
        
        // generate new index of excit. and inhib. neurons? (must say yes if new network is generated)
        strcpy( set_up_new_E_I_dist , "no" );
        
        // set up new weights for network? (must say yes if new network is generated)
        strcpy( set_up_new_weights , "no" );
        
        // set up a new distribution for the M-current?
        strcpy( set_up_new_M_current_dist , "no" );
        
    }
    
    // Set up adjacency matrix that describes connections between neurons
    structN N;
    int M;
    structparam param;
    double pulse_length;
    double t;
    
    
    // number of neurons
    N.ICNS = size_PN;
    N.STELLATE = N.ICNS;
    N.CNS = N.ICNS; // set to 0
    N.PN = N.ICNS;
    M = N.ICNS + N.STELLATE + N.CNS;
    
    
    // probability of connections
    structp p;
    p.ICNS = log( N.ICNS ) / ( N.ICNS - 1. ); // prob. of intra-layer connection in SICNS
    p.STELLATE = log( N.STELLATE ) / ( N.STELLATE - 1. ); // prob. of intra-layer connection in STELLATE
    // Set p.CNS = 0 to agree witn CNS for parasympathetic
    p.CNS = log( N.CNS ) / ( N.CNS - 1. ); // 1/(N.CNS + 1); % prob. of intra-layer connection in CNS - SET TO 0
    p.PN = log( N.PN ) / ( N.PN - 1. ); // prob. of intra-layer connection in PN
    p.efferent_STELLATE = 1.; // prob. of efferent input to STELLATE from CNS
    p.afferent_STELLATE = ( p.STELLATE + p.ICNS ) / 4.; // prob of afferent input to STELLATE from ICNS
    // p.ext_SICNS = 1/2; % prob of external (cardiac/CNS) input to SICNS
    // p.ext_PN = 1/2; % prob of cardiac input to SICNS
    p.efferent_ICNS = 1. / 2. ; // prob. of efferent input to ICNS from STELLATE
    p.S_to_P = 1. / ( N.ICNS + N.PN  + 1. ); // p.ICNS / 4; % prob. of SICNS to PN connections
    p.P_to_S = p.S_to_P; // p.PN / 4; % prob. of PN to SICNS connections
    
    // number of (weak) efferent connections in 1 + n protocol between adjacent layers
    int n_S_int = 3;
    
    // weights of 1 + n connections
    structweights weights;
    weights.strong = 1;
    weights.weak = 1. / n_S_int;
    
    
    // generate adjacency matrix for sympathetic network (A)
    // [A_symp, A_CNS_to_STELLATE] = sympathetic_neural_network(N, weights, n_S, p);
    // double A_symp[size_PN*3][size_PN*3], A_CNS_to_STELLATE[size_PN*3][size_PN*3];
    structM3PN A_symp, A_CNS_to_STELLATE;
    for( int id1=0; id1< M; id1++ ) {
        for( int id2 = 0; id2 < M; id2++ ) {
            A_symp.M[id1][id2] = 0;
            A_CNS_to_STELLATE.M[id1][id2] = 0;
        }
    }
    sympathetic_neural_network( N, weights, n_S_int, p, A_symp, A_CNS_to_STELLATE , generator);
    
   
    
    // weights of vagal input to PN
    weights.PN = 1.;
    
    // build Erdos-Renyi random graph
    // A.PN = Erdos_Renyi(N.PN,p.PN);
    // B.PN = diag(ones(N.PN,1));
    structA A;
    double **APN;
    APN = new double *[N.PN];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        APN[id1] = A.PN[id1];
    }
    Erdos_Renyi( N.PN, p.PN, APN , generator );
    
    for( int id1=0; id1< N.PN; id1++ ) {
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            
        }
       
    }
    
    structB B;
    for( int idi = 0; idi < N.PN; idi++ ) {
        for( int idj = 0; idj < N.PN; idj++ ) {
            B.PN[idi][idj] = 0.;
        }
        B.PN[idi][idi] = 1.;
    }
    
    // % generate weighted adjacency matrix for network (A) and adjacency matrix for
    // % efferent input (B) where the 0.25 and 0.5 are parameters that control
    // % the weights of the connections.
    // % A.PN = A.PN .* (1 - 0.25 * rand(size(A.PN)));
    // % B.PN = B.PN ;%.* 1.5; % (1 + 0.25 * rand(size(B.PN)));
    
    
    // generate the adjacency matrices for the connections between PN and
    // SICNS (PN's and sympathetic intrinsic cardiac NS)
    // A.SICNS_to_PN = couple_two_networks(N.ICNS, N.PN, p.S_to_P);
    double **A_SICNS_to_PN;
    A_SICNS_to_PN = new double *[N.ICNS];
    for( int id1 = 0; id1 < N.ICNS; id1++ ) {
        A_SICNS_to_PN[id1] = A.SICNS_to_PN[id1];
    }
    couple_two_networks( N.ICNS, N.PN, p.S_to_P, A_SICNS_to_PN, generator );
    
    // A.PN_to_SICNS = couple_two_networks(N.PN, N.ICNS, p.P_to_S);
    double **A_PN_to_SICNS;
    A_PN_to_SICNS = new double *[N.PN];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        A_PN_to_SICNS[id1] = A.PN_to_SICNS[id1];
    }
    couple_two_networks( N.PN, N.ICNS, p.P_to_S, A_PN_to_SICNS, generator );
    
    
    
    // Generate which sympathetic neurons are excitatory and inhibitory
    double prob = 0.5;
    int excit_S[ M ];
    double number;
    int excit_P[ N.PN ], SICNS_excit[ N.PN ];
    
    // if( strcmp( set_up_new_E_I_dist, "yes" ) == 0 ){
    // probability threshold, above is excitatory, below is inhibitory
    prob = 0.5;
    
    // generate a new list of excitatory and inhibitory neurons.
    // data = rand(M,1);
    
    for( int id1 = 0; id1 < M; id1 ++ ) {
        number = urand( generator );
        excit_S[id1] = 1;
        if( number < prob ) {
            excit_S[id1] = 0;
            
        }
    }
    
    
    // excit_S = find(data >= prob);
    // inhib_S = find(data < prob);
    
    // data = rand(N.PN,1);
    
    for( int id1 = 0; id1 < N.PN; id1 ++ ) {
        number = urand( generator );
        excit_P[id1] = 1;
        if( number < prob ) {
            excit_P[id1] = 0;
        }
    }
    
    // excit_P = find(data >= prob);
    // inhib_P = find(data < prob);
    // fprintf('New Excitatory and Inhibitory Markers Generated. \n')
    
   
    
    // Generate Synaptic Coupling Strength "distributions" for Sympathetic
    structdistribution avg, std_dev;
    
    // set means for distributions
    avg.EE = 0.5;
    avg.EI = 5.;
    avg.IE = 300.;
    avg.II = 0.5; // 10
    
    // set standard deviations for distributions
    std_dev.EE = 0.1;
    std_dev.EI = 0.1;
    std_dev.IE = 0.1;
    std_dev.II = 0.1;
    std::normal_distribution<double> normrndEE( avg.EE, std_dev.EE );
    std::normal_distribution<double> normrndEI( avg.EI, std_dev.EI );
    std::normal_distribution<double> normrndIE( avg.IE, std_dev.IE );
    std::normal_distribution<double> normrndII( avg.II, std_dev.II );
    structW W;
    
    
    
    for( int id1 = 0; id1 < M; id1++ ) {
        for( int id2 = 0; id2 < M; id2++ ) {
            W.EE[id1][id2] = 0;
            W.EI[id1][id2] = 0;
            W.IE[id1][id2] = 0;
            W.II[id1][id2] = 0;
            
            if( excit_S[id1] == 1 ) {
                if( excit_S[id2] == 1 ) {
                    W.EE[id1][id2] = A_symp.M[id1][id2] * fabs( normrndEE( generator ) );
                } else {
                    W.EI[id1][id2] = A_symp.M[id1][id2] * fabs( normrndEI( generator ) );
                }
            } else {
                if( excit_S[id2] == 1 ) {
                    W.IE[id1][id2] = A_symp.M[id1][id2] * fabs( normrndIE( generator ) );
                } else {
                    W.II[id1][id2] = A_symp.M[id1][id2] * fabs( normrndII( generator ) );
                    
                }
            }
            
            
        }
        
    }
    
    // SICNS_excit = excit_S(excit_S > (N.CNS + N.STELLATE)) - (N.CNS + N.STELLATE);
    // SICNS_inhib = inhib_S(inhib_S > (N.CNS + N.STELLATE)) - (N.CNS + N.STELLATE);
    for( int id1 = (N.CNS + N.STELLATE); id1 < M; id1++ ) {
        int id2 = id1 - ( N.CNS + N.STELLATE );
        SICNS_excit[ id2 ] = 0;
        if( excit_S[id1] == 1 ) {
            SICNS_excit[ id2 ] = 1;
        }
    }
    
    
    // A.PN_to_SICNS_EE(excit_P, SICNS_excit) = A.PN_to_SICNS(excit_P, SICNS_excit) .* abs(normrnd(avg.EE, std_dev.EE, [length(excit_P), length(SICNS_excit)]));
    // A.PN_to_SICNS_EI(excit_P, SICNS_inhib) = A.PN_to_SICNS(excit_P, SICNS_inhib) .* abs(normrnd(avg.EI, std_dev.EI, [length(excit_P), length(SICNS_inhib)]));
    // A.PN_to_SICNS_IE(inhib_P, SICNS_excit) = A.PN_to_SICNS(inhib_P, SICNS_excit) .* abs(normrnd(avg.IE, std_dev.IE, [length(inhib_P), length(SICNS_excit)]));
    // A.PN_to_SICNS_II(inhib_P, SICNS_inhib) = A.PN_to_SICNS(inhib_P, SICNS_inhib) .* abs(normrnd(avg.II, std_dev.II, [length(inhib_P), length(SICNS_inhib)]));
    
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            A.PN_to_SICNS_EE[ id1][ id2 ] = 0;
            A.PN_to_SICNS_EI[ id1][ id2 ] = 0;
            A.PN_to_SICNS_IE[ id1][ id2 ] = 0;
            A.PN_to_SICNS_II[ id1][ id2 ] = 0;
            
            if( excit_P[id1] == 1 ) {
                if( SICNS_excit[id2] == 1 ) {
                    A.PN_to_SICNS_EE[ id1][ id2 ] = A.PN_to_SICNS[ id1][ id2 ] * fabs( normrndEE( generator ) );
                } else {
                    A.PN_to_SICNS_EI[ id1][ id2 ] = A.PN_to_SICNS[ id1][ id2 ] * fabs( normrndEI( generator ) );
                }
            } else {
                if( SICNS_excit[id2] == 1 ) {
                    A.PN_to_SICNS_IE[ id1][ id2 ] = A.PN_to_SICNS[ id1][ id2 ] * fabs( normrndIE( generator ) );
                } else {
                    A.PN_to_SICNS_II[ id1][ id2 ] = A.PN_to_SICNS[ id1][ id2 ] * fabs( normrndII( generator ) );
                }
            }
            
        }
        
    }
   
    
    // A.SICNS_to_PN_EE(SICNS_excit, excit_P) = A.SICNS_to_PN(SICNS_excit,excit_P) .* abs(normrnd(avg.EE, std_dev.EE, [length(SICNS_excit), length(excit_P)]));
    // A.SICNS_to_PN_EI(SICNS_excit, inhib_P) = A.SICNS_to_PN(SICNS_excit, inhib_P) .* abs(normrnd(avg.EI, std_dev.EI, [length(SICNS_excit), length(inhib_P)]));
    // A.SICNS_to_PN_IE(SICNS_inhib, excit_P ) = A.SICNS_to_PN(SICNS_inhib, excit_P) .* abs(normrnd(avg.IE, std_dev.IE, [length(SICNS_inhib), length(excit_P)]));
    // A.SICNS_to_PN_II(SICNS_inhib, inhib_P ) = A.SICNS_to_PN(SICNS_inhib, inhib_P) .* abs(normrnd(avg.II, std_dev.II, [length(SICNS_inhib), length(inhib_P)]));
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            A.SICNS_to_PN_EE[ id1][ id2 ] = 0;
            A.SICNS_to_PN_EI[ id1][ id2 ] = 0;
            A.SICNS_to_PN_IE[ id1][ id2 ] = 0;
            A.SICNS_to_PN_II[ id1][ id2 ] = 0;
            if( SICNS_excit[id1] == 1 ) {
                if( excit_P[id2] == 1 ) {
                    A.SICNS_to_PN_EE[ id1][ id2 ] = A.SICNS_to_PN[ id1][ id2 ] * fabs( normrndEE( generator ) );
                } else {
                    A.SICNS_to_PN_EI[ id1][ id2 ] = A.SICNS_to_PN[ id1][ id2 ] * fabs( normrndEI( generator ) );
                }
            } else {
                if( excit_P[id2] == 1 ) {
                    A.SICNS_to_PN_IE[ id1][ id2 ] = A.SICNS_to_PN[ id1][ id2 ] * fabs( normrndIE( generator ) );
                } else {
                    A.SICNS_to_PN_II[ id1][ id2 ] = A.SICNS_to_PN[ id1][ id2 ] * fabs( normrndII( generator ) );
                }
            }
        }
        
    }
    
    // A.PN_EE(excit_P, excit_P) = A.PN(excit_P,excit_P) .* abs(normrnd(avg.EE, std_dev.EE, [length(excit_P), length(excit_P)]));
    // A.PN_EI(excit_P, inhib_P) = A.PN(excit_P,inhib_P) .* abs(normrnd(avg.EI, std_dev.EI, [length(excit_P), length(inhib_P)]));
    // A.PN_IE(inhib_P, excit_P) = A.PN(inhib_P,excit_P) .* abs(normrnd(avg.IE, std_dev.IE, [length(inhib_P), length(excit_P)]));
    // A.PN_II(inhib_P, inhib_P) = A.PN(inhib_P,inhib_P) .* abs(normrnd(avg.II, std_dev.II, [length(inhib_P), length(inhib_P)]));
    
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        for( int id2 = 0; id2 < N.PN; id2++ ) {
            A.PN_EE[ id1][ id2 ] = 0;
            A.PN_EI[ id1][ id2 ] = 0;
            A.PN_IE[ id1][ id2 ] = 0;
            A.PN_II[ id1][ id2 ] = 0;
            if( excit_P[id1] == 1 ) {
                if( excit_P[id2] == 1 ) {
                    A.PN_EE[ id1][ id2 ] = A.PN[ id1][ id2 ] * fabs( normrndEE( generator ) );
                } else {
                    A.PN_EI[ id1][ id2 ] = A.PN[ id1][ id2 ] * fabs( normrndEI( generator ) );
                }
            } else {
                if( excit_P[id2] == 1 ) {
                    A.PN_IE[ id1][ id2 ] = A.PN[ id1][ id2 ] * fabs( normrndIE( generator ) );
                } else {
                    A.PN_II[ id1][ id2 ] = A.PN[ id1][ id2 ] * fabs( normrndII( generator ) );
                }
            }
        }
        
    }
    
    // fprintf('New Synaptic Coupling Weights Generated \n')
    
    // Set parameters
    
    // Reset potential
    p.v_reset = -68.;
    
    // Threshold potential
    p.v_T = -52.;
    
    // Reversal potentials
    p.EL = -70.;
    p.EK = -90.;
    p.ESyn_E = 0.; // Excitatory
    p.ESyn_I = -80.; // Inhibitory
    p.E_out = p.ESyn_E; // treating external input as excitatory.
    
    // Steady state of n
    p.nss = 0.;
    
    // Time constant of n
    p.tau_n = 75.;
    p.tau_w = 165.;
    p.tau_w_reset = 1250.;
    double tau_SynE = 10.; // synaptic time-constant for vagal input
    double tau_ACh = 1234.;
    
    double tau_ISO = 1234.;
    
    // Conductances
    p.gL = 0.01; // leakage
    p.gK = 100.; // potassium
    p.g_out = 0.00175; // external input conductance limiter
    
    // capacitance
    p.C = 170. / 220.;
    // set max conductance of M-current
    double gx, gy;
    if( strcmp(set_up_new_M_current_dist, "yes") == 0 ) {
        
        for( int id1 = 0; id1 < M; id1++ ) {
            
            gx = gammarand2( generatorg2 );
            gy = gammarand1( generatorg1 );
            p.gM_S[id1] = 10. * ( gx / ( gx + gy ) )  + 0.;
            
        }
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            gx = gammarand2( generatorg2 );
            gy = gammarand1( generatorg1 );
            p.gM_PN[id1] = 10. * ( gx / ( gx + gy ) )  + 0.;
            
        }
        
    }
    // end
    
    p.gSyn_E = 0.0005; // Excitatory synaptic
    p.gSyn_I = 0.0005; // Inhibitory synaptic
    
    // Setting up parameters for system
    // Choose alpha > beta.
    double alpha = 0.5;
    double beta = 0.25;
    
    
    // v_S = zeros(M,  1);
    double v_S[M], v1_S[M], v_S_old[M], k1_v[M], k2_v[M];
    // n_S =  zeros(M, 1);
    double n_S[M], n1_S[M], n_S_old[M];
    // x_S = zeros(2 * M,  1); // x = [ s1 ; s2 ; u1 ; u2]
    double x_S[ 2*M ], x1_S[ 2*M ], x_S_old[ 2*M ];
    // w_S = zeros(M, 1);
    double w_S[ M ], w1_S[ M ], w_S_old[ M ], k1_w[ M ], k2_w[ M ];
    double n_spike[M], vector[M], x_spike[2*M], w_spike[M];
    
    // v_PN = zeros(N.PN, 1);
    double v_PN[ N.PN ], v1_PN[ N.PN ], v_PN_old[ N.PN ], v_PN_old_k1[ N.PN ];
    // n_PN =  zeros(N.PN,  1);
    double n_PN[ N.PN ], n1_PN[ N.PN ], n_PN_old[ N.PN ];
    // x_PN= zeros(2 * N.PN,  1); % x = [ s1 ; s2 ; u1 ; u2]
    double x_PN[ 2*N.PN ], x1_PN[ 2*N.PN ], x_PN_old[ 2*N.PN ] ;
    // w_PN = zeros(N.PN, 1);
    double w_PN[ N.PN ], w1_PN[ N.PN ], w_PN_old[ N.PN ], w_PN_old_k1[ N.PN ];
    
    // gVagal = zeros(N.PN,1); % Vagal conductance
    double gVagal[ N.PN ], gVagal1[ N.PN ];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        gVagal[ id1 ] = 0;
        gVagal1[ id1 ] = 0;
    }
    
    // timer_S = zeros(M,1);
    double timer_S[ M ];
    for( int id1 = 0; id1 < M; id1++ ) {
        timer_S[ id1 ] = 0;
    }
    // timer_PN = zeros(N.PN,1);
    double timer_PN[ N.PN ];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        timer_PN[ id1 ] = 0;
    }
    
    double T_ref = 0.5;
    
    //% Impose initial conditions
    // v_S(:) = p.EL; //%* ones(N,1) + randi([-5,5],N,1);
    for( int id1 = 0; id1 < M; id1++ ) {
        v_S[ id1 ] = p.EL;
        n_S[ id1 ] = 0;
        w_S[ id1 ] = 0;
    }
    //% vectorize s and u to take advantage of liner algebra
    //% x = [s1, s2, ... , sN, u1, u2, ... , uN]
    // x_S(1:M) = 0; % s
    // x_S(M+1:end) = 0; % u
    for( int id1 = 0; id1 < M; id1++ ) {
        x_S[ id1 ] = 0;
        x1_S[ id1 ] = 0;
    }
    for( int id1 = M; id1 < (2*M); id1++ ) {
        x_S[ id1 ] = 0;
        x1_S[ id1 ] = 0;
    }
    
    // v_PN(:) = p.EL; %* ones(N,1) + randi([-5,5],N,1);
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        v_PN[ id1 ] = p.EL;
        n_PN[ id1 ] = 0;
        n1_PN[ id1 ] = 0;
        w_PN[ id1 ] = 0;
    }
    //% vectorize s and u to take advantage of liner algebra
    //% x = [s1, s2, ... , sN, u1, u2, ... , uN]
    // x_PN(1:N.PN) = 0; % s
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        x_PN[ id1 ] = 0;
        x1_PN[ id1 ] = 0;
    }
    // x_PN(N.PN+1:end) = 0; % u
    for( int id1 = N.PN; id1 < (2*N.PN); id1++ ) {
        x_PN[ id1 ] = 0;
        x1_PN[ id1 ] = 0;
    }
    
    // % initalize vectors to log firing times of neurons in PN and sympathetic.
    // fire_count_S = zeros(M,T/dt + 1);
    double fire_count_S[ M ];
    for( int id1 = 0; id1 < M; id1++ ) {
        fire_count_S[ id1 ] = 0;
    }
    //fire_count_PN = zeros(N.PN, T/dt + 1);
    double fire_count_PN[ N.PN ];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        fire_count_PN[ id1 ] = 0;
    }
    //fire_count_CNS = zeros(N.CNS, T/dt + 1);
    double fire_count_CNS[ N.CNS ];
    for( int id1 = 0; id1 < N.CNS; id1++ ) {
        fire_count_CNS[ id1 ] = 0;
    }
    //fire_count_VAGUS = zeros(N.PN, T/dt + 1);
    double fire_count_VAGUS[ N.PN ];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        fire_count_VAGUS[ id1 ] = 0;
    }
    
    // ISO = zeros(M, T/dt + 1);
    double ISO[ M ];
    for( int id1 = 0; id1 < M; id1++ ) {
        ISO[ id1 ] = 0;
    }
    
    param.dISO = 0.5;
    
    //% intialize amount of ACh released
    // ACh = zeros(N.PN, T/dt + 1);
    double ACh[ N.PN ];
    for( int id1 = 0; id1 < N.PN; id1++ ) {
        ACh[ id1 ] = 0;
    }
    // K_ACh_activation = zeros(1, T/dt + 1);
    double K_ACh_activation = 0;
    param.dACh = 0.5;
    
    //%% Generate external input  for SNS - this can be modified to whatever is necessary
    //% uses a poisson process to draw random firing times which generates an
    //% EPSP in the sympathetic portion of the model.
    
    param.ds =  20. * sqrt(dt);
    param.lambda1 =  0.5;
    param.lambda2 =  0.25;
    
    double drop = 0.6;
    
    structtone tone;
    tone.SICNS = 0.1;
    tone.STELLATE = 0.5;
    tone.PICNS = 0;
    
    structP_ICNS P_ICNS;
    P_ICNS.max = 0.8;
    P_ICNS.min = 0.2;
    
    structP_STELLATE P_STELLATE;
    P_STELLATE.max = 0.8;
    P_STELLATE.min = 0.4;
    
    // intensity_function = @(t) intensity_function_vector(t, tone, P_ICNS, P_STELLATE, drop, N, N.PN);
    
    
    //%%%%%%%% Set up vagal input to PN network - probability threshold function
    //% I_Vagal is the probability threshold function (pulse train with non-zero
    //% tone) This is used similar to the symapthetic portion, but the rise time
    //% of the synaptic conductance is instantaneous.
    
    double Vagal_Tone = 8.; // % basal firing threshold for vagus
    
    
    
    double Vaga_tone_start = 0; //20000. / dt;
    double Vaga_tone_end = T / dt; //80000. / dt;
    
    pulse_length = ( Vaga_tone_end - Vaga_tone_start ) + 1.;
    
    //% % Set up I_Vagal
    // I_Vagal = 1/10 * ones(N.PN, T/dt + 1); % initialize I_Vagal
    double I_Vagal; //[ N.PN ];
    
    
    //%% Set up stellate stimulation
    
   
    
    //% Set up I_stellate
    //I_stellate = zeros(M,1);
    double I_stellate[ M ];
    for( int id1 = 0; id1 < M; id1++ ) {
        I_stellate[id1] = 0;
    }
    
    //%%
    //x_out = zeros(2*(M + N.PN),1);
    double x_out[ 2 * ( M + N.PN ) ], x_out_old[ 2 * ( M + N.PN ) ];
    for( int id1 = 0; id1 < 2*(M+N.PN); id1++ ) {
        x_out[id1] = 0;
        x_out_old[id1] = 0;
    }
    
    
    //%% Run simulation
    
    // double t;
    int m;
    double intensities[ M + N.PN ];
    double fire[ M + N.PN ], fire_CNS[ M + N.PN ], t_s[ M + N.PN ], dt1[ M + N.PN ], dt2[ M + N.PN ], dn[ M + N.PN ];
    double delta[ N.PN ], k1_E[ N.PN ], k2_E[ N.PN ];
    double ref_PN[ N.PN ], ref_S[ M ];
    
    double v_reset[ M ];
    for( int id1 = 0; id1 < M; id1++ ) {
        v_reset[ id1 ] = p.v_reset;
    }
    double vk1v[ M ], wk1w[ M ];
    double rand_num[ M + N.PN ];
    double v_S_old_k1vdt[M], w_S_old_k1wdt[M];
    double a, b, AChxPN;
    double sumACh = 0;
    int fire_size = 0;
    double sumISO = 0;
    
    //for m = 1 : T/dt
    for( t = dt; t <= T; t+=dt ) {
        
        m = (int) ( t / dt );
        
        I_Vagal = 0.1;
        if( m >= Vaga_tone_start && m < Vaga_tone_end ) {
            I_Vagal = Vagal_Tone;
        }
        //ACh(:, m + 1) = ACh(:, m) + dt * (-ACh(:, m)/tau_ACh);
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            ACh[id1] = ACh[id1] + dt * ( -ACh[id1] / tau_ACh );
        }
        for( int id1 = 0; id1 < M; id1++ ) {
            //ISO(:, m + 1) = ISO(:,m) + dt * (-ISO(:,m)/tau_ISO);
            ISO[id1] = ISO[id1] + dt * ( -ISO[id1] / tau_ISO );
        }
        
        for( int id1 = 0; id1 < 2 * ( M + N.PN ) ; id1++ ){
            x_out_old[id1] = x_out[id1];
        }
        //%%%%%%%% Cardiac feedback and Central drive %%%%%%%%%%%%%
        //%     % generate random number for pre-synaptic each neuron
        // rand_num = rand(M + N.PN,1);
        for( int id1 = 0; id1 < M + N.PN; id1++ ) {
            rand_num[ id1 ] = urand( generator );
        }
        
        intensity_function_vector( t, tone, P_ICNS, P_STELLATE, drop, N, N.PN, intensities );
        //fire_size = 0;
        for( int id1 = 0; id1 < (M+N.PN); id1++ ) {
            fire[ id1 ] = 0;
            if( rand_num[ id1 ] <= (intensities[ id1 ] * 0.005) ) {
                fire[ id1 ] = 1;
                
            }
        }
        
        for( int id1 = 0; id1 < N.CNS; id1++ ) {
            fire_CNS[ id1 ] = fire[ id1 ];
        }
        
        //fire_count_CNS(fire_CNS,m) = fire_count_CNS(fire_CNS,m) + 1;
        for( int id1 = 0; id1 < N.CNS; id1++ ) {
            if( fire_CNS[ id1 ] > 0) {
                fire_count_CNS[ id1 ] += 1;
            }
        }
        
        //x_out(1: (M + N.PN)) = s_update_stochastic_input(dt, x_out_old(1:M + N.PN), x_out(M + N.PN + 1 : end), param.lambda1, param.lambda2);
        for( int id1 = 0; id1 < (M+N.PN); id1++ ) {
            s_update_stochastic_input( dt, x_out_old[id1], x_out[ id1 + M + N.PN ], param.lambda1, param.lambda2, x_out[id1] );
        }
        //x_out( M + N.PN + 1 : end) = u_update_stochastic_input(dt, x_out_old(1: M + N.PN), x_out( M + N.PN + 1 : end) , param.lambda1, param.lambda2);
        for( int id1 = 0; id1 < (M+N.PN); id1++ ){
            u_update_stochastic_input( dt, x_out_old[ id1 ], x_out[ id1 + M + N.PN ], param.lambda1, param.lambda2, x_out[ id1 + M + N.PN ]);
        }
        //x_out(fire + M + N.PN) = u_update_stochastic_input(dt, x_out_old(fire), x_out(fire + M + N.PN) , param.lambda1, param.lambda2) + param.ds;
        for( int id1 = 0; id1 < (M+N.PN); id1++ ) {
            if( fire[id1] == 1 ){
                u_update_stochastic_input( dt, x_out_old[id1], x_out[ id1 + M + N.PN] , param.lambda1, param.lambda2, x_out[ id1 + M + N.PN ] ) ;
                x_out[ id1 + M + N.PN ] += param.ds;
            }
        }
        //%
        //%%%%%%%%%%%%%%%%%% Vagal input %%%%%%%%%%%%%%%%%%%%%%%%%%
        //% draw from probability threshold to generate vagal input to PN network
        // rand_num = rand(N.PN,1);
        for( int id1 = 0; id1 < N.PN; id1 ++ ) {
            rand_num[ id1 ] = urand( generator );
        }
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            delta[ id1 ] = 0.;
            if( rand_num[ id1 ] <= I_Vagal / 200 ) {
                delta[ id1 ] = 0.75 * 10. * sqrt( dt );
                fire_count_VAGUS[ id1 ] += 1;
            }
        }
        
        //% update the synaptic conductance to PN cell from vagus
       
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            k1_E[ id1 ] = ( -gVagal[ id1 ] + delta[id1 ] ) / tau_SynE;
            k2_E[ id1 ] = ( -gVagal[ id1 ] + k1_E[ id1 ] * dt + delta[ id1 ] ) / tau_SynE;
            gVagal[ id1 ] = gVagal[ id1 ] + dt * (k1_E[ id1 ] + k2_E[ id1 ]) * 0.5;
        }
        
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%% Update PN cells %%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            n_PN_old[ id1 ] = n_PN[ id1 ];
            v_PN_old[ id1 ] = v_PN[ id1 ];
            w_PN_old[ id1 ] = w_PN[ id1 ];
            x_PN_old[ id1 ] = x_PN[ id1 ];
            x_PN_old[ N.PN + id1 ] = x_PN[ N.PN + id1 ];
            ref_PN[ id1 ] = 0;
            if( timer_PN[ id1 ] == 0 ) {
                ref_PN[ id1 ] = 1.;
            };
        }
        
        
        for( int id1 = 0; id1 < M; id1++ ) {
            n_S_old[ id1 ] = n_S[ id1 ];
            v_S_old[ id1 ] = v_S[ id1 ];
            w_S_old[ id1 ] = w_S[ id1 ];
            x_S_old[ id1 ] = x_S[ id1 ];
            x_S_old[ id1 + M ] = x_S[ id1 + M ];
            ref_S[ id1 ] = 0;
            if( timer_S[id1] == 0 ) {
                ref_S[ id1 ] = 1.;
            }
        }
        //
        //        % take standard RK2 step for v, M
        //            k1_v = driving_function_PN_efferent(v_PN_old, n_PN_old, p, gVagal, x_PN_old(1:N.PN), w_PN_old, x_S_old(1:M), x_out_old(M+1:M + N.PN), ref_PN, A, B, N);
        
        driving_function_PN_efferent( v_PN_old, n_PN_old, p, gVagal, x_PN_old, w_PN_old, x_S_old, &x_out_old[M], ref_PN, A, B, N, k1_v );
        
        //        k1_w = 1/p.tau_w * (w_ss(v_PN_old) - w_PN_old);
        //
        //
        //        % find exact solutions at next time step for n, s, u
        //            n_PN = exp(-dt/p.tau_n) * n_PN_old;
        //        x_PN(1:N.PN) = s_update_stochastic_input(dt, x_PN_old(1:N.PN), x_PN_old(N.PN+1:end), alpha, beta);
        //        x_PN(N.PN+1:end) = u_update_stochastic_input(dt, x_PN_old(1:N.PN), x_PN_old(N.PN+1:end), alpha, beta);
        for( int id1 = 0; id1 < N.PN; id1 ++ ) {
            
            k1_w[ id1 ] = 1. / p.tau_w * ( w_ss( v_PN_old[ id1 ] ) - w_PN_old[ id1 ] );
            
            n_PN[ id1 ] = exp( -dt / p.tau_n ) * n_PN_old[ id1 ];
            s_update_stochastic_input(dt, x_PN_old[id1], x_PN_old[N.PN+id1], alpha, beta, x_PN[id1] );
            u_update_stochastic_input(dt, x_PN_old[id1], x_PN_old[N.PN+id1], alpha, beta, x_PN[N.PN+id1] );
            
            v_PN_old_k1[id1] = v_PN_old[ id1 ] + k1_v[ id1 ] * dt;
            w_PN_old_k1[id1] = w_PN_old[ id1 ] + k1_w[ id1 ] * dt;
            
        }
        
        
        
        driving_function_PN_efferent( v_PN_old_k1, n_PN, p, gVagal, x_PN, w_PN_old_k1, x_S, &x_out[M], ref_PN, A, B, N, k2_v );
        
        for( int id1 = 0; id1 < N.PN; id1 ++ ) {
            k2_w[ id1 ] = 1. / p.tau_w * ( w_ss( v_PN_old_k1[ id1 ] ) - w_PN_old_k1[ id1 ] );
            
            v_PN[ id1 ] = v_PN_old[ id1 ] + dt * ( k1_v[ id1 ] + k2_v[ id1 ] ) * 0.5;
            w_PN[ id1 ] = w_PN_old[ id1 ] + dt * ( k1_w[ id1 ] + k2_w[ id1 ] ) * 0.5;
        }
        
        
        //        % find neurons that fired
        //        fire = find(v_PN >= p.v_T);
        //        % log the firing times of each neuron
        //        fire_count_PN(fire, m+1) = fire_count_PN(fire, m+1) + 1;
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            fire[ id1 ] = 0;
            if( v_PN[ id1 ] >= p.v_T ) {
                fire[ id1 ] = 1;
                fire_count_PN[ id1 ] += 1;
                fire_size += 1;
            }
            
        }
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            
            //        % if any neurons fired, update using linear interpolant modified RK2
            //            if (size(fire,1) * size(fire,2) ~= 0)
            //
            //                % set refractory period for neurons that fire
            //                    timer_PN(fire) = T_ref/dt;
            //
            
            n_spike[ id1 ] = n_PN[ id1 ];
            vector[ id1 ] = 0;
            x_spike[ id1 ] = x_PN[ id1 ];
            x_spike[ N.PN + id1 ] = x_PN[ N.PN + id1 ];
            w_spike[ id1 ] = 0;
            
            if( fire[ id1 ] == 1 ){
                timer_PN[ id1 ] = T_ref / dt;
                
                
                //        % spike times for each neuron that fired.
                //            t_s = (p.v_T - (1+m) * v_PN_old(fire) + m * v_PN(fire))*dt ./ (v_PN(fire) - v_PN_old(fire));
                t_s[ id1 ] = ( p.v_T - ( 1 + m ) * v_PN_old[id1] + m * v_PN[id1] ) * dt / ( v_PN[id1] - v_PN_old[id1] );
                
                //        % separate grid points into two intervals - this is for linear
                //            % interpolation to reduce numerical error
                //            dt1 = t_s - m * dt; dt2 = (m+1) * dt - t_s;
                dt1[ id1 ] = t_s[ id1 ] - m * dt;
                dt2[ id1 ] = ( m + 1 ) * dt - t_s[ id1 ];
                
                
                //        % spike-time values of n,s,u for pre-synaptic neurons
                //            n_spike = n_PN;
                //        n_spike(fire) = exp(-dt1/p.tau_n) .* n_PN_old(fire);
                n_spike[ id1 ] = exp( -dt1[ id1 ] / p.tau_n ) * n_PN_old[ id1 ];
                //        dn = (n_spike(fire) - 1) * (exp(-T_ref/p.tau_n)-1);
                dn[ id1 ] = ( n_spike[ id1 ] - 1. ) * ( exp( -T_ref / p.tau_n ) - 1. );
                //        n_spike(fire) = n_spike(fire) + dn;
                n_spike[ id1 ] += dn[ id1 ];
                //
                //        x_spike = x_PN;
                //        x_spike(fire) = s_update_stochastic_input(dt1, x_PN_old(fire), x_PN_old(fire + N.PN), alpha, beta);
                s_update_stochastic_input( dt1[id1], x_PN_old[id1], x_PN_old[id1+N.PN], alpha, beta, x_spike[id1] );
                //        x_spike(fire + N.PN) = u_update_stochastic_input(dt1, x_PN_old(fire), x_PN_old(fire + N.PN), alpha, beta) + 1;
                u_update_stochastic_input( dt1[id1], x_PN_old[id1], x_PN_old[id1 + N.PN], alpha, beta, x_spike[id1+N.PN]);
                x_spike[ id1+N.PN] += 1.;
                //        % spike-time values for m-current activation
                //            vector = zeros(N.PN,1);
                //        vector(fire) = dt1;
                vector[ id1 ] = dt1[ id1 ];
                
                k1_w[ id1 ] = 1. / p.tau_w * ( w_ss(v_PN_old[id1]) - w_PN_old[id1] );
                
                //        k2_w = 1/p.tau_w * (w_ss(v_PN_old + k1_v .* vector) - (w_PN_old + k1_w .* vector));
                k2_w[ id1 ] = 1. / p.tau_w * ( w_ss(v_PN_old[id1] + k1_v[id1] * vector[id1]) - (w_PN_old[id1] + k1_w[id1] * vector[id1] ) );
                
                //        w_spike = w_PN_old + vector/2 .* (k1_w + k2_w) + (w_PN_old - 1) * (exp(-1/p.tau_w_reset) - 1);
                w_spike[ id1 ] = w_PN_old[ id1 ] + vector[ id1] * 0.5 * (k1_w[id1] + k2_w[id1] ) + (w_PN_old[id1] - 1. ) * ( exp( -1. / p.tau_w_reset ) - 1. );
                
                
                //        % calculate new exact value of n after reset dynamics
                //        n_PN(fire) = n_spike(fire) .* exp(-dt2/p.tau_n);
                n_PN[ id1 ] = n_spike[ id1 ] * exp( -dt2[ id1 ] / p.tau_n );
            }
        }
        //        % calulate updated values for v and w after reset dynamics
        //            k1_v = driving_function_PN_efferent(p.v_reset * ones(N.PN,1), n_spike, p, gVagal, x_spike(1:N.PN), w_spike, x_S(1:M), x_out_old(M+1:M+N.PN), ref_PN, A, B, N);
        
        driving_function_PN_efferent( v_reset, n_spike, p, gVagal, x_spike, w_spike, x_S, &x_out_old[M], ref_PN, A, B, N, k1_v );
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            //        k1_w = 1/p.tau_w * (w_ss(p.v_reset*ones(N.PN,1)) - w_spike);
            k1_w[ id1 ] = 1. / p.tau_w * ( w_ss(p.v_reset) - w_spike[ id1 ] );
            //
            vector[ id1 ] = 0;
            if( fire[ id1 ] == 1 ) {
                //        vector = zeros(N.PN,1);
                //        vector(fire) = dt2;
                vector[ id1 ] = dt2[ id1 ];
            }
        }
        for( int id1 = 0; id1 < M; id1++ ) {
            vk1v[ id1 ] = p.v_reset + k1_v[ id1 ] * vector[ id1 ];
        }
        //        k2_v = driving_function_PN_efferent(p.v_reset * ones(N.PN,1) + k1_v .* vector, n_spike, p, gVagal, x_spike(1:N.PN), w_spike, x_S(1:M), x_out(M+1:M+N.PN), ref_PN, A, B, N);
        driving_function_PN_efferent( vk1v, n_spike, p, gVagal, x_spike, w_spike, x_S, &x_out[M], ref_PN, A, B, N, k2_v );
        
        for( int id1 = 0; id1<N.PN; id1++ ) {
            //        k2_w = 1/p.tau_w * (w_ss(p.v_reset*ones(N.PN,1) + k1_v .* vector) - (w_spike + k1_w .* vector));
            k2_w[ id1 ] = 1. / p.tau_w * ( w_ss(vk1v[id1] ) - ( w_spike[id1] + k1_w[id1] * vector[id1] ) );
            //
            if( fire[ id1 ] == 1 ) {
                //        v_PN(fire) = p.v_reset + dt2/2 .* (k1_v(fire) + k2_v(fire));
                v_PN[ id1 ] = p.v_reset + dt2[ id1 ] * 0.5 * ( k1_v[ id1 ] + k2_v[ id1 ] );
                //        w_PN(fire) = w_spike(fire) + dt2/2 .* (k1_w(fire) + k2_w(fire));
                w_PN[ id1 ] = w_spike[ id1 ] + dt2[ id1 ] * 0.5 * ( k1_w[ id1 ] + k2_w[ id1 ] );
            }
        }
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            
            if( fire[ id1 ] == 1 ) {
                //
                //        % update s and u after reset dynamics
                //        x_PN(fire) = s_update_stochastic_input(dt2, x_spike(fire), x_spike(fire + N.PN), alpha, beta);
                s_update_stochastic_input( dt2[id1], x_spike[id1], x_spike[id1 + N.PN], alpha, beta, x_PN[ id1 ] );
                //        x_PN(fire + N.PN) = u_update_stochastic_inputdt2, x_spike(fire), x_spike(fire + N.PN), alpha, beta);
                u_update_stochastic_input( dt2[id1], x_spike[id1], x_spike[id1 + N.PN], alpha, beta, x_PN[id1+N.PN] );
                //
                //        ACh(fire, m + 1) = ACh(fire, m) + param.dACh;
                ACh[ id1 ] += param.dACh;
                //        end
            }
        }
        //
        //        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //        %%% Update sympathetic cells %%%%%
        //        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //
        //        % take standard RK2 step for v, M
        //            k1_v = driving_function_stochastic_input(v_S_old, n_S_old, x_S_old(1:M), w_S_old, p, W, x_out_old(1:M), x_PN_old, I_stellate, ref_S,A, A_CNS_to_STELLATE, N);
        driving_function_stochastic_input( v_S_old, n_S_old, x_S_old, w_S_old, p, W, x_out_old, x_PN_old, I_stellate, ref_S, A, A_CNS_to_STELLATE, N, k1_v );
        

        for( int id1 = 0; id1 < M; id1++ ) {
            
            
            k1_w[ id1 ] = 1. / p.tau_w * ( w_ss(v_S_old[id1]) - w_S_old[id1] );
            
            //        % find new exact solutions at next time step for n, s, u
            //            n_S = exp(-dt/p.tau_n) * n_S_old;
            n_S[ id1 ] = exp( -dt / p.tau_n ) * n_S_old[ id1 ];
            
            //        x_S(1:M) = s_update_stochastic_input(dt, x_S_old(1:M), x_S_old(M+1:end), alpha, beta);
            s_update_stochastic_input( dt, x_S_old[id1], x_S_old[ M+id1 ], alpha, beta, x_S[id1] );
            
            //        x_S(M+1:end) = u_update_stochastic_input(dt, x_S_old(1:M), x_S_old(M+1:end), alpha, beta);
            u_update_stochastic_input( dt, x_S_old[id1], x_S_old[M+id1], alpha, beta, x_S[M+id1] );
            
            v_S_old_k1vdt[ id1 ] = v_S_old[ id1 ] + k1_v[ id1 ] * dt;
            w_S_old_k1wdt[ id1 ] = w_S_old[ id1 ] + k1_w[ id1 ] * dt;
        }
        //cout << endl;
        //
        //
        //        %driving_function_stochastic_input(v, n, s, w, p, W, s_out, x_PN, I_inj, A, N)
        //        k2_v = driving_function_stochastic_input(v_S_old + k1_v * dt, n_S, w_S_old + k1_w * dt, x_S(1:M), p, W, x_out(1:M), x_PN, I_stellate, ref_S, A, A_CNS_to_STELLATE, N);
        
        driving_function_stochastic_input( v_S_old_k1vdt, n_S, w_S_old_k1wdt, x_S, p, W, x_out, x_PN, I_stellate, ref_S, A, A_CNS_to_STELLATE, N, k2_v );
        
        for( int id1 = 0; id1 < M; id1++ ) {
            //        k2_w = 1/p.tau_w * (w_ss(v_S_old + k1_v * dt) - (w_S_old + k1_w * dt));
            k2_w[ id1 ] = 1. / p.tau_w * ( w_ss(v_S_old_k1vdt[ id1 ]) - w_S_old_k1wdt[ id1 ] );
            
            //        v_S = v_S_old + dt/2 * (k1_v + k2_v);
            v_S[ id1 ] = v_S_old[ id1 ] + dt * 0.5 * ( k1_v[ id1 ] + k2_v[ id1 ] );
            
            //        w_S = w_S_old + dt/2 * (k1_w + k2_w);
            w_S[ id1 ] = w_S_old[ id1 ] + dt * 0.5 * ( k1_w[ id1 ] + k2_w[ id1 ] );
        }
        for( int id1 = 0; id1 < M; id1++ ) {
            //
            //        %find neurons that fired (this can be optimized)
            //        fire = find(v_S >= p.v_T);
            fire[ id1 ] = 0;
            t_s[ id1 ] = t + dt;
            
            if( v_S[ id1 ] >= p.v_T ) {
                fire[ id1 ] = 1;
                
                
                //        % log the firing times of each neuron
                //        fire_count_S(fire, m+1) = fire_count_S(fire, m+1) + 1;
                fire_count_S[ id1 ] += 1;
               
                //        %      for j = 1 : length(fire)
                //            %         if isempty(fire) == 0
                //                %         fire_times_S{fire(j)} = [fire_times_S{fire(j)}, m*dt];
                //        %         end
                //        %     end
                //        % if any neurons fired
                //            if (size(fire,1) * size(fire,2) ~= 0)
                //                timer_S(fire) = T_ref/dt;
                timer_S[ id1 ] = T_ref / dt;
                //        %ref_S(fire, m:(m+T_ref/dt)) = 0;
                //
                //        % spike times for each neuron that fired.
                //            t_s = (p.v_T - (1+m) * v_S_old(fire) + m * v_S(fire))*dt ./ (v_S(fire) - v_S_old(fire));
                t_s[ id1 ] = ( p.v_T - ( 1 + m ) * v_S_old[ id1 ] + m * v_S[ id1 ] ) * dt / ( v_S[ id1 ] - v_S_old[ id1 ] );
                //cout << t << "\t" << t_s[id1] << endl;
            }
        }
        
        for( int id1 = 0; id1 < M; id1 ++ ) {
            //        % separate grid points into two intervals
            //        dt1 = t_s - m * dt; dt2 = (m+1) * dt - t_s;
            dt1[ id1 ] = t_s[ id1 ] - m * dt;
            dt2[ id1 ] = ( m + 1 ) * dt - t_s[ id1 ];
            
            //
            //        % spike-time values of n,s,u for pre-synaptic neurons
            //            n_spike = n_S;
            n_spike[ id1 ] = n_S[ id1 ];
            if( fire[ id1 ] == 1 ) {
                //        n_spike(fire) = exp(-dt1/p.tau_n) .* n_S_old(fire);
                n_spike[ id1 ] = exp( -dt1[ id1 ] / p.tau_n ) * n_S_old[ id1 ];
                //        dn = (n_spike(fire) - 1) * (exp(-T_ref/p.tau_n)-1);
                dn[ id1 ] = ( n_spike[ id1 ] - 1. ) * ( exp( -T_ref / p.tau_n ) - 1. );
                //        n_spike(fire) = n_spike(fire) + dn;
                n_spike[ id1 ] = n_spike[ id1 ] + dn[ id1 ];
            }
        }
        for( int id1 = 0; id1 < 2*M; id1++ ) {
            //        x_spike = x_S;
            x_spike[ id1 ] = x_S[ id1 ];
        }
        for( int id1 = 0; id1 < M; id1++ ) {
            if( fire[id1] == 1 ) {
                //        x_spike(fire) = s_update_stochastic_input(dt1, x_S_old(fire), x_S_old(fire + M), alpha, beta);
                s_update_stochastic_input( dt1[ id1 ], x_S_old[id1], x_S_old[id1 + M], alpha, beta, x_spike[id1] );
                //        x_spike(fire + M) = u_update_stochastic_input(dt1, x_S_old(fire), x_S_old(fire + M), alpha, beta) + 1;
                u_update_stochastic_input( dt1[ id1 ], x_S_old[id1], x_S_old[id1 + M], alpha, beta, x_spike[ id1+M ] );
                x_spike[ id1+M ] += 1;
            }
        }
        
        for( int id1 = 0; id1 < M; id1++ ) {
            //        vector = zeros(M,1);
            vector[ id1 ] = 0;
            if( fire[id1] == 1 ) {
                //        vector(fire) = dt1;
                vector[ id1 ] = dt1[ id1 ];
            }
        }
        
        for( int id1 = 0; id1 < M; id1++ ) {
            //        k1_w = 1/p.tau_w * (w_ss(v_S_old) - w_S_old);
            k1_w[ id1 ] = 1. / p.tau_w * ( w_ss( v_S_old[ id1 ] ) - w_S_old[ id1 ] );
            //        k2_w = 1/p.tau_w * (w_ss(v_S_old + k1_v .* vector) - (w_S_old + k1_w .* vector));
            k2_w[ id1 ] = 1. / p.tau_w * ( w_ss( v_S_old[ id1 ] + k1_v[ id1 ] * vector[ id1] ) - ( w_S_old[ id1 ] + k1_w[ id1 ] * vector[ id1 ] ) );
            
            //        w_spike = w_S_old + vector/2 .* (k1_w + k2_w) + (w_S_old - 1) * (exp(-1/p.tau_w_reset) - 1);
            w_spike[ id1 ] = w_S_old[id1] + vector[ id1 ] * 0.5 * ( k1_w[id1] + k2_w[id1] ) + ( w_S_old[ id1 ] - 1. ) * ( exp( -1. / p.tau_w_reset ) - 1. );
        }
        
        for( int id1 = 0; id1 < M; id1++ ) {
            //        vector = zeros(M,1);
            vector[ id1 ] = 0;
            //        vector(fire) = 1;
            if( fire[id1] == 1) {
                vector[ id1 ] = 1;
            }
            //        w_spike = w_spike .* vector;
            w_spike[ id1 ] = w_spike[ id1 ] * vector[ id1 ];
        }
        
        //        % calculate new exact value of n after reset dynamics
        for( int id1 = 0; id1 < M; id1++ ) {
            //        n_S(fire) = n_spike(fire) .* exp(-dt2/p.tau_n);
            if( fire[id1] == 1 ) {
                n_S[id1] = n_spike[id1] * exp( -dt2[id1] / p.tau_n );
            }
        }
        
        //        k1_v = driving_function_stochastic_input(p.v_reset * ones(M,1), n_spike, x_spike(1:M), w_spike, p, W, x_out_old(1:M), x_PN_old, I_stellate, ref_S, A, A_CNS_to_STELLATE, N);
        driving_function_stochastic_input( v_reset, n_spike, x_spike, w_spike, p, W, x_out_old, x_PN_old, I_stellate, ref_S, A, A_CNS_to_STELLATE, N, k1_v );
        
        for( int id1 = 0; id1 < M; id1++ ) {
            //        k1_w = 1/p.tau_w * (w_ss(p.v_reset*ones(M,1)) - w_spike);
            k1_w[ id1 ] = 1. / p.tau_w * ( w_ss( p.v_reset ) - w_spike[ id1 ] );
            //        vector = zeros(M,1);
            vector[ id1 ] = 0;
            //        vector(fire) = dt2;
            if( fire[ id1 ] == 1 ) {
                vector[ id1 ] = dt2[id1];
            }
        }
        //        k2_v = driving_function_stochastic_input(p.v_reset * ones(M,1) + k1_v .* vector, n_S, x_S(1:M), w_spike + k1_w .* vector, p, W, x_out(1:M), x_PN, I_stellate, ref_S, A,A_CNS_to_STELLATE, N);
        for( int id1 = 0; id1 < M; id1++ ) {
            vk1v[ id1 ] = p.v_reset + k1_v[ id1 ] * vector[ id1 ];
            wk1w[ id1 ] = w_spike[ id1 ] + k1_w[ id1 ] * vector[ id1 ];
        }
        driving_function_stochastic_input( vk1v, n_S, x_S, wk1w, p, W, x_out, x_PN, I_stellate, ref_S, A,A_CNS_to_STELLATE, N, k2_v);
        
        //        k2_w = 1/p.tau_w * (w_ss(p.v_reset*ones(M,1) + k1_v .* vector) - (w_spike + k1_w .* vector));
        for( int id1 = 0; id1 < M; id1++ ) {
            k2_w[ id1 ] = 1. / p.tau_w * ( w_ss( p.v_reset + k1_v[id1] * vector[id1] ) - ( w_spike[id1] + k1_w[id1] * vector[id1] ) );
            
            if( fire[id1] == 1 ) {
                //        v_S(fire) = p.v_reset + dt2/2 .* (k1_v(fire) + k2_v(fire));
                v_S[ id1 ] = p.v_reset + dt2[id1] * 0.5 * ( k1_v[id1] + k2_v[id1] );
                //        w_S(fire) = w_spike(fire) + dt2/2 .* (k1_w(fire) + k2_w(fire));
                w_S[ id1 ] = w_spike[id1] + dt2[id1] * 0.5 * ( k1_w[id1] + k2_w[id1] );
                
                //        % update s and u after reset dynamics
                //        x_S(fire) = s_update_stochastic_input(dt2, x_spike(fire), x_spike(fire + M), alpha, beta);
                s_update_stochastic_input( dt2[id1], x_spike[id1], x_spike[id1 + M], alpha, beta, x_S[id1] );
                //        x_S(fire + M) = u_update_stochastic_input(dt2, x_spike(fire), x_spike(fire + M), alpha, beta);
                u_update_stochastic_input( dt2[id1], x_spike[id1], x_spike[id1+M], alpha, beta, x_S[id1+M] );
                
                //        ISO(fire, m + 1) = ISO(fire,m) + param.dISO;
                ISO[ id1 ] = ISO[ id1 ] + param.dISO;
            }
        }
        //        end
        //        %     %%%%%%%%%%% Potassium ACh Channel Activation update %%%%%%
        //        %     % Concentration of ACh binding to muscarinic receptors.. Note that only
        //        %     % the PN cells are releasing ACh onto the SAN in this model.
        //        % ACh(m+1) = sum(x_PN(1:N.PN));% use nM in Behar-Yaniv model % * 10^(-6); to mM
        //        %     % time constants of activation of KACh current
        //        %     % rate of closing
        
        a = 9.96;
        //        % rate of opening
        AChxPN = 0;
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            AChxPN += x_PN[ id1 ];
        }
        AChxPN = ACh[ m%50 ];
        
        b = 12.32 * AChxPN / ( AChxPN + 4.2E-6 );
        
        K_ACh_activation = b / ( a + b );
        
        for( int id1 = 0; id1 < N.PN; id1++ ) {
            //        timer_PN = timer_PN - dt;
            timer_PN[ id1 ] -= dt;
            //        timer_PN(timer_PN <= 0) = 0;
            if( timer_PN[ id1 ] <= 0 ) {
                timer_PN[ id1 ] = 0;
            }
        }
        for( int id1 = 0; id1 < M; id1++ ) {
            //        timer_S = timer_S - dt;
            timer_S[ id1 ] -= dt;
            //        timer_S(timer_S <= 0) = 0;
            if( timer_S[ id1 ] <= 0 ) {
                timer_S[ id1 ] = 0;
            }
        }
        
        if( fmod(t,dt) == 0 ) {
            sumACh = 0;
            sumISO = 0;
            for( int id1 = 0; id1 < N.PN; id1++ ) {
                sumACh += ACh[id1];
            }
            for( int id1 = 0; id1 < M; id1++ ) {
                sumISO += ISO[id1];
            }
            fprintf( output, "%-16.14g\t%-10.8g\t%-10.8g\n",t, sumACh, sumISO);
        }
        
        if( fmod(t,100) == 0 ) {
            cout << t * 0.001 << "sec\t" << fire_size << endl;

        }
    }
    fclose( output );
}


// %%
// function y = intensity_function_vector(t, tone, P_ICNS, P_STELLATE, alpha, N, K)
void intensity_function_vector( double t, structtone &tone, structP_ICNS &P_ICNS, structP_STELLATE &P_STELLATE, double alpha, structN &N, int K, double * y ) {
    int M = N.CNS + N.STELLATE + N.ICNS;
    double t1 = fmod( t, 1000. );
    // y = zeros(M + K,1);
    using namespace std;
    
    //% STELLATE intensity function
    //for j = 1 : N.CNS + N.STELLATE
    double yj;
    yj = tone.STELLATE;
    if( t1 <= 50. ) {
        yj = yj +  ( ( P_STELLATE.max - P_STELLATE.min ) / 50. * t1 + P_STELLATE.min );
    } else if( t1 > 50 & t1 <= 250. ) {
        yj = yj + ( P_STELLATE.max * ( alpha - 1. ) / 200. * ( t1 - 50. ) + P_STELLATE.max );
    } else if ( t1 > 250 & t1 <= 500. ) {
        yj = yj + ( ( P_STELLATE.min - alpha * P_STELLATE.max ) / 250. * ( t1 - 250. ) + alpha * P_STELLATE.max );
    } else {
        yj = yj + P_STELLATE.min;
    }
    
    
    
    for( int j = 0; j < ( N.CNS + N.STELLATE ); j++ ) {
        //y(j) = ((P_STELLATE.max - P_STELLATE.min)/50 .* mod(t,1000) + P_STELLATE.min) .* (mod(t,1000) <= 50)...
        //+ (P_STELLATE.max * (alpha - 1)/200 .* (mod(t,1000) - 50) + P_STELLATE.max) .* (mod(t,1000) > 50 & mod(t,1000) <= 250) ...
        //+ ((P_STELLATE.min - alpha * P_STELLATE.max)/250 .* (mod(t,1000) - 250) + alpha * P_STELLATE.max) .*( mod(t,1000) > 250 & mod(t,1000) <= 500)...
        //+P_STELLATE.min .* (mod(t,1000) > 500) + tone.STELLATE;
        y[ j ] = yj;
        //end
    }
    //% SICNS intensity function
    // for j = N.CNS + N.STELLATE + 1 : M
    yj = tone.SICNS;
    if( t1 <= 50. ) {
        yj += (( P_ICNS.max - P_ICNS.min ) / 50. * t1 + P_ICNS.min );
    } else if( t1 > 50 & t1 <= 250. ) {
        yj += ( P_ICNS.max * ( alpha - 1. ) / 200. * ( t1 - 50 ) + P_ICNS.max ) ;
    } else if( t1 > 250 & t1 <= 500 ) {
        yj +=  ( P_ICNS.min - alpha * P_ICNS.max ) / 250. * ( t1 - 250 ) + alpha * P_ICNS.max ;
    } else {
        yj += P_ICNS.min;
    }
    
    // cout << yj << "\t";
    
    for( int j = ( N.CNS + N.STELLATE ); j < M ; j++ ) {
        //%   y(j) = 0; % No Sympathetic ICNS
        // y(j) = ((P_ICNS.max - P_ICNS.min)/50 .* mod(t,1000) + P_ICNS.min) .* (mod(t,1000) <= 50)...
        // + (P_ICNS.max * (alpha - 1)/200 .* (mod(t,1000) - 50) + P_ICNS.max) .* (mod(t,1000) > 50 & mod(t,1000) <= 250) ...
        // + ((P_ICNS.min - alpha * P_ICNS.max)/250 .* (mod(t,1000) - 250) + alpha * P_ICNS.max) .*( mod(t,1000) > 250 & mod(t,1000) <= 500)...
        // +P_ICNS.min .* (mod(t,1000) > 500) + tone.SICNS;
        y[ j ] = yj;
        //end
    }
    
    // % PICNS intensity function
    // for j = M + 1 : M + K
    //%     y(j) = 0; % No Parasympathetic ICNS
    yj = tone.PICNS;
    if( t1 <= 50. ) {
        yj += (( P_ICNS.max - P_ICNS.min ) / 50. * t1 + P_ICNS.min );
    } else if( t1 > 50 & t1 <= 250. ) {
        yj += ( P_ICNS.max * ( alpha - 1. ) / 200. * ( t1 - 50. ) + P_ICNS.max ) ;
    } else if( t1 > 250 & t1 <= 500. ) {
        yj += ( P_ICNS.min - alpha * P_ICNS.max ) / 250. * ( t1 - 250. ) + alpha * P_ICNS.max ;
    } else {
        yj += P_ICNS.min;
    }
    
    // cout << yj << "\n";
    
    for ( int j = M; j < M + K; j++ ) {
        // y(j) = ((P_ICNS.max - P_ICNS.min)/50 .* mod(t,1000) + P_ICNS.min) .* (mod(t,1000) <= 50)...
        // + (P_ICNS.max * (alpha - 1)/200 .* (mod(t,1000) - 50) + P_ICNS.max) .* (mod(t,1000) > 50 & mod(t,1000) <= 250) ...
        // + ((P_ICNS.min - alpha * P_ICNS.max)/250 .* (mod(t,1000) - 250) + alpha * P_ICNS.max) .*( mod(t,1000) > 250 & mod(t,1000) <= 500)...
        // +P_ICNS.min .* (mod(t,1000) > 500) + tone.PICNS;
        y[ j ] = yj;
        //end
    }
    //end
    
}
