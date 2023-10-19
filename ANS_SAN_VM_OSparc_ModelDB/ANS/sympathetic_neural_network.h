//% Build a random network with two levels - ICNS and Stellate ganglia. 
//
//% n determines how many "weak" synaptic couplings there are in the
//% efferent connections between adjacent layers. We assume this is uniform across the
//% layers.
//
//% weights is a structure which contains weights.strong (n x 1) and
//% weights.weak (1 x 1)
//% These determine the strength of the efferent coupling betweem descending
//% layers.
//
//% N is a structure which contains N.CNS, N.STELLATE, and N.ICNS. These give
//% the number of neurons in each layer. NOTE
//
//% p is a structure which contains p.CNS, P.STELLATE, p.CNS,
//% p.effferent_STELLATe, and p.efferent_ICNS. p.x determines the probabilty
//% of neurons being coupled intra-layer, and p.efferent_x determines the
//% probability that a neuron receives effernt input from the layer
//% immediately above.
//
//
//function [A, A_CNS_to_STELLATE] = sympathetic_neural_network(N,weights,n,p)

void sympathetic_neural_network(structN &N, structweights &weights, int n, structp &p, structM3PN &A, structM3PN &A_CNS_to_STELLATE, std::default_random_engine &generator);

void randperm( int nn, int k, int *C , std::default_random_engine &generator ) {
    std::uniform_real_distribution<double> urand(0.0,1.0);
    
    int id1, id2, a, b;
    int B[nn];
    for( id1 = 0; id1 < nn; id1++ ){
        B[id1] = id1;
    }
    
    for( id1 = 0; id1 < k; id1++ ){
        id2 = floor( urand( generator ) * nn );
        a = B[ id2 ];
        B[ id2 ] = B[ id1 ];
        B[ id1 ] = a;
    }
    
    for( id1 = 0; id1 < k; id1++ ) {
        C[ id1 ] = B[ id1 ];
    }
}

void sympathetic_neural_network(structN &N, structweights &weights, int n, structp &p, structM3PN &A, structM3PN &A_CNS_to_STELLATE, std::default_random_engine &generator){
    using namespace std;
    //cout.precision(2);
    
    //std::default_random_engine generator(time(0));

//
//
//% Set up total network
//A = zeros(N.ICNS + N.STELLATE + N.CNS);
    
//
//% Set up sub-networks
//A_CNS = Erdos_Renyi(N.CNS, p.CNS);
    double **A_CNS;
    A_CNS = new double *[N.CNS];
    for( int id1=0; id1<N.CNS; id1++ ) {
        A_CNS[id1] = new double[N.CNS];
    }
    Erdos_Renyi( N.CNS, p.CNS, A_CNS , generator );
    
    
//A_STELLATE = Erdos_Renyi(N.STELLATE, p.STELLATE);
    double **A_STELLATE;
    A_STELLATE = new double *[N.STELLATE];
    for( int id1=0; id1<N.STELLATE; id1++ ) {
        A_STELLATE[id1] = new double[N.STELLATE];
    }
    Erdos_Renyi( N.STELLATE, p.STELLATE, A_STELLATE, generator );
    
//A_ICNS = Erdos_Renyi(N.ICNS, p.ICNS);
    double **A_ICNS;
    A_ICNS = new double *[N.ICNS];
    for( int id1=0; id1<N.ICNS; id1++ ) {
        A_ICNS[id1] = new double[N.ICNS];
    }
    Erdos_Renyi( N.ICNS, p.ICNS, A_ICNS, generator );

//A(1:N.CNS, 1:N.CNS) = A_CNS;
    for( int id1=0; id1<N.CNS; id1++ ) {
        for( int id2=0; id2<N.CNS; id2++ ) {
            A.M[id1][id2] = A_CNS[id1][id2];
            //cout << A.M[id1][id2] << "  ";
        }
        //cout << endl << endl;
    }

//A(N.CNS + 1 : N.CNS + N.STELLATE, N.CNS + 1 : N.CNS + N.STELLATE) = A_STELLATE;
    int sum1 = N.CNS + N.STELLATE;
    for( int id1 = N.CNS; id1 < sum1; id1++ ) {
        for( int id2 = N.CNS; id2 < sum1; id2++ ) {
            A.M[id1][id2] = A_STELLATE[id1-N.CNS][id2-N.CNS];
        }
    }

    //A(N.CNS + N.STELLATE + 1 : end, N.CNS + N.STELLATE + 1 : end) = A_ICNS;
    int sum2 = sum1 + N.ICNS;
    for( int id1 = sum1; id1 < sum2; id1++ ) {
        for( int id2 = sum1; id2 < sum2; id2++ ) {
            A.M[id1][id2] = A_ICNS[id1-sum1][id2-sum1];
        }
    }

    // std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> urand(0.0,1.0);

//% Set n+1 weak coupling from CNS to STELLATE and from STELLATE to ICNS
//% Generate neurons that receive efferent input
//idx_STELLATE = rand(N.STELLATE, 1);
    double rand_STELLATE;
    int idx_STELLATE[N.STELLATE];
    int length_idx_STELLATE = 0;
    for( int id1 = 0; id1<N.STELLATE; id1++ ) {
        rand_STELLATE = urand( generator );
        // cout << rand_STELLATE << "\t" << p.efferent_STELLATE << endl;
        if( rand_STELLATE < p.efferent_STELLATE ){
            idx_STELLATE[length_idx_STELLATE] = id1;
            length_idx_STELLATE += 1;
            // cout << length_idx_STELLATE << endl;

        }
    }
    
//idx_ICNS = rand(N.ICNS, 1);
    double rand_ICNS;
    int idx_ICNS[N.ICNS];
    int length_idx_ICNS = 0;
    for( int id1 = 0; id1<N.ICNS; id1++ ) {
        rand_ICNS = urand( generator );
        if( rand_ICNS < p.efferent_ICNS ) {
            idx_ICNS[length_idx_ICNS] = id1;
            length_idx_ICNS += 1;
        }
    }
    
//% Both layers have probability p.efferent_STELLATE/ICNS of receiving efferent input from
//% the layer immediatly above them i.e. CNS to STELLATE, STELLATE to ICNS.
//% Here we determine the index of the neurons in each layer that receive efferent input
//% from the layer immediately above.
//
//idx_STELLATE = find(idx_STELLATE < p.efferent_STELLATE);
//idx_ICNS = find(idx_ICNS < p.efferent_ICNS);
//
//if (isempty(idx_STELLATE) || isempty(idx_ICNS))
//    fprintf('No intra-layer connections, rebuild network \n')
//    return
//end
    
    if( length_idx_ICNS == 0 || length_idx_STELLATE == 0 ) {
        cout << "No intra-layer connections, rebuild network \n";
        return;
    }

    
//% Here we determine the neurons in the layer immediately above which synapse
//% onto the STELLATE and ICNS.
//
//% initalizing values
//CNS_to_STELLATE = zeros(n + 1, length(idx_STELLATE));
    int **CNS_to_STELLATE;
    CNS_to_STELLATE = new int *[n+1];
    
//STELLATE_to_ICNS = zeros(n + 1,  length(idx_ICNS));
    int **STELLATE_to_ICNS;
    STELLATE_to_ICNS = new int *[n+1];
    
    for( int id1 = 0; id1 < (n+1); id1++ ) {
        CNS_to_STELLATE[id1] = new int [N.STELLATE];
        STELLATE_to_ICNS[id1] = new int [N.ICNS];
    }
    
//
//% n + 1 rows corresponding to n + 1 connections,  idx_x columns corresponding to
//% number of neurons in layer x which receive efferent input.
//A_CNS_to_STELLATE = zeros(N.CNS, N.STELLATE);
    for( int id1 = 0; id1 < N.CNS; id1++ ) {
        for( int id2 = 0; id2 < N.STELLATE; id2++ ) {
            A_CNS_to_STELLATE.M[id1][id2] = 0.;
        }
    }
//
    int Vperm[n+1];
//for j = 1 : length(idx_STELLATE)
    for( int j = 0; j < length_idx_STELLATE; j++ ) {
//    % choose the nuerons in the CNS which synapse on the STELLATE neurons
//    % which receives the n + 1 synaptic currents
//    CNS_to_STELLATE(:,j) = randperm(N.CNS, n + 1);
        randperm( N.CNS, (n+1), Vperm, generator );
        for( int id1 = 0; id1 < (n+1); id1++ ) {
            CNS_to_STELLATE[id1][j] = Vperm[ id1 ];
        }
//    A_CNS_to_STELLATE(CNS_to_STELLATE(1:(end-1),j), idx_STELLATE(j)) = weights.weak/2;
        for( int id1 = 0; id1 < n; id1++ ) {
            A_CNS_to_STELLATE.M[ CNS_to_STELLATE[id1][j] ][ idx_STELLATE[j] ] = weights.weak * 0.5;
        }
//    A_CNS_to_STELLATE(CNS_to_STELLATE(end, j), idx_STELLATE(j)) = weights.strong/2;
        A_CNS_to_STELLATE.M[ CNS_to_STELLATE[n][j] ][ idx_STELLATE[j] ] = weights.strong * 0.5;
//end
    }
    
    for( int j = 0; j < length_idx_ICNS; j++ ) {
//for j = 1 : length(idx_ICNS)
//    % choose the nuerons in the STELLATE which synapse on the ICNS neurons
//    % which receives the n + 1 synaptic currents
//    STELLATE_to_ICNS(:,j) = randperm(N.STELLATE, n + 1);
        randperm( N.STELLATE, (n+1), Vperm, generator );
        for( int id1 = 0; id1 < (n+1); id1++ ) {
            STELLATE_to_ICNS[id1][j] = Vperm[id1];
        }
//    A(N.CNS + STELLATE_to_ICNS(1:(end-1),j), N.CNS + N.STELLATE + idx_ICNS(j)) = weights.weak/2;
        for( int id1 = 0; id1 < n; id1++ ) {
            A.M[(N.CNS + STELLATE_to_ICNS[id1][j])][(N.CNS + N.STELLATE + idx_ICNS[j])] = weights.weak * 0.5;
        }
        
        //????????????????????????????????
        //  weights.strong?
        //????????????????????????????????
//    A(N.CNS + STELLATE_to_ICNS(end,j), N.CNS + N.STELLATE + idx_ICNS(j)) = weights.weak/2;
        A.M[(N.CNS + STELLATE_to_ICNS[n][j])][(N.CNS + N.STELLATE + idx_ICNS[j])] = weights.weak * 0.5;
//end
    }
//
//
//% set up afferent input to STELLATE from ICNS
//idx_STELLATE = rand(N.STELLATE, 1);

//idx_STELLATE = find(idx_STELLATE < p.afferent_STELLATE);
    // double rand_STELLATE;
    // int idx_STELLATE[N.STELLATE];
    length_idx_STELLATE = 0;
    for( int id1 = 0; id1<N.STELLATE; id1++ ) {
        rand_STELLATE = urand( generator );
        if( rand_STELLATE < p.efferent_STELLATE ) {
            idx_STELLATE[length_idx_STELLATE] = id1;
            length_idx_STELLATE += 1;
        }
    }
    
//ICNS_to_STELLATE = zeros(1, length(idx_STELLATE));
    int ICNS_to_STELLATE[ N.STELLATE ];
//for j = 1 : length(idx_STELLATE)
//    ICNS_to_STELLATE(1,j) = randperm(N.ICNS, 1);
//    A(N.CNS + N.STELLATE + ICNS_to_STELLATE(1,j), N.CNS + idx_STELLATE(j)) = 1;
//end
    for( int j = 0; j < length_idx_STELLATE; j++ ) {
        ICNS_to_STELLATE[ j ] = floor( urand( generator ) * N.ICNS ) ;
        A.M[(N.CNS + N.STELLATE + ICNS_to_STELLATE[j])][(N.CNS + idx_STELLATE[j])] = 1;
        
    }
}
