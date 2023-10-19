//function A = couple_two_networks(N1, N2, p)
void couple_two_networks( int N1, int N2, double p, double ** A, std::default_random_engine &generator );
void couple_two_networks( int N1, int N2, double p, double ** A, std::default_random_engine &generator ){
    // std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> urand(0.0,1.0);
    
//% here N1 and N2 are the respective sizes of the two networks, and p is the
//% probability that a neuron in network 1 is connected to a neuron in
//% network 2. Note here that we are only generating a directed network FROM
//% network 1 TO network 2. This is for ease of coding.
//
//A = rand(N1, N2);
//A(A < p) = 1;
//A(A ~= 1) = 0;
//
//end
    double myrand;
    for( int id1 = 0; id1 < N1; id1++ ) {
        for( int id2 = 0; id2 < N2; id2++ ) {
            myrand = urand( generator );
            A[id1][id2] = 0;
            if( myrand < p ) {
                A[id1][id2] = 1;
            }
        }
    }
}
