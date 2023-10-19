void Erdos_Renyi( int n, double p, double **A , std::default_random_engine &generator );

//#include <fstream>
//#include <iostream>
//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#include <string.h>
//#include <random>
//function A = Erdos_Renyi(n,p)

void Erdos_Renyi( int n, double p, double **A, std::default_random_engine &generator ){
    // std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> urand(0.0,1.0);
    using namespace std;
    
    //A = zeros(n,n);
    //for j = 1 : n
    //    for k = 1 : n
    //        prob = rand(1);
    //        if prob < p
    //            A(j,k) = 1;
    //        end
    //    end
    //end
    //A = A - diag(diag(A));
    //end
    double prob;
    
    for( int j = 0; j < n; j++ ) {
        for( int k = 0; k < n; k++ ) {
            A[j][k] = 0;
            prob = urand( generator );
            if( prob < p ) {
                A[j][k] = 1;
            }
            // cout << A[j][k] << "\t";
        }
        // cout << "\n";
    }
    for( int j = 0; j < n; j++ ) {
        A[j][j] = 0;
    }
}
