#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <random>
#include <new>
#include <stdlib.h>
#include <cstdlib>

#include <unistd.h>


// for mkdir
#include <sys/types.h>
#include <sys/stat.h>


int main( ) {
    using namespace std;
    cout.precision(6);
    using namespace std;
    char filename[200];

    char parameters[100];
    double datasample;
    int datas_int;
    double datas_double;
    
    sprintf( filename, "%s", "ANS/stim_param.txt");
    FILE *output1;
    output1 = fopen( filename, "w");

    cout << "Enter simulation time length (sec):" ;
    cin >> datasample;
    fprintf( output1, "%-16.14g", datasample );
    
    fclose( output1 );
    

    sprintf( filename, "%s", "SAN/stim_param.txt");
    output1 = fopen( filename, "w");

    datas_int = 2;
    
    while( datas_int<0 || datas_int>1 ) {
        cout << "ICNS On or Off:\n" ;
        cout << "Enter \"1\" for On, \"0\" for Off.\n";
        cin >> datas_int;
        if( datas_int == 1 ) {
            cout << "Enter start time(sec): ";
            cin >> datas_double;
            fprintf( output1, "%-16.14g\n", datas_double );
            
            cout << "Enter end time(sec): ";
            cin>> datas_double;
            fprintf( output1, "%-16.14g\n", datas_double );
        } else {
            fprintf( output1, "%-d\n%-d\n", 0, 0 );
        }
    }
    
    datas_int = 2;
    
    while( datas_int<0 || datas_int>1 ) {
        cout << "ACh On or Off:\n" ;
        cout << "Enter \"1\" for On, \"0\" for Off.\n";
        cin >> datas_int;
        if( datas_int == 1 ) {
            cout << "Enter start time(sec): ";
            cin >> datas_double;
            fprintf( output1, "%-16.14g\n", datas_double );
            
            cout << "Enter end time(sec): ";
            cin>> datas_double;
            fprintf( output1, "%-16.14g\n", datas_double );
        } else {
            fprintf( output1, "%-d\n%-d\n", 0, 0 );
        }
    }
    fclose( output1 );

    
    
    chdir("ANS");
    //system("ls");
    
    sprintf( filename, "%s", "Autonomic_Layered_Network_Model.cpp" );
    string str = "g++ -O2 ";
    str = str + " -o a.out " + filename;
    
    // Convert string to const char * as system requires // parameter of type const char *
    const char *command = str.c_str();
    //cout << "Compiling file using " << command << endl;
    system(command);
    

    system("./a.out stim_param.txt outputANSFolder");
    chdir("..");
    //system("ls");
    
    chdir("SAN");
    //system("ls");
    system("g++ -O2 -o a.out main_ICNS_isoach.cpp");
    system("./a.out ../ANS/outputANSFolder/outputfile.txt stim_param.txt");
    system("python3 plotall.py");
    chdir("..");
    //system("ls");
    
    
    chdir("VM");
    system("ls");
    system("g++ -O2 -o a.out soltis_biophysJ2010_masterCompute.cpp");
    system("./a.out ../SAN/cycles_icnsach.txt ../ANS/outputANSFolder/outputfile.txt ../ANS/stim_param.txt ../SAN/stim_param.txt");
    system("python3 plotall.py");
    chdir("..");
    //system("ls");
    
     
    
    return 0;
}
