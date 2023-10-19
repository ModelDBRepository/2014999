#include <complex>
#include <cmath>

double update_PKA( double x0 );

double update_PKA( double x0 ) {
//function PKA = update_PKA(x)
//% x: cAMP concentration (pmol/mg protein)
//
//PKA =.33333333333333333333333333333333e-4./(300.+15.*x+x.^2).*...
//    (-468891428493892500.*x.^2-46425728249607375.*...
//    x.^3-14479420376562975.*x.^4-666054712166085.*...
//    x.^5-47282992935913.*x.^6-741463359651000000.-111219503947650000.*...
//    x+9000000.*x.*(-1643644331734924815.*x.^6-19572532067205238230.*...
//    x.^5-393461126229536661375.*x.^4-80438302732946971500.*...
//    x.^3-805328393113024740000.*x.^2-22688778805320600000.*...
//    x-151258525368804000000.).^(1/2)+450000.*x.^2.*(-1643644331734924815.*...
//    x.^6-19572532067205238230.*x.^5-393461126229536661375.*...
//    x.^4-80438302732946971500.*x.^3-805328393113024740000.*...
//    x.^2-22688778805320600000.*x-151258525368804000000.).^(1/2)+30000.*...
//    x.^3.*(-1643644331734924815.*x.^6-19572532067205238230.*...
//    x.^5-393461126229536661375.*x.^4-80438302732946971500.*...
//    x.^3-805328393113024740000.*x.^2-22688778805320600000.*...
//    x-151258525368804000000.).^(1/2)).^(1/3)+.33333333333333333333333333333333e-4.*...
//    (349407388425.*x.^2+17367968670.*x.^3+1548762289.*x.^4+819206010000.+81920601000.*x)...
//    ./(300.+15.*x+x.^2)./(-468891428493892500.*x.^2-46425728249607375.*...
//    x.^3-14479420376562975.*x.^4-666054712166085.*x.^5-47282992935913.*...
//    x.^6-741463359651000000.-111219503947650000.*x+9000000.*x.*...
//    (-1643644331734924815.*x.^6-19572532067205238230.*x.^5-393461126229536661375.*...
//    x.^4-80438302732946971500.*x.^3-805328393113024740000.*...
//    x.^2-22688778805320600000.*x-151258525368804000000.).^(1/2)+450000.*...
//    x.^2.*(-1643644331734924815.*x.^6-19572532067205238230.*x.^5-393461126229536661375.*...
//    x.^4-80438302732946971500.*x.^3-805328393113024740000.*x.^2-22688778805320600000.*...
//    x-151258525368804000000.).^(1/2)+30000.*x.^3.*(-1643644331734924815.*...
//    x.^6-19572532067205238230.*x.^5-393461126229536661375.*x.^4-80438302732946971500.*...
//    x.^3-805328393113024740000.*x.^2-22688778805320600000.*x-151258525368804000000.).^(1/2)).^(1/3)...
//    -.33333333333333333333333333333333e-4.*(905100.+45255.*x+23017.*x.^2)./(300.+15.*x+x.^2);
//
//    PKA = abs(PKA);

    std::complex<double> x = x0;
    
    std::complex<double> mycomplex = -(x*1.5085+(x*x)*7.672333333333333E-1+3.017E1)/(x*1.5E1+x*x+3.0E2)+(pow(x*-1.1121950394765E17+(x*x)*sqrt(x*-2.26887788053206E19-(x*x)*8.053283931130247E20-(x*x*x)*8.043830273294696E19-(x*x*x*x)*3.934611262295367E20-(x*x*x*x*x)*1.957253206720524E19-(x*x*x*x*x*x)*1.643644331734925E18-1.51258525368804E20)*4.5E5+(x*x*x)*sqrt(x*-2.26887788053206E19-(x*x)*8.053283931130247E20-(x*x*x)*8.043830273294696E19-(x*x*x*x)*3.934611262295367E20-(x*x*x*x*x)*1.957253206720524E19-(x*x*x*x*x*x)*1.643644331734925E18-1.51258525368804E20)*3.0E4+x*sqrt(x*-2.26887788053206E19-(x*x)*8.053283931130247E20-(x*x*x)*8.043830273294696E19-(x*x*x*x)*3.934611262295367E20-(x*x*x*x*x)*1.957253206720524E19-(x*x*x*x*x*x)*1.643644331734925E18-1.51258525368804E20)*9.0E6-(x*x)*4.688914284938925E17-(x*x*x)*4.642572824960738E16-(x*x*x*x)*1.447942037656298E16-(x*x*x*x*x)*6.66054712166085E14-(x*x*x*x*x*x)*4.7282992935913E13-7.414633596510001E17,1.0/3.0)*3.333333333333333E-5)/(x*1.5E1+x*x+3.0E2)+((x*2.7306867E6+(x*x)*1.16469129475E7+(x*x*x)*5.78932289E5+(x*x*x*x)*5.162540963333333E4+2.7306867E7)*1.0/pow(x*-1.1121950394765E17+(x*x)*sqrt(x*-2.26887788053206E19-(x*x)*8.053283931130247E20-(x*x*x)*8.043830273294696E19-(x*x*x*x)*3.934611262295367E20-(x*x*x*x*x)*1.957253206720524E19-(x*x*x*x*x*x)*1.643644331734925E18-1.51258525368804E20)*4.5E5+(x*x*x)*sqrt(x*-2.26887788053206E19-(x*x)*8.053283931130247E20-(x*x*x)*8.043830273294696E19-(x*x*x*x)*3.934611262295367E20-(x*x*x*x*x)*1.957253206720524E19-(x*x*x*x*x*x)*1.643644331734925E18-1.51258525368804E20)*3.0E4+x*sqrt(x*-2.26887788053206E19-(x*x)*8.053283931130247E20-(x*x*x)*8.043830273294696E19-(x*x*x*x)*3.934611262295367E20-(x*x*x*x*x)*1.957253206720524E19-(x*x*x*x*x*x)*1.643644331734925E18-1.51258525368804E20)*9.0E6-(x*x)*4.688914284938925E17-(x*x*x)*4.642572824960738E16-(x*x*x*x)*1.447942037656298E16-(x*x*x*x*x)*6.66054712166085E14-(x*x*x*x*x*x)*4.7282992935913E13-7.414633596510001E17,1.0/3.0))/(x*1.5E1+x*x+3.0E2);
    
    //cout << x << "\t" << mycomplex << endl;
    
    double PKA = std::abs( mycomplex );
    
    return PKA;

}

