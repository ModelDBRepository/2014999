//% This function updates the u state-variable
//function u = u_update_stochastic_input(t,s,u,alpha,beta)
//u = alpha * beta/(alpha - beta) * (exp(-alpha * t) - exp(-beta * t)) .* s + ...
//    1/(alpha - beta) * (alpha * exp(-alpha * t) - beta * exp(-beta * t)) .* u;
//end

void u_update_stochastic_input(double t, double s, double u, double alpha,double beta, double &s1 );
void u_update_stochastic_input(double t, double s, double u, double alpha,double beta, double &s1 ) {
    s1 = ( alpha * beta / ( alpha - beta) * (exp(-alpha * t) - exp(-beta * t)) * s
          + 1. / ( alpha - beta ) * ( alpha * exp( -alpha * t ) - beta * exp(-beta * t) ) * u );
}
