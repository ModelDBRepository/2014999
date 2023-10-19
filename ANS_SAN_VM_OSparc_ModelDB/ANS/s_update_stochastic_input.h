//% here t is the time to update s by,
//% i.e. t = t_(k+1) - t_k
//function s = s_update_stochastic_input(t,s,u,alpha,beta)
//s = 1/(alpha - beta) * (alpha * exp(-beta * t) - beta * exp(-alpha * t)) .* s + ...
//    1/(alpha - beta) * (exp(-beta * t) - exp(-alpha * t)) .* u;
//end

void s_update_stochastic_input(double t, double s, double u, double alpha, double beta, double &s1);

void s_update_stochastic_input(double t, double s, double u, double alpha, double beta, double &s1) {

    s1 = ( 1. / ( alpha - beta ) * ( alpha * exp(-beta * t) - beta * exp(-alpha * t) ) * s
          + 1. / ( alpha - beta ) * (exp(-beta * t) - exp(-alpha * t)) * u );

}
