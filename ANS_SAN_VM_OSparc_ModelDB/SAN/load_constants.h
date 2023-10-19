
struct cell_con {
    double I_CaL;
    double I_NCX;
    double I_KACh;
    double I_f;
    
    
    double ISO;
    double CCh;
    int ALL_VAR;
    double cAB_Conc;
    double cAB_on_rate;
    double cAB_off_rate;
    
    double Cao; //  2; % [mM] extracellular Ca2+
    double Ko; //  5.4; % [mM] extracellular K+
    double Ki; //  140; % [mM] intracellular K+
    double Nao; //  140; % [mM] extracellular Na+
    double Nai; //  10; % [mM] intracellular Na+
    double Mgi; //  2.5; % [mM] intracellular Mg+
    
    // %% Cell compartments
    double C; //  32; % [pF] cell electric capacitance
    double L_cell; //  70; % [um] cell length
    double R_cell; //  4; % [um] cell radius
    double L_sub; //  0.02; % [um] distrance between JSR and surface membrane (submembrane space)
    double V_cell; //  pi*con.R_cell^2*con.L_cell*10^-3; % [pL] cell volume (3.5185838)
    double V_sub; //  2*pi*con.L_sub*(con.R_cell-con.L_sub/2)*con.L_cell*10^-3; % [pL] submembrane cell volume (0.035097874) -> 10^-3 to make it pL
    double V_jSR_part; //  0.0012; % [pL] part of cell volume occupied by junctional SR
    double V_jSR; //  con.V_jSR_part*con.V_cell; % [pL] volume of junctional SR (Ca2+ release store)
    double V_i_part; //  0.46; % part of the cell volume occupied with myoplasm
    double V_i; //  con.V_i_part*con.V_cell-con.V_sub; % [pL] myoplasmic volume
    double V_nSR_part; //  0.0116; % part of the cell volume occupied by network SR
    double V_nSR; //  con.V_nSR_part*con.V_cell; % [pL] volume of network SR (Ca2+ uptake store)
    double V_myto; //  0.6334; % [pL] mitochondrial volume
    
    // %% The Nernst equation
    double F; //  96485; % C/M is the Faraday constant
    double T; //  310.15; % K is the absolute emperature for 37*C
    double R; //  8.3144; % J/(M K) is the universal gaz constant
    double E_T; //  1000*(con.R*con.T/con.F); % [mV] is the 'RT/F' factor
    
    // %% Electric potential
    double E_Na; //  con.E_T*log(con.Nao/con.Nai); % [mV] equilibrium potential for Na+
    double E_K; //  con.E_T*log(con.Ko/con.Ki); % [mV] equilibrium potential for K+
    double E_Ks; //  con.E_T*log((con.Ko+0.12*con.Nao)/(con.Ki+0.12*con.Nai)); % [mV] reversal potential of I_Ks
    double E_CaL; //  45; % [mV] apparent reversal potential of I_CaL
    double E_CaT; //  45; % [mV] apparent reversal potential of I_CaT
    double E_st; //  37.4; % [mV] apparent reversal potential of I_st
    
    // %% Sarcolemmal Ion currents and their normalized conductances (g_x)
    double g_CaL; //  0.9; %0.5720; % [nS/pF] normalized conductance for I_CaL channels
    double g_CaT; //  0.1832; % [nS/pF] normalized conductance for I_CaT channels
    double g_If; //  0.455; %0.07; % [nS/pF] normalized conductance for I_f channels
    double g_st; //  0.00001; % [nS/pF] normalized conductance for I_st channels
    double g_Kr; //  0.1; % 0.08113973; % [nS/pF] normalized conductance for I_Kr channels
    double g_Ks; //  0.013; % [nS/pF] normalized conductance for I_Ks channels
    double g_to; //  0.252; % [nS/pF] normalized conductance for I_to channels
    double g_sus; //  0.02; % [nS/pF] normalized conductance for I_sus channels
    double g_bCa; //  0.0006; % [nS/pF] normalized conductance for I_bCa channels
    double g_bNa; //  0.00486; % [nS/pF] normalized conductance for I_bNa channels
    double g_KACh_max; //  2.2*0.14241818; % [nS/pF] normalized conductance for I_KACh channels
    double I_NaKmax; //  2.88; % [pA/pF] Maximal Na/K pump current conductance
    double k_NCX; //  225; % [pA/pF] Maximal Na+/Ca2+ exchanger current conductance
    double V_If12; //  -64; % -43.14; % [mV] Half activation voltage for I_f current in the basal state
    
    // %% Modulation of sarcolemmal ion current by ions
    double K_mfCa; //  0.00035; % [mM] dissociation constant of Ca2+ dependent I_CaL inactivation
    double K_mKp; //  1.4; % [mM] half-maximal Ko for I_NaK
    double K_mNap; //  14; % [mM] half-maximal Nai for I_NaK
    double beta_fCa; //  60; % [mM-1.ms-1] Ca2+ association rate constant for I_CaL
    double alpha_fCa; //  0.021; % [ms-1] Ca2+ dissociation rate constant for I_CaL
    
    // %% NCX function
    double K_1ni; //  395.3; % [mM] Nai binding to the first site on NCX
    double K_2ni; //  2.289; % [mM] Nai binding to second site on NCX
    double K_3ni; //  26.44; % [mM] Nai binding to third site on NCX
    double K_1no; //  1628; % [mM] Nao binding to first site on NCX
    double K_2no; //  561.4; % [mM] Nao binding to second site on NCX
    double K_3no; //  4.663; % [mM] Nao binding to third site on NCX
    double K_ci; //  0.0207; % [mM] Cai binding on NCX transporter
    double K_co; //  3.663; % [mM] Cao binding to NCX transporter
    double K_cni; //  26.44; % [mM] Nai and Cai simultaneous binding to NCX
    double Q_ci; //  0.1369; % Cai occlusion reaction of NCX
    double Q_co; //  0; % Cao occlusion reaction of NCX
    double Q_n; //  0.4315; % Na occlusion of NCX
    
    // %% Ca2+ diffusion
    double tho_difCa; //  0.04; % [ms] Time constant of Ca diffusion from the submembrane to myoplasm
    double tho_tr; //  40; % [ms] Time constant for Ca Ca transport from the network to junctional SR
    
    // %% SR Ca2+ ATPase function
    double K_up; //  0.6*10^-3; % Half-maximal Cai for Ca uptake in the network SR [mM]
    double P_up; //  0.8*0.012; % Rate constant for Ca uptake by the Ca pump in the network SR [mM/ms]
    double k_m_up; //  0.01;
    double k_i_up_ADP; //  5.1; % [mM]
    
    // %% RyR function
    double k_oCa_max; //  10; % [1/(mM^2*ms^1)]
    double k_om; //  0.06; % [1/ms]
    double k_iCa; //  0.5; % [1/(mM*ms)]
    double k_im; //  0.005; % [1/ms]
    double EC_50SR; //  0.45; % [mM]
    double k_s; //  400e3;
    double MaxSR; //  13;
    double MinSR; //  1;
    double HSR; //  3;
    // % RyR phosphorylation constants
    double RyR_min; //  0.0127; % derived from PP1 activity
    double RyR_max; //  0.02;
    double n_RyR; //  9.773;
    double k_05Ry; //  0.7;
    
    // %% Ca2+ and Mg2+ buffering
    double k_bCM; //  0.542; % [1/ms] Ca dissociation constant for calmodulin
    double k_bCQ; //  0.445; % [1/ms] Ca dissociation constant for calsequestrin
    double k_bTC; //  0.446; % [1/ms] Ca dissociation constant for the troponin-Ca site
    double k_bTMC; //  0.00751; % [1/ms] Ca dissociation constant for the troponin-Mg site
    double k_bTMM; //  0.751; % [1/ms] Mg dissociation constant for the troponin-Mg site
    double k_fCM; //  227.7; % [1/mM*ms] Ca association constant for calmodulin
    double k_fCQ; //  0.534; % [1/mM*ms] Ca association constant for calsequestrin
    double k_fTC; //  88.8; % [mm/ms] Ca association constant for troponin
    double k_fTMC; //  227.7; % [mM/ms] Ca association constant for the troponin-Mg site
    double k_fTMM; //  2.277; % [mM/ms] Mg association constant for the troponin-Mg site
    double TC_tot; //  0.042; % [mM] Total concentration of the troponin-Ca site
    // %con.TMC_tot; //  0.062; % [mM] Total concentration of the troponin-Mg site
    double CQ_tot; //  10; % [mM] Total calsequestrin concentration
    double CM_tot; //  0.045; % [mM] Total calmodium concentration
    
    // %% cAMP equations constants
    // double k_bCM; //  0.5420; % [ms-1] Ca2+ dissociation constant for calmodium
    // double k_fCM; //  227.7; % [mM-1 ms-1] Ca2+ associaton constant for calmodium
    double k_iso; //  0.1; % [1/min] Maximal AC activity
    double k_05_iso; //  3.34; % [nM] Half-maximal AC activation
    double n_iso; //  0.68; % Hill coefficient
    double K_ACI; //  0.016; % [min-1] Non-Ca2+ AC activity
    double K_AC; //  0.0735; % [min-1] Non-Ca2+ AC activity
    double K_Ca; //  0.000178; % [mM] Maximal Ca2+ AC activation
    double K_AC_Ca; //  0.000024; % [mM] Half-maximal Ca2+ AC activation
    double k_PDE; //  98500; % [mg protein/nmol/min] - Maximal PDE activity
    double k_PKA; //  9000; % [pmol/protein/min] - Maximal PKA activity
    double n_PKA; //  5; % Hill coefficient
    double k_PKA_cAMP; //  284.5; % [pmol/protein] - Half-maximal PKA activation
    double ATP_max; //  2.533; % [mM]
    double cAMPb; //  20; % [pmol/mg]
    double ADPm; //  0.276; % [mM]
    
    // %% PKA constants
    double PKI_tot; //  0.3; % [pmol/protin] - Total amount of PKA inhibitor
    double PKA_tot; //  1; % [pmol/protein] - Total amount of PKA
    
    // %% PLB phosphorylation constants
    double PP1; //  0.89; % [uM] PP1 concentration
    double k_PLBp; //  52.25; % [1/min] Maximal PLB phosphorylation
    double n_PLB; //  1; % Hill coefficient
    double k_PKA_PLB; //  1.651; % Half-maximal PLB phosphorylation
    double k_PP1; //  23.575; % [1/uM/min] - Maximal PP1 activity
    double k_pp1_PLB; //  0.06967; % Half maximal PP1 activity
    
    // %% I_f cAMP activity constants
    // %%
    double K_if; //  24.4; % [mV] - maximal I_f activation by cAMP
    double n_if; //  9.281; % Hill coefficient
    double K_05if; //  17.57; % [pmol/mg protein] - Half-maximal I_f activation by cAMP
    
    // %% Mithocondrial Ca2+ fluxes
    // %%
    double beta_Ca; //  0.1; % [] fraction of Ca2+ that binds to Ca2+ buffers in the mitochondria
    double P_Ca; //  100*25.6e-4; % [1/ms] uniporter Ca2+ permeability
    double phi_m; //  154.226; % [mV] mitochontrial membrane potential
    double alpha_m; //  0.2; % [] mitochondrial activity coefficients
    double alpha_e; //  0.341; % [] extracmitochondrial activity coefficients
    double Q_mo; //  2.4*25.6e-4; % [mM/ms] Na+-Ca2+ exchanger maximal velocity
    double K_Cam; //  0.003; % [mM] Na+ - Ca2+ exchanger Ca2+ affinity
    
    // %% I_Ks phosphorilation
    // %%
    double V_smin; //  257.1429; % [mV] Minimal Iks activation shift by PKA
    double V_smax; //  400; % [mV] Maximal Iks activation shift by PKA
    double g_kmax; //  25; % [] Maximal Iks activation by PKA
    double g_K05; //  0.5; % [] Half-maximal Iks activation by PKA
    double k_05k; //  0.4; % [] Half-maximal Iks activation shift by PKA
    double g_kmin; //  14.7541; % [] Minimal Iks activation by PKA
    
    // %% ATP-ADP
    // %%
    double kATP; //  61.42*100; % []
    double k_ATP05; //  6724; % []
    // double cAMPb; //  20; % [pmol/mg protein] - Baseline cAMP
    double n_ATP; //  3.36; % []
    double K_ATPmin; //  6034; % []
    double CATPi; // 2.6; % [mM] Total nucleotide concentrations
    double k_i_up; //  0.14; % [mM]
    double k_NaK_ATP; //  8*1e-3; % [mM]
    double k_NaK_ADP; //  0.1; % [mM]
    // double CATPi; //  2.6; % [mM]
    // %con.CATPm; //  1.5; % [mM]
    
    // %% Force
    // %%
    double SL_lo; //  0.8e-6; % [m] A constant coefﬁcient that describes the effect of the actin- and myosin-ﬁlament lengths on the single overlap length.
    double Nc; //  2e13; % [1/mm^2] The SAN cross-section area
    double FK0; //  350;  % [1/mM] The cross-bridge independent coefﬁcient of calcium afﬁnity
    double FK1; //  3e3; % [1/mM] The cooperativity coefﬁcient. Describes the dependence of calcium afﬁnity on the number of strong cross-bridges.
    double FN; //  3.5; % Hill coefficient
    double FK05; //  2.5e9; % [1/mm^3] Half-maximal cross-bridge Ca2+ affinity
    double Fkl; //  60; % [1/mM/ms] The rate constant of calcium binding to troponin low-afﬁnity sites
    double Ff; //  40e-3; % [1/ms] The cross-bridge turnover rate from the weak to the strong conformation
    double Fg0; //  30e-3; % [1/ms] The cross-bridge weakening rate at isometric regime
    double Fg1; //  4.4e6; % [1/m] The mechanical-feedback coefﬁcient. Describes the dependence of the XB weakening rate on the shortening velocity
    double Fxb; //  2e-9; % [mN] The unitary force per cross-bridge at isometric regime
    
    // %% ATP utilizers
    // %%
    double Max_ATP; //  12e-4; % [mM]
    double K_M_ATP; //  0.03; % [mM]
    double K_M_ADP; //  0.26; % [mM]
    
    // %% L-type phosphorilation constants
    // % ==> all constants currently in the function itself
};

void load_constants( cell_con * con );
//function con = load_constants()
void load_constants( cell_con * con ) {
    
    
    // %% Fixed ion concentration
    con->Cao = 2; // % [mM] extracellular Ca2+
    con->Ko = 5.4; // % [mM] extracellular K+
    con->Ki = 140; // % [mM] intracellular K+
    con->Nao = 140; // % [mM] extracellular Na+
    con->Nai = 10; // % [mM] intracellular Na+
    con->Mgi = 2.5; // % [mM] intracellular Mg+
    
    // %% Cell compartments
    con->C = 32; // % [pF] cell electric capacitance
    con->L_cell = 70; // % [um] cell length
    con->R_cell = 4; // % [um] cell radius
    con->L_sub = 0.02; // % [um] distrance between JSR and surface membrane (submembrane space)
    con->V_cell = pi * pow( con->R_cell, 2) * con->L_cell * 0.001; //10^-3; // % [pL] cell volume (3.5185838)
    con->V_sub = 2. * pi * con->L_sub * ( con->R_cell - con->L_sub / 2 ) * con->L_cell *0.001; //10^-3; // % [pL] submembrane cell volume (0.035097874) -> 10^-3 to make it pL
    con->V_jSR_part = 0.0012; // % [pL] part of cell volume occupied by junctional SR
    con->V_jSR = con->V_jSR_part * con->V_cell; // % [pL] volume of junctional SR (Ca2+ release store)
    con->V_i_part = 0.46; // % part of the cell volume occupied with myoplasm
    con->V_i = con->V_i_part * con->V_cell - con->V_sub; // % [pL] myoplasmic volume
    con->V_nSR_part = 0.0116; // % part of the cell volume occupied by network SR
    con->V_nSR = con->V_nSR_part*con->V_cell; // % [pL] volume of network SR (Ca2+ uptake store)
    con->V_myto = 0.6334; // % [pL] mitochondrial volume
    
    // %% The Nernst equation
    con->F = 96485; // % C/M is the Faraday constant
    con->T = 310.15; // % K is the absolute emperature for 37*C
    con->R = 8.3144; // % J/(M K) is the universal gaz constant
    con->E_T = 1000*(con->R*con->T/con->F); // % [mV] is the 'RT/F' factor
    
    // %% Electric potential
    con->E_Na = con->E_T*log(con->Nao/con->Nai); // % [mV] equilibrium potential for Na+
    con->E_K = con->E_T*log(con->Ko/con->Ki); // % [mV] equilibrium potential for K+
    con->E_Ks = con->E_T*log((con->Ko+0.12*con->Nao)/(con->Ki+0.12*con->Nai)); // % [mV] reversal potential of I_Ks
    con->E_CaL = 45; // % [mV] apparent reversal potential of I_CaL
    con->E_CaT = 45; // % [mV] apparent reversal potential of I_CaT
    con->E_st = 37.4; // % [mV] apparent reversal potential of I_st
    
    // %% Sarcolemmal Ion currents and their normalized conductances (g_x)
    con->g_CaL = 0.9; // %0.5720;% [nS/pF] normalized conductance for I_CaL channels
    con->g_CaT = 0.1832; // % [nS/pF] normalized conductance for I_CaT channels
    con->g_If = 0.455;// * 0.2;//0.455; // %0.07; // % [nS/pF] normalized conductance for I_f channels
    con->g_st = 0.00001; // % [nS/pF] normalized conductance for I_st channels
    con->g_Kr = 0.1; // % 0.08113973; // % [nS/pF] normalized conductance for I_Kr channels
    con->g_Ks = 0.013; // % [nS/pF] normalized conductance for I_Ks channels
    con->g_to = 0.252; // % [nS/pF] normalized conductance for I_to channels
    con->g_sus = 0.02; // % [nS/pF] normalized conductance for I_sus channels
    con->g_bCa = 0.0006; // % [nS/pF] normalized conductance for I_bCa channels
    con->g_bNa = 0.00486; // % [nS/pF] normalized conductance for I_bNa channels
    con->g_KACh_max = 2.2*0.14241818; // % [nS/pF] normalized conductance for I_KACh channels
    con->I_NaKmax = 2.88; // % [pA/pF] Maximal Na/K pump current conductance
    con->k_NCX = 225; // % [pA/pF] Maximal Na+/Ca2+ exchanger current conductance
    con->V_If12 = -64; // % -43.14; // % [mV] Half activation voltage for I_f current in the basal state
    
    // %% Modulation of sarcolemmal ion current by ions
    con->K_mfCa = 0.00035; // % [mM] dissociation constant of Ca2+ dependent I_CaL inactivation
    con->K_mKp = 1.4; // % [mM] half-maximal Ko for I_NaK
    con->K_mNap = 14; // % [mM] half-maximal Nai for I_NaK
    con->beta_fCa = 60; // % [mM-1.ms-1] Ca2+ association rate constant for I_CaL
    con->alpha_fCa = 0.021; // % [ms-1] Ca2+ dissociation rate constant for I_CaL
    
    // %% NCX function
    con->K_1ni = 395.3; // % [mM] Nai binding to the first site on NCX
    con->K_2ni = 2.289; // % [mM] Nai binding to second site on NCX
    con->K_3ni = 26.44; // % [mM] Nai binding to third site on NCX
    con->K_1no = 1628; // % [mM] Nao binding to first site on NCX
    con->K_2no = 561.4; // % [mM] Nao binding to second site on NCX
    con->K_3no = 4.663; // % [mM] Nao binding to third site on NCX
    con->K_ci = 0.0207; // % [mM] Cai binding on NCX transporter
    con->K_co = 3.663; // % [mM] Cao binding to NCX transporter
    con->K_cni = 26.44; // % [mM] Nai and Cai simultaneous binding to NCX
    con->Q_ci = 0.1369; // % Cai occlusion reaction of NCX
    con->Q_co = 0; // % Cao occlusion reaction of NCX
    con->Q_n = 0.4315; // % Na occlusion of NCX
    
    // %% Ca2+ diffusion
    con->tho_difCa = 0.04; // % [ms] Time constant of Ca diffusion from the submembrane to myoplasm
    con->tho_tr = 40; // % [ms] Time constant for Ca Ca transport from the network to junctional SR
    
    // %% SR Ca2+ ATPase function
    con->K_up = 0.6* 0.001; //10^-3; // % Half-maximal Cai for Ca uptake in the network SR [mM]
    con->P_up = 0.8*0.012; // % Rate constant for Ca uptake by the Ca pump in the network SR [mM/ms]
    con->k_m_up = 0.01;
    con->k_i_up_ADP = 5.1; // % [mM]
    
    // %% RyR function
    con->k_oCa_max = 10; // % [1/(mM^2*ms^1)]
    con->k_om = 0.06; // % [1/ms]
    con->k_iCa = 0.5; // % [1/(mM*ms)]
    con->k_im = 0.005; // % [1/ms]
    con->EC_50SR = 0.45; // % [mM]
    con->k_s = 400e3;
    con->MaxSR = 13;
    con->MinSR = 1;
    con->HSR = 3;
    // % RyR phosphorylation constants
    con->RyR_min = 0.0127; // % derived from PP1 activity
    con->RyR_max = 0.02;
    con->n_RyR = 9.773;
    con->k_05Ry = 0.7;
    
    // %% Ca2+ and Mg2+ buffering
    con->k_bCM = 0.542; // % [1/ms] Ca dissociation constant for calmodulin
    con->k_bCQ = 0.445; // % [1/ms] Ca dissociation constant for calsequestrin
    con->k_bTC = 0.446; // % [1/ms] Ca dissociation constant for the troponin-Ca site
    con->k_bTMC = 0.00751; // % [1/ms] Ca dissociation constant for the troponin-Mg site
    con->k_bTMM = 0.751; // % [1/ms] Mg dissociation constant for the troponin-Mg site
    con->k_fCM = 227.7; // % [1/mM*ms] Ca association constant for calmodulin
    con->k_fCQ = 0.534; // % [1/mM*ms] Ca association constant for calsequestrin
    con->k_fTC = 88.8; // % [mm/ms] Ca association constant for troponin
    con->k_fTMC = 227.7; // % [mM/ms] Ca association constant for the troponin-Mg site
    con->k_fTMM = 2.277; // % [mM/ms] Mg association constant for the troponin-Mg site
    con->TC_tot = 0.042; // % [mM] Total concentration of the troponin-Ca site
    // %con->TMC_tot = 0.062; // % [mM] Total concentration of the troponin-Mg site
    con->CQ_tot = 10; // % [mM] Total calsequestrin concentration
    con->CM_tot = 0.045; // % [mM] Total calmodium concentration
    
    // %% cAMP equations constants
    con->k_bCM = 0.5420; // % [ms-1] Ca2+ dissociation constant for calmodium
    con->k_fCM = 227.7; // % [mM-1 ms-1] Ca2+ associaton constant for calmodium
    con->k_iso = 0.1; // % [1/min] Maximal AC activity
    con->k_05_iso = 3.34; // % [nM] Half-maximal AC activation
    con->n_iso = 0.68; // % Hill coefficient
    con->K_ACI = 0.016; // % [min-1] Non-Ca2+ AC activity
    con->K_AC = 0.0735; // % [min-1] Non-Ca2+ AC activity
    con->K_Ca = 0.000178; // % [mM] Maximal Ca2+ AC activation
    con->K_AC_Ca = 0.000024; // % [mM] Half-maximal Ca2+ AC activation
    con->k_PDE = 98500; // % [mg protein/nmol/min] - Maximal PDE activity
    con->k_PKA = 9000; // % [pmol/protein/min] - Maximal PKA activity
    con->n_PKA = 5; // % Hill coefficient
    con->k_PKA_cAMP = 284.5; // % [pmol/protein] - Half-maximal PKA activation
    con->ATP_max = 2.533; // % [mM]
    con->cAMPb = 20; // % [pmol/mg]
    con->ADPm = 0.276; // % [mM]
    
    // %% PKA constants
    con->PKI_tot = 0.3; // % [pmol/protin] - Total amount of PKA inhibitor
    con->PKA_tot = 1; // % [pmol/protein] - Total amount of PKA
    
    // %% PLB phosphorylation constants
    con->PP1 = 0.89; // % [uM] PP1 concentration
    con->k_PLBp = 52.25; // % [1/min] Maximal PLB phosphorylation
    con->n_PLB = 1; // % Hill coefficient
    con->k_PKA_PLB = 1.651; // % Half-maximal PLB phosphorylation
    con->k_PP1 = 23.575; // % [1/uM/min] - Maximal PP1 activity
    con->k_pp1_PLB = 0.06967; // % Half maximal PP1 activity
    
    // %% I_f cAMP activity constants
    // %%
    con->K_if = 24.4; // % [mV] - maximal I_f activation by cAMP
    con->n_if = 9.281; // % Hill coefficient
    con->K_05if = 17.57; // % [pmol/mg protein] - Half-maximal I_f activation by cAMP
    
    // %% Mithocondrial Ca2+ fluxes
    // %%
    con->beta_Ca = 0.1; // % [] fraction of Ca2+ that binds to Ca2+ buffers in the mitochondria
    con->P_Ca = 100*25.6e-4; // % [1/ms] uniporter Ca2+ permeability
    con->phi_m = 154.226; // % [mV] mitochontrial membrane potential
    con->alpha_m = 0.2; // % [] mitochondrial activity coefficients
    con->alpha_e = 0.341; // % [] extracmitochondrial activity coefficients
    con->Q_mo = 2.4*25.6e-4; // % [mM/ms] Na+-Ca2+ exchanger maximal velocity
    con->K_Cam = 0.003; // % [mM] Na+ - Ca2+ exchanger Ca2+ affinity
    
    // %% I_Ks phosphorilation
    // %%
    con->V_smin = 257.1429; // % [mV] Minimal Iks activation shift by PKA
    con->V_smax = 400; // % [mV] Maximal Iks activation shift by PKA
    con->g_kmax = 25; // % [] Maximal Iks activation by PKA
    con->g_K05 = 0.5; // % [] Half-maximal Iks activation by PKA
    con->k_05k = 0.4; // % [] Half-maximal Iks activation shift by PKA
    con->g_kmin = 14.7541; // % [] Minimal Iks activation by PKA
    
    // %% ATP-ADP
    // %%
    con->kATP = 61.42*100; // % []
    con->k_ATP05 = 6724; // % []
    con->cAMPb = 20; // % [pmol/mg protein] - Baseline cAMP
    con->n_ATP = 3.36; // % []
    con->K_ATPmin = 6034; // % []
    con->CATPi= 2.6; // % [mM] Total nucleotide concentrations
    con->k_i_up = 0.14; // % [mM]
    con->k_NaK_ATP = 8*1e-3; // % [mM]
    con->k_NaK_ADP = 0.1; // % [mM]
    con->CATPi = 2.6; // % [mM]
    // %con->CATPm = 1.5; // % [mM]
    
    // %% Force
    // %%
    con->SL_lo = 0.8e-6; // % [m] A constant coefﬁcient that describes the effect of the actin- and myosin-ﬁlament lengths on the single overlap length.
    con->Nc = 2e13; // % [1/mm^2] The SAN cross-section area
    con->FK0 = 350;  // % [1/mM] The cross-bridge independent coefﬁcient of calcium afﬁnity
    con->FK1 = 3e3; // % [1/mM] The cooperativity coefﬁcient. Describes the dependence of calcium afﬁnity on the number of strong cross-bridges.
    con->FN = 3.5; // % Hill coefficient
    con->FK05 = 2.5e9; // % [1/mm^3] Half-maximal cross-bridge Ca2+ affinity
    con->Fkl = 60; // % [1/mM/ms] The rate constant of calcium binding to troponin low-afﬁnity sites
    con->Ff = 40e-3; // % [1/ms] The cross-bridge turnover rate from the weak to the strong conformation
    con->Fg0 = 30e-3; // % [1/ms] The cross-bridge weakening rate at isometric regime
    con->Fg1 = 4.4e6; // % [1/m] The mechanical-feedback coefﬁcient. Describes the dependence of the XB weakening rate on the shortening velocity
    con->Fxb = 2e-9; // % [mN] The unitary force per cross-bridge at isometric regime
    
    // %% ATP utilizers
    // %%
    con->Max_ATP = 12e-4; // % [mM]
    con->K_M_ATP = 0.03; // % [mM]
    con->K_M_ADP = 0.26; // % [mM]
    
    // %% L-type phosphorilation constants
    // % ==> all constants currently in the function itself
    
    // end
}
