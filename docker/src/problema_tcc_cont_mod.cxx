//////////////////////////////////////////////////////////////////////////
//////////////////        problema_tcc.cxx         ///////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           Roger  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//// Title:         Inserção em órbita com tempo mínimo      /////////////
//////// Last modified: 04 February 2021                   ///////////////
//////// Reference:     My TCC                            ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Roger M. Sehnem, 2021          ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"

//////////////////////////////////////////////////////////////////////////
/////////  Declare an auxiliary structure to hold local constants  ///////
//////////////////////////////////////////////////////////////////////////

struct Continua {
    adouble aero;
    adouble r0;
    adouble r0_init;
};

typedef struct Continua Continua_;

struct Otimo {
    adouble T1;
    adouble T2;
    adouble T3;
    adouble k1;
    adouble k2;
    adouble k3;
    adouble m0;
    adouble m12;
    adouble m23;
    adouble s12;
    adouble s23;
    adouble r0;
    adouble sigma0;
    adouble u0;
    adouble v0;
    adouble rf;
    adouble uf;
    adouble vf;
    bool aero;
    Continua_ Continua;
};

typedef struct Otimo Otimo_;

struct Foguete {
  adouble Sref;
  adouble i;
  adouble mue;
  Otimo_ Otimo;
};

typedef struct Foguete Foguete_;

//////////////////////////////////////////////////////////////////////////
///////////////////  Declare auxiliary functions  ////////////////////////
//////////////////////////////////////////////////////////////////////////

void aerodin(const adouble &time, adouble  *states, adouble  *controls,
             const int iphase, const Foguete_ &Foguete, adouble &FA, 
             adouble &FN);

adouble Vsom(const adouble &hkm);

adouble rho(const adouble &hm);

adouble CD0(const adouble &Mach, const int ie);

adouble CNa(const adouble &Mach, const int ie);

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return (0.0);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    return  1.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
/////// This is the "orbitaode.m" function in the Matlab program /////////
//////////////////////////////////////////////////////////////////////////
/* A transformação tau=2t/(tf-t0)-(tf+t0)/(tf-t0) é aplicada em cada fase
automáticamente, nada precisa ser feito para a conversão. Desta forma o 
tempo de cada fase é t_i^(j) <= t <= t_f^(j)*/

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

    Foguete_ &Foguete = *( (Foguete_ *) workspace->problem->user_data );

    // extrai estados
    adouble r = states[0], sigma = states[1], u = states[2];
    adouble v = states[3], m = states[4];

    //extrai controles
    adouble sin_beta = controls[0], cos_beta = controls[1];

    //calcula forças aerodinamicas
    adouble FA = 0, FN = 0;

    // resolver com aerodinamica completa?
    if (Foguete.Otimo.aero)
    {
        aerodin(time,states,controls,iphase,Foguete,FA,FN); // chama a função aerodinamica para calcular as forças aerodinamicas
        // usa apenas percentual de FA e FN dado um passo
        FA *= Foguete.Otimo.Continua.aero;
        FN *= Foguete.Otimo.Continua.aero;
    }

    const adouble mi = Foguete.mue;

    adouble g=-mi/pow(r,2);

    adouble T = 0, ki = 0;

    if (iphase == 1)
    {
        T=Foguete.Otimo.T1;
        ki=Foguete.Otimo.k1;
    }
    else if (iphase == 2)
    {
        T=Foguete.Otimo.T2;
        ki=Foguete.Otimo.k2;
    }
    else if (iphase == 3)
    {
        T=Foguete.Otimo.T3;
        ki=Foguete.Otimo.k3;
    }

    adouble rDot, sigmaDot, uDot, vDot, mDot;

    // Dinamica do foguete cos_beta
    rDot = u;
    sigmaDot = v/r;
    uDot = (T*sin_beta-FA*sin_beta+FN*cos_beta)/m+g+pow(v,2)/r;
    vDot = (T*cos_beta-FA*cos_beta-FN*sin_beta)/m-(u*v)/r;
    mDot = -ki;

    derivatives[0] = rDot; derivatives[1] = sigmaDot; derivatives[2] = uDot;
    derivatives[3] = vDot; derivatives[4] = mDot;

    path[ 0 ] = dot( controls, controls, 2);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
    Foguete_ &Foguete = *( (Foguete_ *) workspace->problem->user_data );

    adouble r0 = Foguete.Otimo.r0,              rf = Foguete.Otimo.rf;
    adouble sigma0 = Foguete.Otimo.sigma0,      uf = Foguete.Otimo.uf;
    adouble u0 = Foguete.Otimo.u0,              vf = Foguete.Otimo.vf;
    adouble v0 = Foguete.Otimo.v0;
    adouble m0 = Foguete.Otimo.m0;
   
    if (iphase == 1)
    {   
        adouble ri = initial_states[ 0 ];
        adouble sigmai = initial_states[ 1 ];
        adouble ui = initial_states[ 2 ];
        adouble vi = initial_states[ 3 ];
        adouble mi = initial_states[ 4 ];
        // continuação em r0
        e[ 0 ] = ri-(r0*Foguete.Otimo.Continua.r0+Foguete.Otimo.Continua.r0_init*(1-Foguete.Otimo.Continua.r0)); // r0_init p/ cont.r0=0 e r0 p/ cont.r0=1
        e[ 1 ] = sigmai-sigma0;
        e[ 2 ] = ui-u0;
        e[ 3 ] = vi-v0;
        e[ 4 ] = mi-m0;
    }
    else if (iphase == 3)
    {   
        adouble rfi = final_states[ 0 ];
        adouble ufi = final_states[ 2 ];
        adouble vfi = final_states[ 3 ];
        e[ 0 ] = rfi-rf;
        e[ 1 ] = ufi-uf;
        e[ 2 ] = vfi-vf;
        //for (int i=0; i != 3; ++i)
        //   std::cout <<"e["<<i<<"]" << e[i] << std::endl;
    }

    //int dummy;
    //for (int i=0; i != 5; ++i)
    //   std::cout <<"initial states:" << initial_states[i] << std::endl;
    //std::cin >> dummy;
}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    // funcao orbitabc do matlab (exceto pelos estados finais e iniciais)
    
    adouble xf_p1[5], xi_p2[5], xf_p2[5], xi_p3[5], tf_p1, ti_p2, tf_p2, ti_p3;

    get_final_states(xf_p1, xad, 1, workspace); // a fase é iphase
    get_final_states(xf_p2, xad, 2, workspace);
    get_initial_states(xi_p2, xad, 2, workspace);
    get_initial_states(xi_p3, xad, 3, workspace);
    tf_p1 = get_final_time(xad, 1, workspace);
    tf_p2 = get_final_time(xad, 2, workspace);
    ti_p2 = get_initial_time(xad, 2, workspace);
    ti_p3 = get_initial_time(xad, 3, workspace);

    // pega estrutura Foguete
    Foguete_ &Foguete = *( (Foguete_ *) workspace->problem->user_data );

    // extrai dados das massas e saltos na massa

    // massa inicial
    adouble m0=Foguete.Otimo.m0;
    // massa gasta até a troca do primeiro-segundo estágio
    adouble m12=Foguete.Otimo.m12;
    // segundo-terceiro
    adouble m23=Foguete.Otimo.m23;
    // Salto na massa primeiro-segundo
    adouble s12=Foguete.Otimo.s12;
    // segundo-terceiro
    adouble s23=Foguete.Otimo.s23;

    // avalia links psis nos estados
    /*for (int j = 0; j != 4; j++) // estados físicos continuos, exceto a massa
    {
        linkages[j] = xf_p1[j]-xi_p2[j];
    }
    //linkages[4] = xf_p1[4]-(m0-m12); // massa em que o salto na massa ocorre
    linkages[4] = xf_p1[4]-s12-xi_p2[4];

    //int dummy;
    //std::cout << xf_p1[4] << " <-> " << xi_p2[4] << linkages[5] << std::endl;
    //std::cin >> dummy;

    for (int j = 0; j != 4; j++) // estados físicos continuos, exceto a massa
    {
        linkages[j+5] = xf_p2[j]-xi_p3[j];
    }
    //linkages[10] = xf_p2[4]-(m0-m12-s12-m23); // massa em que o salto na massa ocorre
    linkages[9] = xf_p2[4]-s23-xi_p3[4];

    // continuidade e links no tempo
    linkages[10] = tf_p1-ti_p2;
    linkages[11] = tf_p2-ti_p3;*/

    int index=0;

    auto_link(linkages, &index, xad, 1, 2, workspace );
    linkages[index-2] -= s12;
    auto_link(linkages, &index, xad, 2, 3, workspace );
    linkages[index-2] -= s23;
    // linkages[12] = xf_p1[4]-(m0-m12);
    // linkages[13] = xf_p2[4]-(m0-m12-s12-m23);

    //int dummy;
    //for (int i=0; i != 14; ++i)
    //   std::cout << linkages[i] << std::endl;
    //std::cin >> dummy;
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Problem do TCC";
    problem.outfilename         = "out_file.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 3;
    problem.nlinkages           = 12;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    for (int j=1; j != 4; j++)
    { 
        problem.phases(j).nstates   		= 5;
        problem.phases(j).ncontrols 		= 2;
        problem.phases(j).npath     		= 1;
    }

    problem.phases(1).nevents   		= 5;
    problem.phases(2).nevents   		= 0;
    problem.phases(3).nevents   		= 3;

    problem.phases(1).nodes             = (RowVectorXi(2) <<15,18).finished(); // numeros de nós em sequencia
    problem.phases(2).nodes             = (RowVectorXi(2) <<15,18).finished();
    problem.phases(3).nodes             = (RowVectorXi(2) <<15,18).finished();

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

////////////////////////////////////////////////////////////////////////////
//////////  Declare an instance of Foguete and Otimo structure /////////////
////////////////////////////////////////////////////////////////////////////

    
    Foguete_ Foguete; Otimo_ Otimo; Continua_ Continua;
    
    problem.user_data = (void*) &Foguete;

////////////////////////////////////////////////////////////////////////////
///////////////////  Initialize CONSTANTS and //////////////////////////////
///////////////////  declare local variables  //////////////////////////////
////////////////////////////////////////////////////////////////////////////
    
    Continua.r0 = 1;
    Continua.aero = 1;

    Otimo.T1 = 1169300;
    Otimo.T2 = 314710;
    Otimo.T3 = 206940;
    Otimo.k1 = 458.5397;
    Otimo.k2 = 115.8710;
    Otimo.k3 = 5; // 20 deve ser 76.7586 (errado tbm no matlab)
    Otimo.m12 = 28888;
    Otimo.m23 = 7184;
    Otimo.s12 = 5922;
    Otimo.s23 = 1664;
    Otimo.Continua = Continua;

    Foguete.Sref = 0.7964;
    Foguete.i = 0.0402; //0.1745;
    Foguete.mue = 3.986012000000000e+14;
    Foguete.Otimo=Otimo;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double rMin = 6378145, sigmaMin = 0, uMin = 0, vMin = 0, mMin = 0;
    double rMax = 10e6, sigmaMax = 10000, uMax = 100000, vMax = 100000, mMax = 5*50591;

    for (int iphase = 1; iphase != 4; iphase++)
    {
        problem.phases(iphase).bounds.lower.states << rMin, sigmaMin, uMin, vMin, mMin;
        problem.phases(iphase).bounds.upper.states << rMax, sigmaMax, uMax, vMax, mMax;
        problem.phases(iphase).bounds.lower.controls << -1, -1;
        problem.phases(iphase).bounds.upper.controls << 1, 1;
        problem.phases(iphase).bounds.lower.path << 1;
        problem.phases(iphase).bounds.upper.path << 1;
    }

    problem.phases(1).bounds.lower.StartTime    = 0;
    problem.phases(1).bounds.upper.StartTime    = 0;
    problem.phases(1).bounds.lower.EndTime    = 28888/458.5397;
    problem.phases(1).bounds.upper.EndTime    = 528888/458.5397;

    problem.phases(2).bounds.lower.StartTime    = 28888/458.5397;
    problem.phases(2).bounds.upper.StartTime    = 28888/458.5397;
    problem.phases(2).bounds.lower.EndTime    = 28888/458.5397+7184/115.8710;
    problem.phases(2).bounds.upper.EndTime    = 28888/458.5397+7184/115.8710;

    problem.phases(3).bounds.lower.StartTime    = 28888/458.5397+7184/115.8710;
    problem.phases(3).bounds.upper.StartTime    = 28888/458.5397+7184/115.8710;
    problem.phases(3).bounds.lower.EndTime    = 28888/458.5397+7184/115.8710;
    problem.phases(3).bounds.upper.EndTime    = 500;

    double r0_init = 7e6, r0 = 6378146;
    double sigma0 = 0, u0 = 0, v0 = 463.8318, m0 = 50591;
    //double r0 = 7e6, sigma0 = 0, u0 = 0, v0 = 7.5e3, m0 = 5*50591;
    double rf = 7071000, uf = 0, vf = 7.5081e3;
    // double rf = 7.5e6, uf = 0, vf = 7.5081e3;
    //double rf = 9e6, uf = 0, vf = sqrt(3.986012000000000e+14/9e6);

    // Define condições de contorno
    Foguete.Otimo.r0 = r0; Foguete.Otimo.sigma0 = sigma0; Foguete.Otimo.u0 = u0;
    Foguete.Otimo.v0 = v0; Foguete.Otimo.m0 = m0; Foguete.Otimo.Continua.r0_init=r0_init;

    Foguete.Otimo.rf = rf; Foguete.Otimo.uf = uf; Foguete.Otimo.vf = vf;

    // events deve ser igual a zero
    for (int i = 0; i != 5; i++)
    {
        problem.phases(1).bounds.lower.events(i)=0.0;
        problem.phases(1).bounds.upper.events(i)=0.0;
    }
    for  (int i = 0; i != 3;  i++)
    {
        problem.phases(3).bounds.lower.events(i)=0.0;
        problem.phases(3).bounds.upper.events(i)=0.0;
    }


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    for (int iphase = 1; iphase != 4; iphase++)
    {
        problem.phases(iphase).guess.states = zeros(5,15);

        problem.phases(iphase).guess.states.row(1) = linspace( sigma0, sigma0, 15);
        problem.phases(iphase).guess.states.row(2) = linspace( u0, uf, 15);

        problem.phases(iphase).guess.controls = zeros(2,15);
    }

    problem.phases(1).guess.states.row(0) = linspace( r0_init , r0_init+ (rf-r0_init)*1/3, 15);
    problem.phases(2).guess.states.row(0) = linspace(  r0_init+ (rf-r0_init)*1/3 ,  r0_init+ (rf-r0_init)*2/3, 15);
    problem.phases(3).guess.states.row(0) = linspace(  r0_init+ (rf-r0_init)*2/3 , rf, 15);

    problem.phases(1).guess.states.row(3) = linspace( v0, vf/3, 15);
    problem.phases(2).guess.states.row(3) = linspace( vf/3, vf*2/3, 15);
    problem.phases(3).guess.states.row(3) = linspace( vf*2/3, vf, 15);

    problem.phases(1).guess.states.row(4) = linspace( m0, 50591-28888, 15);
    problem.phases(2).guess.states.row(4) = linspace( 50591-28888-5922, 50591-28888-5922-7184, 15);
    problem.phases(3).guess.states.row(4) = linspace( 50591-28888-5922-7184, 1000, 15);

    double  t0 = 0, t1=63, t2=63+62, t3=350;

    problem.phases(1).guess.time = linspace(t0,t1, 15);
    problem.phases(2).guess.time = linspace(t1,t2, 15);
    problem.phases(3).guess.time = linspace(t2,t3, 15);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    // algorithm.derivatives                 = "numerical";
    algorithm.mesh_refinement             = "automatic";
    // algorithm.mr_max_iterations           = 15;
    algorithm.collocation_method          = "Chebyshev";
    // algorithm.collocation_method          = "Hermite-Simpson";
    // algorithm.collocation_method          = "trapezoidal";
    // algorithm.defect_scaling = "jacobian-based";
    //algorithm.mr_min_extrapolation_points = 6;
    // algorithm.ipopt_linear_solver         = "ma27";
    // algorithm.diff_matrix                 = "central-differences";
    // algorithm.hessian                     = "exact";
    // algorithm.nlp_tolerance               = 1e-9;
    // algorithm.ode_tolerance               = 1e-9;
    algorithm.nsteps_error_integration    = 30;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    // resolve sem aero primeiro
    Foguete.Otimo.aero = true;

    psopt(solution, problem, algorithm);

    // std::cout << std::endl <<"******************************************************************************"<< std::endl <<
    //              "******************************************************************************"<< std::endl <<
    //              "Solução sem aerodinâmica encontrada" << std::endl<< 
    //              "******************************************************************************"<< std::endl <<
    //              "******************************************************************************"<<std::endl<<std::endl;

    // std::cout << '\a';


    // change algorithm options
    // algorithm.collocation_method          = "Legendre";
    // algorithm.collocation_method          = "Hermite-Simpson";
    // algorithm.nlp_tolerance               = 1e-9;
    // algorithm.ode_tolerance               = 1e-9;
    // algorithm.mesh_refinement             = "automatic";
    // algorithm.mr_max_iterations           = 15;

    // passos para de continuação
    double steps_aero = 2, steps_r0 = 10; 

     // resolve com aero
    Foguete.Otimo.aero = true;

    // loops da continuação, vai de 0 a 1
//     for (double i = 0; i != steps_aero+1; i++)
//     {
//         // usa solução como chute inicial
//         for (int iphase = 1; iphase != 4; iphase++)
//         {   
//             problem.phases(iphase).guess.states = solution.get_states_in_phase(iphase);
//             problem.phases(iphase).guess.controls = solution.get_controls_in_phase(iphase);
//             problem.phases(iphase).guess.time = solution.get_time_in_phase(iphase);
//             problem.phases(iphase).guess.parameters = solution.get_parameters_in_phase(iphase);
//             // problem.phases(iphase).nodes = (RowVectorXi(1) << solution.nodes[iphase-1].size()).finished(); // usa o número de nós anterior
//             problem.phases(iphase).nodes = (RowVectorXi(1) << 18).finished();
//         }

//         Foguete.Otimo.Continua.aero = i/(steps_aero); // vai de 0 a 1 em steps passos
//         psopt(solution, problem, algorithm);
//     }

//     std::cout << std::endl <<"******************************************************************************"<< std::endl <<
//                  "******************************************************************************"<< std::endl <<
//                  "Solução com aerodinâmica encontrada" << std::endl<< 
//                  "******************************************************************************"<< std::endl <<
//                  "******************************************************************************"<<std::endl<<std::endl;
//     std::cout << '\a';

//     for (double i = 0; i != steps_r0+1; i++)
//     {
//         // usa solução como chute inicial
//         for (int iphase = 1; iphase != 4; iphase++)
//         {   
//             problem.phases(iphase).guess.states = solution.get_states_in_phase(iphase);
//             problem.phases(iphase).guess.controls = solution.get_controls_in_phase(iphase);
//             problem.phases(iphase).guess.time = solution.get_time_in_phase(iphase);
//             problem.phases(iphase).guess.parameters = solution.get_parameters_in_phase(iphase);
//             // problem.phases(iphase).nodes = (RowVectorXi(1) << solution.nodes[iphase-1].size()).finished();
//             problem.phases(iphase).nodes = (RowVectorXi(1) << 18).finished();
//         }

//         Foguete.Otimo.Continua.r0 = i/(steps_r0); // vai de 0 a 1 em steps passos
//         psopt(solution, problem, algorithm);
        
//     }

//    std::cout << '\a';
   

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x_ph1, x_ph2, x_ph3, u_ph1, u_ph2, u_ph3, lambda_ph1, lambda_ph2, lambda_ph3;
    MatrixXd t_ph1, t_ph2, t_ph3, H_ph1, H_ph2, H_ph3;
    
    x_ph1 = solution.get_states_in_phase(1);
    x_ph2 = solution.get_states_in_phase(2);
    x_ph3 = solution.get_states_in_phase(3);
    
    u_ph1 = solution.get_controls_in_phase(1);
    u_ph2 = solution.get_controls_in_phase(2);
    u_ph3 = solution.get_controls_in_phase(3); 
    
    t_ph1 = solution.get_time_in_phase(1);
    t_ph2 = solution.get_time_in_phase(2);
    t_ph3 = solution.get_time_in_phase(3);    

    lambda_ph1 = solution.get_dual_costates_in_phase(1);
    lambda_ph2 = solution.get_dual_costates_in_phase(2);
    lambda_ph3 = solution.get_dual_costates_in_phase(3);

    H_ph1 = solution.get_dual_hamiltonian_in_phase(1);
    H_ph2 = solution.get_dual_hamiltonian_in_phase(2);
    H_ph3 = solution.get_dual_hamiltonian_in_phase(3);

    x.resize(5, x_ph1.cols()+ x_ph2.cols()+ x_ph3.cols());
    u.resize(2, u_ph1.cols()+ u_ph2.cols()+ u_ph3.cols());
    t.resize(1, t_ph1.cols()+ t_ph2.cols()+ t_ph3.cols());
    lambda.resize(5, lambda_ph1.cols()+ lambda_ph2.cols()+ lambda_ph3.cols());
    H.resize(1, H_ph1.cols()+ H_ph2.cols()+ H_ph3.cols());

    x << x_ph1, x_ph2, x_ph3; 
    u << u_ph1, u_ph2, u_ph3;
    t << t_ph1, t_ph2, t_ph3;
    lambda << lambda_ph1, lambda_ph2, lambda_ph3;
    H  << H_ph1, H_ph2, H_ph3;

    MatrixXd beta = u.row(0);

    for (int i = 0; i =! u.row(0).size(); i++)
        beta(i) = atan2(u(0,i),u(1,i));


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    Save(x, "x.dat");
    Save(u,"u.dat");
    Save(t,"t.dat");
    Save(lambda,"lambda.dat");
    Save(H,"H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////


    // plot(t,x.row(0),problem.name+": state x", "time (s)", "x","x");

    // bool plotar = false;

    // std::cout << std::endl << "Plotar o resto das Figuras? " << std::endl << "1 p/ S 0 p/ N" << std::endl;
    // std::cin >> plotar;
    // if (plotar)
    // {
    //     plot(t,x.row(1),problem.name+": state sigma", "time (s)", "sigma","sigma");
    //     plot(t,x.row(2),problem.name+": state u", "time (s)", "u","u");
    //     plot(t,x.row(3),problem.name+": state v", "time (s)", "v","v");
    //     plot(t,x.row(4),problem.name+": state m", "time (s)", "m","m");

    //     plot(t,beta,problem.name+": controls","time (s)", "controls","beta");

    //     plot(t,lambda.row(0),problem.name+": coestate lx", "time (s)", "lx","lx");
    //     plot(t,lambda.row(1),problem.name+": coestate lsigma", "time (s)", "lsigma","lsigma");
    //     plot(t,lambda.row(2),problem.name+": coestate lu", "time (s)", "lu","lu");
    //     plot(t,lambda.row(3),problem.name+": coestate lv", "time (s)", "lv","lv");
    //     plot(t,lambda.row(4),problem.name+": coestate lm", "time (s)", "lm","lm");
        
    //     plot(t,H,problem.name+": Hamiltoniana", "time (s)", "H[]","H");
    // }

    // plot(t,x.row(0),problem.name+": state x", "time (s)", "x","x", "pdf", "estado_x.pdf");
    // plot(t,x.row(1),problem.name+": state sigma", "time (s)", "sigma","sigma", "pdf", "estado_y.pdf");
    // plot(t,x.row(2),problem.name+": state u", "time (s)", "u","u", "pdf", "estado_u.pdf");
    // plot(t,x.row(3),problem.name+": state v", "time (s)", "v","v", "pdf", "estado_v.pdf");
    // plot(t,x.row(4),problem.name+": state m", "time (s)", "m","m", "pdf", "estado_m.pdf");

    // plot(t,beta,problem.name+": controls","time (s)", "controls", "beta", "pdf", "beta.pdf");

    // plot(t,lambda.row(0),problem.name+": coestate lx", "time (s)", "lx","lx", "pdf", "coestado_lx.pdf");
    // plot(t,lambda.row(1),problem.name+": coestate lsigma", "time (s)", "lsigma","lsigma", "pdf", "coestado_ls.pdf");
    // plot(t,lambda.row(2),problem.name+": coestate lu", "time (s)", "lu","lu", "pdf", "coestado_lu.pdf");
    // plot(t,lambda.row(3),problem.name+": coestate lv", "time (s)", "lv","lv", "pdf", "coestado_lv.pdf");
    // plot(t,lambda.row(4),problem.name+": coestate lm", "time (s)", "lm","lm", "pdf", "coestado_lm.pdf");
    
    // plot(t,H,problem.name+": Hamiltoniana", "time (s)", "H[]","H","pdf","hamil.pdf");

    // plot(t,x,problem.name+": states", "time (s)", "states","x sigma u v m",
    //                         "pdf", "brymr_states.pdf");

    //plot(t,u,problem.name+": controls","time (s)", "controls", "beta",
    //                         "pdf", "brymr_states.pdf");
}

////////////////////////////////////////////////////////////////////////////
////////////////// Define auxiliary functions //////////////////////////////
////////////////////////////////////////////////////////////////////////////

void aerodin(const adouble &time, adouble  *states, adouble  *controls,
             const int iphase, const Foguete_ &Foguete, adouble &FA, 
             adouble &FN)
{
    // extrai controles
    adouble beta = atan2(controls[0],controls[1]);
    // extrai estados
    adouble r = states[0];
    adouble sigma = states[1];
    adouble u = states[2];
    adouble v = states[3];
    adouble m = states[4];

    adouble i = Foguete.i;
    const adouble WE = 7.292115856e-5;
    adouble aux = (v-r*WE*cos(i));
    adouble gama = atan2(u,aux); // funciona  do mesmo modo que o atan2 do matlab
    adouble alfa=beta-gama;
    const adouble RE=6378145; 
    adouble h=r-RE;
    if (h < 0) // garantindo que a altitude é positiva para n ferrar os modelos
        h=0;
    adouble VR=sqrt((pow(u,2)+pow((v-r*WE*cos(i)),2)));

    adouble Mach=VR/Vsom(h);

    // parte do FAero

    const adouble S = Foguete.Sref;
    adouble rhoe = rho(h);
    adouble q = .5*rhoe*pow(VR,2);

    // calcula FA
    adouble  Cx = CD0(Mach,iphase);
    FA = q*S*Cx;

    // calcula FN
    adouble CNx = CNa(Mach,iphase);
    FN = q*S*CNx*alfa;
}

adouble Vsom(const adouble &hm)
{
    // Função para a velocidade do som em função da altitude (Modelo do Paulo
    // Sérgio, Pg. 158.)
    //
    // Entradas: h [m] (altitude)
    const adouble A0=340;
    const adouble A1=-7.59813873;
    const adouble A2=4.41739218e-1;
    const adouble A3=-1.12391220e-2;
    const adouble A4=1.70917283e-4;
    const adouble A5=-1.59888891e-6;
    const adouble A6=6.51961005e-9;

    adouble h=hm/1e3; // transforma m p/ km
    adouble vs = A0+A1*h+A2*pow(h,2)+A3*pow(h,3)+A4*pow(h,4)+
                A5*pow(h,5)+A6*pow(h,6);
    return vs;
}

adouble rho(const adouble &hm)
{
    // Função que calcula a densidade do ar como função da altitude (Paulo
    // Sérgio. Pg 156
    //
    // Entradas: h [m] (altitude)

    const adouble rho0=1.225;
    const adouble k1=1.02280550;
    const adouble k2=1.21226930e-1;
    const adouble k3=-3.48643241e-2;
    const adouble k4=3.50991865e-3;
    const adouble k5=-8.33000533e-5;
    const adouble k6=1.15219733e-6;
    const adouble k7=2.204504;
    const adouble k8=3.12550644e-3;
    const adouble k9=5.82501032e-4;
    const adouble k10=7.44266792e-6;
    const adouble k11=3.88398609e-7;

    adouble  h=hm/1e3; // m p/ km
    adouble rh = 0;

    // if (h <= 65)
    // {
        adouble phi = k1*exp(-(k3*h+k4*pow(h,2)+k5*pow(h,3)+k6*pow(h,4)));
        rh = rho0*exp(-(k1+k2*h-phi));
    // }    

    return rh;
}

adouble CD0(const adouble &Mach, const int ie)
{   
    const adouble E = 2.718281828459046;
    adouble CD = 0; 
    if (ie==1)
    {
        CD=2.1469002185913078 -0.12985835302355783*(-(1/(1 +pow(E,32.6071161645013
            -6.375163653461174 * Mach))) +1/(1 + pow(E,-4.271647777532277
            +12.498619352150442 * Mach))) -3.3492312615775073*(-(1/(1 +pow(E,-3.4365836650245463
            +1.690203772082327 * Mach))) +1/(1 + pow(E,-12.030758568673956 +13.04861114710195 * Mach)));
    }
    else if (ie==2)
    {
        CD=0.9464291281418414 -0.33570516081004914*(-(1/(1 +pow(E,12.264303954215212
         -6.081521881854877 * Mach))) +1/(1 + pow(E,6.944630189974732 -1.8850103540724115 * Mach))) 
         -0.5061172842895649*(1/(1 + pow(E,-5.09511067121503 +0.8218810971026423 * Mach))
         -1/(1 + pow(E,-11.188493920637118 +5.201042464440416 * Mach)));
    }

    return CD;
    
}

adouble CNa(const adouble &Mach, const int ie)
{   
    const adouble E = 2.718281828459046;
    adouble CD = 0; 
    if (ie==1)
    {
        CD=21.408256872962873 +3.046370389879771/(1 + pow(E,-18.82015342286602 +5.411680324379552 * Mach))
            - 3.2166758427052713*(1/(1 + pow(E,0.7160708433853873 +3.488000598084749 * Mach)) -1/(1 + pow(E,-8.808790572851724
            +4.070534916988994 * Mach))) +(-160.278400329585*pow(Mach,3) +831.0490437163307*pow(Mach,4) -1385.4497280418225*pow(Mach,5) 
            +295.77505936694376*pow(Mach,6) +1134.4545098395427*pow(Mach,7) -715.6296555451819*pow(Mach,8))/(1 
            + pow(E,-15.173236414749201 +16.833623217058175 * Mach));
    }
    else if (ie==2)
    {
        CD=3.335031207941182 +1.1338365259835768/(1 + pow(E,30.966226472324983 -30.47085810182766 * Mach))
            - 0.19683687178648548*(-(1/(1 +pow(E,14.88955509269116 +23.361955543192877 * Mach))) 
            +1/(1 + pow(E,-17.560612011448253 +30.087414446944024 * Mach))) +(-0.19548929306102572*pow(Mach,3)
             +18.043217225374175*pow(Mach,4) -65.87808649373738*pow(Mach,5) +75.93541706025756*pow(Mach,6) 
             -34.390552397751755*pow(Mach,7) +5.492980102821696*pow(Mach,8))/(1 + pow(E,-2.5901792078381107 
             +3.7748175249034936 * Mach));
    }

    return CD;
    
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
