//////////////////////////////////////////////////////////////////////////
////////////////        Inveted Pendulum.cxx         /////////////////////
////////////////           Roger M. Sehnem           /////////////////////
////////////////          Kenedy M. Portella         /////////////////////
////// Title: Optimal control to zero origin of inverted pendulum  ///////


#include "psopt.h"

//////////////////////////////////////////////////////////////////////////
/////////  Declare an auxiliary structure to hold local constants  ///////
//////////////////////////////////////////////////////////////////////////

// Holds the constants of the inverted pendulum.
struct Constants_pendulum {
    adouble rod_length;
    adouble rod_inertia;
    adouble rod_mass;
    adouble cart_mass;
    adouble translation_friction_coefficient;
    adouble rotation_friction_coefficient;
};

typedef struct Constants_pendulum Constants_pendulum_;

// Holds the initial conditions of the inverted pendulum.
struct initial_conditions {
    adouble x1;
    adouble x2;
    adouble x3;
    adouble x4;
};

typedef struct initial_conditions initial_conditions_;

// Holds the constants for the optimal control formulation.
struct Optimal_settings {
    double r;
    double q;
};

typedef struct Optimal_settings Optimal_settings_;

// Holds the constants for the optimal control formulation.
struct User_Data {
    Optimal_settings_ Optimal_settings;
    initial_conditions_ initial_conditions;
    Constants_pendulum_ Constants_pendulum;
};

typedef struct User_Data User_Data_;


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
    // Declare structure that holds parameters for the optimal control
    // formulation. Pure C++ whitchcraft, basically  magic that I once did
    // and don't remember what it means.
    User_Data_ &User_Data = *( (User_Data_ *) workspace->problem->user_data);

    // The constants
    double r = User_Data.Optimal_settings.r;
    double q = User_Data.Optimal_settings.q;

    // This is the cost function of the optimal control problem
    return states[0]*q*states[0] + \
           states[1]*q*states[1] + \
           states[2]*q*states[2] + \
           states[3]*q*states[3] + \
           controls[0]*r*controls[0] + \
           q*time;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
///////////// The dynamics of the pendulum is defined here ///////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   // Structure that holds parameters of the pendulum
   User_Data_ &User_Data = *( (User_Data_ *) workspace->problem->user_data);

   // Extract parameters from the structure
   adouble l = User_Data.Constants_pendulum.rod_length;
   adouble I = User_Data.Constants_pendulum.rod_inertia;
   adouble mb = User_Data.Constants_pendulum.rod_mass;
   adouble mc = User_Data.Constants_pendulum.cart_mass;
   adouble at = User_Data.Constants_pendulum.translation_friction_coefficient;
   adouble ar = User_Data.Constants_pendulum.rotation_friction_coefficient;

   // Defining helper variables
   adouble g = 9.8; // The acceleration of gravity
   adouble C1 = l*mb;
   adouble C2 = I+pow(l,2)*mb;
   adouble C3 = mb+mc; // Those are just helper constants

   // Unpacking states
   adouble x1 = states[ 0 ];
   adouble x2 = states[ 1 ];
   adouble x3 = states[ 2 ];
   adouble x4 = states[ 3 ];

   // Unpacking control
   adouble u = controls[ 0 ];
   adouble x1_dot, x2_dot, x3_dot, x4_dot;

   x1_dot = x2;
   x2_dot = (g*pow(C1,2)*cos(x3)*sin(x3)+C2*(u-at*x2)-ar*C1*cos(x3)*x4-C1*C2*sin(x3)*pow(x4,2))/(C2*C3-pow(C1,2)*pow(cos(x3),2));
   x3_dot = x4;
   x4_dot = (g*C1*C3*sin(x3)+C1*cos(x3)*(u-at*x2)-ar*C3*x4-pow(C1,2)*cos(x3)*sin(x3)*pow(x4,2))/(C2*C3-pow(C1,2)*pow(cos(x3),2));

   derivatives[ 0 ] = x1_dot;
   derivatives[ 1 ] = x2_dot;
   derivatives[ 2 ] = x3_dot;
   derivatives[ 3 ] = x4_dot;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   // Structure that holds initial conditions of the pendulum
   User_Data_ &User_Data = *( (User_Data_ *) workspace->problem->user_data);
   adouble x1_init = User_Data.initial_conditions.x1;
   adouble x2_init = User_Data.initial_conditions.x2;
   adouble x3_init = User_Data.initial_conditions.x3;
   adouble x4_init = User_Data.initial_conditions.x4;

   e[ 0 ] = x1_init - initial_states[0];
   e[ 1 ] = x2_init - initial_states[1];
   e[ 2 ] = x3_init - initial_states[2];
   e[ 3 ] = x4_init - initial_states[3];
   e[ 4 ] = final_states[0];
   e[ 5 ] = final_states[1];
   e[ 6 ] = final_states[2];
   e[ 7 ] = final_states[3];

}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
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

    problem.name        		= "Inverted pendulum problem";
    problem.outfilename                 = "inv_pend_prob.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 4;
    problem.phases(1).ncontrols 		= 1;
    problem.phases(1).nevents   		= 8;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes         << 50;

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare MatrixXd objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

////////////////////////////////////////////////////////////////////////////
//////////  Declare an instance structs/////////////
////////////////////////////////////////////////////////////////////////////

    initial_conditions_ initial_conditions;
    Constants_pendulum_ Constants_pendulum;
    Optimal_settings_ Optimal_settings;
    User_Data_ User_Data;
    
    problem.user_data = (void*) &User_Data;

////////////////////////////////////////////////////////////////////////////
///////////////////  Initialize CONSTANTS and //////////////////////////////
///////////////////  declare local variables  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Constants_pendulum.rod_length = 0.3;
    Constants_pendulum.rod_inertia = 2;
    Constants_pendulum.rod_mass = 1;
    Constants_pendulum.cart_mass = 3;
    Constants_pendulum.translation_friction_coefficient = 0.2;
    Constants_pendulum.rotation_friction_coefficient = 0.2;

    double x1_init = 0.1;
    double x2_init = 0.1;
    double x3_init = 0.1;
    double x4_init = 0.02;

    initial_conditions.x1 = x1_init;
    initial_conditions.x2 = x2_init;
    initial_conditions.x3 = x3_init;
    initial_conditions.x4 = x4_init;

    Optimal_settings.r = 1;
    Optimal_settings.q = 1;

    User_Data.Optimal_settings = Optimal_settings;
    User_Data.initial_conditions = initial_conditions;
    User_Data.Constants_pendulum = Constants_pendulum;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    double x1_lower = -10.0;
    double x2_lower = -100.0;
    double x3_lower = -6.28;
    double x4_lower = -6.28;

    double x1_upper = 10.0;
    double x2_upper = 100.0;
    double x3_upper = 6.28;
    double x4_upper = 6.28;   

    double u_lower = -5000;
    double u_upper = 5000;

    problem.phases(1).bounds.lower.states(0) = x1_lower;
    problem.phases(1).bounds.lower.states(1) = x2_lower;
    problem.phases(1).bounds.lower.states(2) = x3_lower;
    problem.phases(1).bounds.lower.states(4) = x4_lower;

    problem.phases(1).bounds.upper.states(0) = x1_upper;
    problem.phases(1).bounds.upper.states(1) = x2_upper;
    problem.phases(1).bounds.upper.states(2) = x3_upper;
    problem.phases(1).bounds.upper.states(4) = x4_upper;

    problem.phases(1).bounds.lower.controls(0) = u_lower;
    problem.phases(1).bounds.upper.controls(0) = u_upper;

    problem.phases(1).bounds.lower.events(0) = x1_init;
    problem.phases(1).bounds.lower.events(1) = x2_init;
    problem.phases(1).bounds.lower.events(2) = x3_init;
    problem.phases(1).bounds.lower.events(3) = x4_init;

    problem.phases(1).bounds.upper.events(0) = x1_init;
    problem.phases(1).bounds.upper.events(1) = x2_init;
    problem.phases(1).bounds.upper.events(2) = x3_init;
    problem.phases(1).bounds.upper.events(3) = x4_init;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 30.0;
    problem.phases(1).bounds.upper.EndTime      = 30.0;

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

    int nnodes    			            = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd x_guess    =  zeros(nstates,nnodes);

    x_guess.row(0)  = x1_init*ones(1,nnodes);
    x_guess.row(1)  = x2_init*ones(1,nnodes);
    x_guess.row(2)  = x3_init*ones(1,nnodes);
    x_guess.row(3)  = x4_init*ones(1,nnodes);


    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0,15.0,nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
//    algorithm.mesh_refinement             = "automatic";
    algorithm.collocation_method          = "trapezoidal";
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-6;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


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

//    plot(t,x,problem.name+": states", "time (s)", "states","x y v");

//     plot(t,u,problem.name+": controls","time (s)", "controls", "u_1 u_2");

//     plot(t,x,problem.name+": states", "time (s)", "states","x y v",
//                              "pdf", "brymr_states.pdf");

//     plot(t,u,problem.name+": controls","time (s)", "controls", "u_1 u_2",
//                              "pdf", "brymr_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
