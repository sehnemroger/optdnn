//////////////////////////////////////////////////////////////////////////
////////////////        Inveted Pendulum.cxx         /////////////////////
////////////////           Roger M. Sehnem           /////////////////////
////////////////          Kenedy M. Portella         /////////////////////
////// Title: Optimal control to zero origin of inverted pendulum  ///////


#include "psopt.h"

// Include this to convert floats to strings properly
#include <iostream>
#include <sstream>

// See https://stackoverflow.com/questions/2125880/convert-float-to-stdstring-in-c
template < typename T > std::string to_str (const T& t)
{
  std::ostringstream os;
  os << t;
  return os.str ();
}

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
    adouble r;
    adouble q;
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
   return (pow(final_states[0],2)+pow(final_states[1],2)+pow(final_states[2],2)+pow(final_states[3],2));
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
    adouble r = User_Data.Optimal_settings.r;
    adouble q = User_Data.Optimal_settings.q;

    // This is the cost function of the optimal control problem
    return states[0]*q*states[0] + states[1]*q*states[1] + states[2]*q*states[2] + states[3]*q*states[3] + controls[0]*r*controls[0] + q*time;
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

    if ( iphase == 1) {
        e[ 0 ] = x1_init - initial_states[0];
        e[ 1 ] = x2_init - initial_states[1];
        e[ 2 ] = x3_init - initial_states[2];
        e[ 3 ] = x4_init - initial_states[3];
    }
    else if ( iphase == 3){
        e[ 0 ] = final_states[0];
        e[ 1 ] = final_states[1];
        e[ 2 ] = final_states[2];
        e[ 3 ] = final_states[3];
    }
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
    MSdata msdata;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Inverted pendulum problem";
    problem.outfilename                 = "inv_pend_prob.txt";

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    msdata.nsegments = 3;
    msdata.nstates  = 4;
    msdata.ncontrols = 1;
    msdata.nparameters = 0;
    msdata.npath = 0;
    msdata.ninitial_events = 4;
    msdata.nfinal_events = 4;
    msdata.nodes << 20; // nodes per segment

    multi_segment_setup(problem, algorithm, msdata);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////     Declare an instance structs      //////////////////
    ////////////////////////////////////////////////////////////////////////////

    initial_conditions_ initial_conditions;
    Constants_pendulum_ Constants_pendulum;
    Optimal_settings_ Optimal_settings;
    User_Data_ User_Data;
    
    problem.user_data = (void*) &User_Data;


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    double x1_lower = -15.0;
    double x2_lower = -50.0;
    double x3_lower = -6.28;
    double x4_lower = -6.28;

    double x1_upper = 15.0;
    double x2_upper = 50.0;
    double x3_upper = 6.28;
    double x4_upper = 6.28;   

    double u_lower = -1000;
    double u_upper = 1000;

    problem.phases(1).bounds.lower.states << x1_lower, x2_lower, x3_lower, x4_lower;
    problem.phases(1).bounds.upper.states << x1_upper, x2_upper, x3_upper, x4_upper;
    problem.phases(1).bounds.lower.controls << u_lower;
    problem.phases(1).bounds.upper.controls << u_upper;

    // all events must be zero
    for (int i = 0; i != 4; i++)
    {
        problem.phases(1).bounds.lower.events(i)=0.0;
        problem.phases(1).bounds.upper.events(i)=0.0;
    }

    problem.phases(3).bounds.lower.states << x1_lower, x2_lower, x3_lower, x4_lower;
    problem.phases(3).bounds.upper.states << x1_upper, x2_upper, x3_upper, x4_upper;
    problem.phases(3).bounds.lower.controls << u_lower;
    problem.phases(3).bounds.upper.controls << u_upper;

    // all events must be zero
    for (int i = 0; i != 4; i++)
    {
        problem.phases(3).bounds.lower.events(i)=0.0;
        problem.phases(3).bounds.upper.events(i)=0.0;
    }

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(3).bounds.lower.EndTime      = .4;
    problem.phases(3).bounds.upper.EndTime      = 300.0;

    // problem.bounds.lower.times.resize(1,4);
    // problem.bounds.upper.times.resize(1,4);

    // problem.bounds.lower.times << 0.0, 0.1, 0.2, 0.3;
    // problem.bounds.upper.times << 0.0, 20.0, 21.0, 22.0;

    auto_phase_bounds(problem);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.mesh_refinement             = "automatic";
    // algorithm.collocation_method          = "Hermite-Simpson";
    // algorithm.collocation_method          = "Chebyshev";
    // algorithm.defect_scaling = "jacobian-based";
    // algorithm.nsteps_error_integration    = 30;
    // algorithm.mr_max_increment_factor     = 0.8;
    // algorithm.hessian                     = "exact";
    // algorithm.nlp_tolerance               = 1e-9;
    algorithm.ode_tolerance               = 1.e-5;
    // algorithm.mr_min_extrapolation_points = 4;
    // algorithm.mr_initial_increment        = 20;
    // algorithm.mr_max_iterations           = 20;
    // algorithm.jac_sparsity_ratio          = 1;
    // algorithm.hess_sparsity_ratio         = 1;
   
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
    double x3_init = 3.1415;
    double x4_init = 0.1;

    initial_conditions.x1 = x1_init;
    initial_conditions.x2 = x2_init;
    initial_conditions.x3 = x3_init;
    initial_conditions.x4 = x4_init;

    Optimal_settings.r = 1.0;
    Optimal_settings.q = 1.0;

    User_Data.Optimal_settings = Optimal_settings;
    User_Data.initial_conditions = initial_conditions;
    User_Data.Constants_pendulum = Constants_pendulum;
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int nnodes    			            = problem.phases(1).nodes(0);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    MatrixXd state_guess    =  zeros(nstates,nnodes);

    state_guess.row(0) = linspace(x1_init, 0.0, nnodes);
    state_guess.row(1) = linspace(x2_init, 0.0, nnodes);
    state_guess.row(2) = linspace(x3_init, 0.0, nnodes);
    state_guess.row(3) = linspace(x4_init, 0.0, nnodes);

    MatrixXd control_guess = zeros(ncontrols,nnodes);
    MatrixXd time_guess = linspace(0.0,5.0,nnodes);
    MatrixXd param_guess;

    auto_phase_guess(problem, control_guess, state_guess, param_guess, time_guess);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   /////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t, x_ph1, u_ph1, t_ph1, x_ph2, u_ph2, t_ph2;
    MatrixXd x_ph3, u_ph3, t_ph3;


    x_ph1      = solution.get_states_in_phase(1);
    u_ph1      = solution.get_controls_in_phase(1);
    t_ph1      = solution.get_time_in_phase(1);

    x_ph2      = solution.get_states_in_phase(2);
    u_ph2      = solution.get_controls_in_phase(2);
    t_ph2      = solution.get_time_in_phase(2);

    x_ph3      = solution.get_states_in_phase(3);
    u_ph3      = solution.get_controls_in_phase(3);
    t_ph3      = solution.get_time_in_phase(3);

    x.resize(4, length(t_ph1) + length(t_ph2) + length(t_ph3) );
    u.resize(1, length(t_ph1) + length(t_ph2) + length(t_ph3) );
    t.resize(1, length(t_ph1) + length(t_ph2) + length(t_ph3) );

    x << x_ph1, x_ph2, x_ph3;
    u << u_ph1, u_ph2, u_ph3;
    t << t_ph1, t_ph2, t_ph3;
    // lambda = solution.get_dual_costates_in_phase(1);
    // H      = solution.get_dual_hamiltonian_in_phase(1);

    ////////////////////////////////////////////////////////////////////////////
    /////////  Save solution with the name of the initial conditions    ////////
    ////////////////////////////////////////////////////////////////////////////

    // Define the output matrix to save to file,
    // First collumn is time, second to fifth are the states and the sixth 
    // is the control. 
    MatrixXd out;

    // The data comes from the solver in the line form (as opposed to the
    // collumn form)

    out.resize(6, t.cols()); // get in the correct shape
    out << t, x, u; // pass in the data

    // Create name
    std::string x_name = "q_" + to_str(1) + "r_" + to_str(1);
    x_name += "x1i_" + to_str(x1_init);
    x_name += "x2i_" + to_str(x2_init);
    x_name += "x3i_" + to_str(x3_init);
    x_name += "x4i_" + to_str(x4_init) + ".dat";

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Save solution data to files if desired ////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Save(out.transpose(), x_name.c_str());
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
