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

    //////////////////////////////////////////////////////////////////////////
    /////////////////  Declare MatrixXd objects to store results //////////////
    //////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd lambda, H;

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

    problem.phases(1).bounds.lower.states << x1_lower, x2_lower, x3_lower, x4_lower;

    problem.phases(1).bounds.upper.states << x1_upper, x2_upper, x3_upper, x4_upper;

    problem.phases(1).bounds.lower.controls << u_lower;
    problem.phases(1).bounds.upper.controls << u_upper;

    // all events must be zero
    for (int i = 0; i != 8; i++)
    {
        problem.phases(1).bounds.lower.events(i)=0.0;
        problem.phases(1).bounds.upper.events(i)=0.0;
    }

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = .0001;
    problem.phases(1).bounds.upper.EndTime      = 300.0;

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
    algorithm.nlp_tolerance               = 1.e-4;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.mesh_refinement             = "automatic";
    // algorithm.collocation_method          = "trapezoidal";
    algorithm.collocation_method          = "Chebyshev";
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-6;
    algorithm.nsteps_error_integration    = 30;

    ////////////////////////////////////////////////////////////////////////////
    //////////////////  Define Initial Conditions  /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Carefull, the total number of solutions will be in the power of 4.
    const int n_per_dim = 2;
    MatrixXd x1is = zeros(1,n_per_dim);
    MatrixXd x2is = zeros(1,n_per_dim);
    MatrixXd x3is = zeros(1,n_per_dim);
    MatrixXd x4is = zeros(1,n_per_dim);

    // Define the minimum and maximum values for the initial conditions
    double x1_lower_ci = -1, x2_lower_ci = -10, x3_lower_ci = -3.1415, x4_lower_ci = - 5;
    double x1_upper_ci = 1 , x2_upper_ci = 10 , x3_upper_ci = 3.1415 , x4_upper_ci = 5;
    x1is = linspace(x1_lower_ci, x1_upper_ci, n_per_dim);
    x2is = linspace(x2_lower_ci, x2_upper_ci, n_per_dim);
    x3is = linspace(x3_lower_ci, x3_upper_ci, n_per_dim);
    x4is = linspace(x4_lower_ci, x4_upper_ci, n_per_dim);
    
    // Now compute the solutions for every possible combination of the initial
    // states defined above.
    for (int i = 0; i != n_per_dim; i++)
    {   
        for (int j = 0; j != n_per_dim; j++)
        {   
            for (int k = 0; k != n_per_dim; k++)
            {   
                for (int l = 0; l != n_per_dim; l++)
                {   
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

                    double x1_init = x1is.coeff(0,i);
                    double x2_init = x2is.coeff(0,j);
                    double x3_init = x3is.coeff(0,k);
                    double x4_init = x4is.coeff(0,l);

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

                    MatrixXd x_guess    =  zeros(nstates,nnodes);

                    x_guess.row(0) = linspace(x1_init, 0, nnodes);
                    x_guess.row(1) = linspace(x2_init, 0, nnodes);
                    x_guess.row(2) = linspace(x3_init, 0, nnodes);
                    x_guess.row(3) = linspace(x4_init, 0, nnodes);
                


                    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
                    problem.phases(1).guess.states         = x_guess;
                    problem.phases(1).guess.time           = linspace(0.0,1.0,nnodes);
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
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
