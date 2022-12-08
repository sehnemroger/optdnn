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
template <typename T>
std::string to_str(const T &t)
{
    std::ostringstream os;
    os << t;
    return os.str();
}

//////////////////////////////////////////////////////////////////////////
/////////  Declare an auxiliary structure to hold local constants  ///////
//////////////////////////////////////////////////////////////////////////

// Holds the constants of the inverted pendulum.
struct Constants_pendulum
{
    adouble rod_length;
    adouble rod_inertia;
    adouble rod_mass;
    adouble cart_mass;
    adouble translation_friction_coefficient;
    adouble rotation_friction_coefficient;
};

typedef struct Constants_pendulum Constants_pendulum_;

// Holds the initial conditions of the inverted pendulum.
struct initial_conditions
{
    adouble x1;
    adouble x2;
    adouble x3;
    adouble x4;
};

typedef struct initial_conditions initial_conditions_;

// Holds the constants for the optimal control formulation.
struct Optimal_settings
{
    adouble r;
    adouble q;
    int phases;
};

typedef struct Optimal_settings Optimal_settings_;

// Holds the constants for the optimal control formulation.
struct User_Data
{
    Optimal_settings_ Optimal_settings;
    initial_conditions_ initial_conditions;
    Constants_pendulum_ Constants_pendulum;
};

typedef struct User_Data User_Data_;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble *initial_states, adouble *final_states,
                      adouble *parameters, adouble &t0, adouble &tf,
                      adouble *xad, int iphase, Workspace *workspace)
{
    return 0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble *states, adouble *controls, adouble *parameters,
                       adouble &time, adouble *xad, int iphase, Workspace *workspace)
{
    // Declare structure that holds parameters for the optimal control
    // formulation. Pure C++ whitchcraft, basically  magic that I once did
    // and don't remember what it means.
    User_Data_ &User_Data = *((User_Data_ *)workspace->problem->user_data);

    // The constants
    adouble r = User_Data.Optimal_settings.r;
    adouble q = User_Data.Optimal_settings.q;

    // This is the cost function of the optimal control problem
    return states[0] * q * states[0] + states[1] * q * states[1] + states[2] * q * states[2] + states[3] * q * states[3] + controls[0] * r * controls[0] + q * time;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
///////////// The dynamics of the pendulum is defined here ///////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble *derivatives, adouble *path, adouble *states,
         adouble *controls, adouble *parameters, adouble &time,
         adouble *xad, int iphase, Workspace *workspace)
{
    // Structure that holds parameters of the pendulum
    User_Data_ &User_Data = *((User_Data_ *)workspace->problem->user_data);

    // Extract parameters from the structure
    adouble l = User_Data.Constants_pendulum.rod_length;
    adouble I = User_Data.Constants_pendulum.rod_inertia;
    adouble mb = User_Data.Constants_pendulum.rod_mass;
    adouble mc = User_Data.Constants_pendulum.cart_mass;
    adouble at = User_Data.Constants_pendulum.translation_friction_coefficient;
    adouble ar = User_Data.Constants_pendulum.rotation_friction_coefficient;

    // Defining helper variables
    adouble g = 9.8; // The acceleration of gravity
    adouble C1 = l * mb;
    adouble C2 = I + pow(l, 2) * mb;
    adouble C3 = mb + mc; // Those are just helper constants

    // Unpacking states
    adouble x1 = states[0];
    adouble x2 = states[1];
    adouble x3 = states[2];
    adouble x4 = states[3];

    // Unpacking control
    adouble u = controls[0];
    adouble x1_dot, x2_dot, x3_dot, x4_dot;

    x1_dot = x2;
    x2_dot = (g * pow(C1, 2) * cos(x3) * sin(x3) + C2 * (u - at * x2) - ar * C1 * cos(x3) * x4 - C1 * C2 * sin(x3) * pow(x4, 2)) / (C2 * C3 - pow(C1, 2) * pow(cos(x3), 2));
    x3_dot = x4;
    x4_dot = (g * C1 * C3 * sin(x3) + C1 * cos(x3) * (u - at * x2) - ar * C3 * x4 - pow(C1, 2) * cos(x3) * sin(x3) * pow(x4, 2)) / (C2 * C3 - pow(C1, 2) * pow(cos(x3), 2));

    derivatives[0] = x1_dot;
    derivatives[1] = x2_dot;
    derivatives[2] = x3_dot;
    derivatives[3] = x4_dot;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble *e, adouble *initial_states, adouble *final_states,
            adouble *parameters, adouble &t0, adouble &tf, adouble *xad,
            int iphase, Workspace *workspace)
{
    // Structure that holds initial conditions of the pendulum
    User_Data_ &User_Data = *((User_Data_ *)workspace->problem->user_data);

    adouble x1_init = User_Data.initial_conditions.x1;
    adouble x2_init = User_Data.initial_conditions.x2;
    adouble x3_init = User_Data.initial_conditions.x3;
    adouble x4_init = User_Data.initial_conditions.x4;
    int phases = User_Data.Optimal_settings.phases;

    if (iphase == 1)
    {
        e[0] = x1_init - initial_states[0];
        e[1] = x2_init - initial_states[1];
        e[2] = x3_init - initial_states[2];
        e[3] = x4_init - initial_states[3];
    }
    else if (iphase == phases)
    {
        e[0] = final_states[0];
        e[1] = final_states[1];
        e[2] = final_states[2];
        e[3] = final_states[3];
        // forces the state derivatives to be zero at the final time
        adouble state_derivative_final_1;
        adouble state_derivative_final_2;
        adouble state_derivative_final_3;
        adouble state_derivative_final_4;
        get_state_derivative(&state_derivative_final_1, 0, iphase, tf, xad, workspace);
        get_state_derivative(&state_derivative_final_2, 1, iphase, tf, xad, workspace);
        get_state_derivative(&state_derivative_final_3, 2, iphase, tf, xad, workspace);
        get_state_derivative(&state_derivative_final_4, 3, iphase, tf, xad, workspace);
        e[4] = state_derivative_final_1;
        e[5] = state_derivative_final_2;
        e[6] = state_derivative_final_3;
        e[7] = state_derivative_final_4;
        // forces the final control to be zero
        // adouble controls;
        // get_final_controls(&controls, xad, iphase, workspace);
        // e[4] = controls;
    }
    // if ( iphase == 1 ) {
    //     e[ 0 ] = x1_init - initial_states[0];
    //     e[ 1 ] = x2_init - initial_states[1];
    //     e[ 2 ] = x3_init - initial_states[2];
    //     e[ 3 ] = x4_init - initial_states[3];
    //     e[ 4 ] = final_states[0] - x1_init*2/3;
    //     e[ 5 ] = final_states[1] - x2_init*2/3;
    //     e[ 6 ] = final_states[2] - x3_init*2/3;
    //     e[ 7 ] = final_states[3] - x4_init*2/3;
    // }
    // else if ( iphase == 2) {
    //     e[ 0 ] = initial_states[0] - x1_init*2/3;
    //     e[ 1 ] = initial_states[1] - x2_init*2/3;
    //     e[ 2 ] = initial_states[2] - x3_init*2/3;
    //     e[ 3 ] = initial_states[3] - x4_init*2/3;
    //     e[ 4 ] = final_states[0] -  x1_init*1/3;
    //     e[ 5 ] = final_states[1] -  x2_init*1/3;
    //     e[ 6 ] = final_states[2] -  x3_init*1/3;
    //     e[ 7 ] = final_states[3] -  x4_init*1/3;
    // }
    // else if ( iphase == 3) {
    //     e[ 0 ] = initial_states[0] - x1_init*1/3;
    //     e[ 1 ] = initial_states[1] - x2_init*1/3;
    //     e[ 2 ] = initial_states[2] - x3_init*1/3;
    //     e[ 3 ] = initial_states[3] - x4_init*1/3;
    //     e[ 4 ] = final_states[0];
    //     e[ 5 ] = final_states[1];
    //     e[ 6 ] = final_states[2];
    //     e[ 7 ] = final_states[3];
    // }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages(adouble *linkages, adouble *xad, Workspace *workspace)
{

    // Structure that holds initial conditions of the pendulum
    User_Data_ &User_Data = *((User_Data_ *)workspace->problem->user_data);

    int index = 0;
    adouble tf_p, tf; // final time of the phases
    int phases = User_Data.Optimal_settings.phases;

    for (int iphase = 1; iphase != phases; iphase++)
    {
        auto_link_2(linkages, &index, xad, iphase, iphase + 1, workspace);
    }

    // auto_link_2(linkages, &index, xad, 1, 2, workspace);
    // auto_link_2(linkages, &index, xad, 2, 3, workspace);
    // auto_link_2(linkages, &index, xad, 3, 4, workspace);

    tf = get_final_time(xad, phases, workspace); // final time

    for (int iphase = 1; iphase != phases; iphase++)
    {
        tf_p = get_final_time(xad, iphase, workspace);
        linkages[index] = tf_p - tf * iphase / phases;
        index++;
    }

    // linkages[index] = tf_p1 - tf_p4 / 4; // makes the final time of the first phase 1/4 of the final time
    // index++;
    // linkages[index] = tf_p2 - tf_p4 * 2 / 4; // makes the final time of the second phase 2/4 of the final time
    // index++;
    // linkages[index] = tf_p3 - tf_p4 * 3 / 4; // makes the final time 3/4 of the final time
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Alg algorithm;
    Sol solution;
    Prob problem;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.name = "Inverted pendulum problem";
    problem.outfilename = "inv_pend_prob.txt";

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    int phases = 6;
    int nstates = 4;
    int ncontrols = 1;

    problem.nphases = phases;
    // problem.nlinkages = 21;
    problem.nlinkages = (nstates + ncontrols + 1) * (phases - 1) + phases - 1;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    int nodes = 15;

    for (int iphases = 1; iphases != phases + 1; iphases++)
    {
        problem.phases(iphases).nstates = nstates;
        problem.phases(iphases).ncontrols = ncontrols;
        problem.phases(iphases).npath = 0;
        // problem.phases(iphases).nodes             = (RowVectorXi(2) << 15, 18).finished(); // numeros de nós em sequencia
        // problem.phases(iphases).nodes = (RowVectorXi(3) << 15, 20, 25).finished(); // numeros de nós em sequencia
        problem.phases(iphases).nodes << nodes;
        problem.phases(iphases).nevents = 0;
    }

    problem.phases(1).nevents = 4;
    problem.phases(phases).nevents = 8;

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////     Declare an instance structs      //////////////////
    ////////////////////////////////////////////////////////////////////////////

    initial_conditions_ initial_conditions;
    Constants_pendulum_ Constants_pendulum;
    Optimal_settings_ Optimal_settings;
    User_Data_ User_Data;

    problem.user_data = (void *)&User_Data;
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

    double x1_init = -1.0;
    double x2_init = 1.0;
    double x3_init = 3.1415;
    double x4_init = -1.0;

    initial_conditions.x1 = x1_init;
    initial_conditions.x2 = x2_init;
    initial_conditions.x3 = x3_init;
    initial_conditions.x4 = x4_init;

    Optimal_settings.r = 1.0;
    Optimal_settings.q = 1.0;
    Optimal_settings.phases = phases;

    User_Data.Optimal_settings = Optimal_settings;
    User_Data.initial_conditions = initial_conditions;
    User_Data.Constants_pendulum = Constants_pendulum;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    double x1_lower = -20.0;
    double x2_lower = -50.0;
    double x3_lower = -6.28;
    double x4_lower = -6.28;

    double x1_upper = 20.0;
    double x2_upper = 50.0;
    double x3_upper = 6.28;
    double x4_upper = 6.28;

    double u_lower = -1000;
    double u_upper = 1000;

    for (int iphase = 1; iphase != phases + 1; iphase++)
    {
        problem.phases(iphase).bounds.lower.states << x1_lower, x2_lower, x3_lower, x4_lower;
        problem.phases(iphase).bounds.upper.states << x1_upper, x2_upper, x3_upper, x4_upper;
        problem.phases(iphase).bounds.lower.controls << u_lower;
        problem.phases(iphase).bounds.upper.controls << u_upper;
    }

    // all events must be zero
    for (int i = 0; i != 4; i++)
    {
        problem.phases(1).bounds.lower.events(i) = 0.0;
        problem.phases(1).bounds.upper.events(i) = 0.0;
    }
    // all events must be zero
    for (int i = 0; i != 8; i++)
    {
        problem.phases(phases).bounds.lower.events(i) = 0.0;
        problem.phases(phases).bounds.upper.events(i) = 0.0;
    }

    double start_time = 0.0, end_time = 30.0;

    problem.phases(1).bounds.lower.StartTime = start_time;
    problem.phases(1).bounds.upper.StartTime = start_time;
    problem.phases(1).bounds.lower.EndTime = start_time;
    problem.phases(1).bounds.upper.EndTime = end_time;

    for (int iphase = 2; iphase != phases; iphase++)
    {
        problem.phases(iphase).bounds.lower.StartTime = start_time;
        problem.phases(iphase).bounds.upper.StartTime = end_time;
        problem.phases(iphase).bounds.lower.EndTime = start_time;
        problem.phases(iphase).bounds.upper.EndTime = end_time;
    }

    problem.phases(phases).bounds.lower.StartTime = start_time;
    problem.phases(phases).bounds.upper.StartTime = end_time;
    problem.phases(phases).bounds.lower.EndTime = start_time;
    problem.phases(phases).bounds.upper.EndTime = end_time;

    // problem.bounds.lower.times.resize(1,4);
    // problem.bounds.upper.times.resize(1,4);

    // problem.bounds.lower.times << 0.0, 0.1, 0.2, 0.3;
    // problem.bounds.upper.times << 0.0, 20.0, 21.0, 22.0;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae;
    problem.events = &events;
    problem.linkages = &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max = 2000;
    algorithm.nlp_method = "IPOPT";
    algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic";
    // algorithm.mesh_refinement = "automatic";
    // algorithm.collocation_method          = "Hermite-Simpson";
    // algorithm.collocation_method          = "Chebyshev";
    // algorithm.defect_scaling = "jacobian-based";
    // algorithm.nsteps_error_integration = 30;
    // algorithm.mr_max_increment_factor = 0.8;
    // algorithm.hessian                     = "exact";
    algorithm.nlp_tolerance = 1e-6;
    algorithm.ode_tolerance = 1.e-6;
    // algorithm.mr_min_extrapolation_points = 4;
    // algorithm.mr_initial_increment        = 20;
    // algorithm.mr_max_iterations = 20;
    // algorithm.jac_sparsity_ratio          = 1;
    // algorithm.hess_sparsity_ratio         = 1;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int nnodes = problem.phases(1).nodes(0);
    // int ncontrols = problem.phases(1).ncontrols;
    // int nstates = problem.phases(1).nstates;

    for (int iphase = 1; iphase != phases + 1; iphase++)
    {
        problem.phases(iphase).guess.states = zeros(nstates, nnodes);
        problem.phases(iphase).guess.controls = zeros(ncontrols, nnodes);

        float factor_init = (phases - iphase + 1) / phases;
        float factor_final = (phases - iphase) / phases;
        problem.phases(iphase).guess.states.row(0) = linspace(x1_init * factor_init, x1_init * factor_final, nnodes);
        problem.phases(iphase).guess.states.row(1) = linspace(x2_init * factor_init, x2_init * factor_final, nnodes);
        problem.phases(iphase).guess.states.row(2) = linspace(x3_init * factor_init, x3_init * factor_final, nnodes);
        problem.phases(iphase).guess.states.row(3) = linspace(x4_init * factor_init, x4_init * factor_final, nnodes);
    }

    // problem.phases(1).guess.states.row(0) = linspace(x1_init, x1_init * 3 / 4, nnodes);
    // problem.phases(1).guess.states.row(1) = linspace(x2_init, x2_init * 3 / 4, nnodes);
    // problem.phases(1).guess.states.row(2) = linspace(x3_init, x3_init * 3 / 4, nnodes);
    // problem.phases(1).guess.states.row(3) = linspace(x4_init, x4_init * 3 / 4, nnodes);

    // problem.phases(2).guess.states.row(0) = linspace(x1_init * 3 / 4, x1_init * 2 / 4, nnodes);
    // problem.phases(2).guess.states.row(1) = linspace(x2_init * 3 / 4, x2_init * 2 / 4, nnodes);
    // problem.phases(2).guess.states.row(2) = linspace(x3_init * 3 / 4, x3_init * 2 / 4, nnodes);
    // problem.phases(2).guess.states.row(3) = linspace(x4_init * 3 / 4, x4_init * 2 / 4, nnodes);

    // problem.phases(3).guess.states.row(0) = linspace(x1_init * 2 / 4, x1_init * 1 / 4, nnodes);
    // problem.phases(3).guess.states.row(1) = linspace(x2_init * 2 / 4, x2_init * 1 / 4, nnodes);
    // problem.phases(3).guess.states.row(2) = linspace(x3_init * 2 / 4, x3_init * 1 / 4, nnodes);
    // problem.phases(3).guess.states.row(3) = linspace(x4_init * 2 / 4, x4_init * 1 / 4, nnodes);

    // problem.phases(4).guess.states.row(0) = linspace(x1_init * 1 / 4, 0, nnodes);
    // problem.phases(4).guess.states.row(1) = linspace(x2_init * 1 / 4, 0, nnodes);
    // problem.phases(4).guess.states.row(2) = linspace(x3_init * 1 / 4, 0, nnodes);
    // problem.phases(4).guess.states.row(3) = linspace(x4_init * 1 / 4, 0, nnodes);

    double tf = end_time;
    // double tf = 10;

    for (int iphase = 1; iphase != phases + 1; iphase++)
    {
        problem.phases(iphase).guess.time = linspace((iphase - 1) / phases, tf * iphase / phases, nnodes);
    }

    // problem.phases(1).guess.time = linspace(0, tf / 4, nnodes);
    // problem.phases(2).guess.time = linspace(tf / 4, tf * 2 / 4, nnodes);
    // problem.phases(3).guess.time = linspace(tf * 2 / 4, tf * 3 / 4, nnodes);
    // problem.phases(4).guess.time = linspace(tf * 3 / 4, tf, nnodes);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   /////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    MatrixXd x, u, t;
    MatrixXd x_ph, u_ph, t_ph;

    // x_ph1 = solution.get_states_in_phase(1);
    // u_ph1 = solution.get_controls_in_phase(1);
    // t_ph1 = solution.get_time_in_phase(1);

    // x_ph2 = solution.get_states_in_phase(2);
    // u_ph2 = solution.get_controls_in_phase(2);
    // t_ph2 = solution.get_time_in_phase(2);

    // x_ph3 = solution.get_states_in_phase(3);
    // u_ph3 = solution.get_controls_in_phase(3);
    // t_ph3 = solution.get_time_in_phase(3);

    // x_ph4 = solution.get_states_in_phase(4);
    // u_ph4 = solution.get_controls_in_phase(4);
    // t_ph4 = solution.get_time_in_phase(4);

    x.resize(4, nodes * phases);
    u.resize(1, nodes * phases);
    t.resize(1, nodes * phases);

    for (int iphase = 1; iphase != phases + 1; iphase++)
    {
        x_ph = solution.get_states_in_phase(iphase);
        u_ph = solution.get_controls_in_phase(iphase);
        t_ph = solution.get_time_in_phase(iphase);
        x << x_ph;
        u << u_ph;
        t << t_ph;
    }

    // x << x_ph1, x_ph2, x_ph3, x_ph4;
    // u << u_ph1, u_ph2, u_ph3, u_ph4;
    // t << t_ph1, t_ph2, t_ph3, t_ph4;
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
    out << t, x, u;          // pass in the data

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
