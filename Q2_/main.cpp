#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>

#include"DS289.h"

using namespace std;

// Using the double as a real number
typedef double real;

// The differential equation given: y' = yt^2 -1.1y

// This function calculates y'
double y_prime(real y, real t)
{
    return  ((y*t*t) - (1.1*y));
}

// This function calculates the analytical solution the ODE
void analytical_solution_y(vector<real> &sol,vector<real> &t)
{
        for (int i = 0;i <t.size() ;i++)
    {
        sol[i] = exp((t[i]*t[i]*t[i]/3.0) - (1.1*t[i]));
    }
}

int main()
{
    // dt = Time Step
    // T0 = Start of the interval
    // T = End of the interval
    // Initial Condition: y(0) = 1
    real dt,T,T0,y0;

    // Vector to read from the input files
    vector<real> input, Values_of_Step_size;

    // Vector to store the names of all the output files generated
    vector<string> Output_file_names;

    // "Input.txt" has all the parameters except the step size(s)
    input = input_parameters("Input.txt");

    //  "Values_of_step_size.txt" has the values for different step size(s)
    Values_of_Step_size = input_parameters("Values_of_step_size.txt");

    T0 = input[0];
    T = input[1];
    y0 = input[2];

    // This for loop will iterate over different values of the step size(s)
    for (int i = 0;i<Values_of_Step_size.size();i++)
    {
        // Value of the step size(s)
        dt = Values_of_Step_size[i];

        // num_interval_analytical = Number of intervals for the analytical solution
        // num_analytical_grid_points = Number of grid points for the analytical solution
        int num_interval_analytical,num_analytical_grid_points;

        num_interval_analytical = int((T-T0)/dt);

        num_analytical_grid_points = num_interval_analytical+1;

        // t = Location of the grid points
        vector<real> t(num_analytical_grid_points);

        // Function t_interval(address to the vector of the grid points (t), step size(dt), start of interval(T0), end of interval(T))
        // Writes the location of the grid points to the vector of the grid points (t)
        t_interval(t,dt,T0,T);

        // Vector corresponding to analytical solution
        vector<real> analytical_sol_y(num_analytical_grid_points);

        // Function analytical_solution_y solves the ODE analytically
        analytical_solution_y(analytical_sol_y,t);

        // Writes the grid point location in a file
        write_to_file(t,"Q2_a_Grid_Points_Analytical_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_a_Grid_Points_Analytical_Solution_Step_Size_"+to_string(dt)+".csv");

        // Writes the analytical solution at different grid point location in a file
        write_to_file(analytical_sol_y,"Q2_a_Analytical_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_a_Analytical_Solution_Step_Size_"+to_string(dt)+".csv");

        // num_interval_euler_explicit = Number of intervals for the explicit euler method
        // num_euler_explicit_grid_points = Number of grid points for the explicit euler method
        int num_interval_euler_explicit,num_euler_explicit_grid_points;

        num_interval_euler_explicit = int((T-T0)/dt);

        num_euler_explicit_grid_points = num_interval_euler_explicit+1;

        // t_euler = Location of the grid points
        vector<real> t_euler(num_euler_explicit_grid_points);

        // Function t_interval(address to the vector of the grid points (t_euler), step size(dt), start of interval(T0), end of interval(T))
        // Writes the location of the grid points to the vector of the grid points (t_euler)
        t_interval(t_euler,dt,T0,T);

        // Vector corresponding to explicit euler solution
        vector<real> sol_euler(num_euler_explicit_grid_points);

        // Initial condition
        sol_euler[0] = y0;

        // Function euler_explicit solves the ODE using explicit euler method
        euler_explicit(sol_euler,t_euler, &y_prime);

        // Writes the grid point location in a file
        write_to_file(t_euler,"Q2_b_i_Grid_Points_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_b_i_Grid_Points_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        // Writes the numerical solution using explicit euler method at different grid point location in a file
        write_to_file(sol_euler,"Q2_b_i_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_b_i_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        // Vectors corresponding to average error and max error for explicit euler method
        vector<real> Average_Error_Euler_Explicit,Max_Error_Euler_Explicit;

        // Function avg_error_vector(analytical_sol,numerical_sol)
        // Calculates the average error between the analytical solution and numerical solution
        Average_Error_Euler_Explicit.push_back(avg_error_vector(analytical_sol_y,sol_euler));

        // Writes the average error in a file
        write_to_file(Average_Error_Euler_Explicit,"Q2_c_Average_error_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_c_Average_error_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        // Function max_error_vector(analytical_sol,numerical_sol)
        // Calculates the max error between the analytical solution and numerical solution
        Max_Error_Euler_Explicit.push_back(max_error_vector(analytical_sol_y,sol_euler));

        // Writes the max error in a file
        write_to_file(Max_Error_Euler_Explicit,"Q2_c_Max_Error_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_c_Max_Error_Euler_Explicit_Solution_Step_Size_"+to_string(dt)+".csv");

        // num_interval_Adam = Number of intervals for the second-order Adams-Bashforth method
        // num_Adam_grid_points = Number of grid points for the second-order Adams-Bashforth method
        int num_interval_Adam,num_Adam_grid_points;

        num_interval_Adam = int((T-T0)/dt);

        num_Adam_grid_points = num_interval_Adam+1;

        // t_Adam = Location of the grid points
        vector<real> t_Adam(num_Adam_grid_points);

        // Function t_interval(address to the vector of the grid points (t_Adam), step size(dt), start of interval(T0), end of interval(T))
        // Writes the location of the grid points to the vector of the grid points (t_Adam)
        t_interval(t_Adam,dt,T0,T);

        // Vector corresponding to second-order Adams-Bashforth solution
        vector<real> sol_Adam(num_Adam_grid_points);

        // Initial condition
        sol_Adam[0] = y0;

        // Solving the 1st step using euler explicit method
        sol_Adam[1] = euler_explicit(sol_Adam[0],1,t_Adam,&y_prime);

        // Function adam solves the ODE using second-order Adams-Bashforth method
        adam(sol_Adam,t_Adam, &y_prime);

        // Writes the grid point location in a file
        write_to_file(t_Adam,"Q2_b_ii_Grid_Points_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_b_ii_Grid_Points_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        // Writes the numerical solution using second-order Adams-Bashforth method at different grid point location in a file
        write_to_file(sol_Adam,"Q2_b_ii_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_b_ii_Adam_Solution_Step_Size_"+to_string(dt)+".csv");
        // Vectors corresponding to average error and max error for second-order Adams-Bashforth method
        vector<real> Average_Error_Adam,Max_Error_Adam;

        // Function avg_error_vector(analytical_sol,numerical_sol)
        // Calculates the average error between the analytical solution and numerical solution
        Average_Error_Adam.push_back(avg_error_vector(analytical_sol_y,sol_Adam));

        // Writes the average error in a file
        write_to_file(Average_Error_Adam,"Q2_c_Average_error_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_c_Average_error_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        // Function max_error_vector(analytical_sol,numerical_sol)
        // Calculates the max error between the analytical solution and numerical solution
        Max_Error_Adam.push_back(max_error_vector(analytical_sol_y,sol_Adam));

        // Writes the max error in a file
        write_to_file(Max_Error_Adam,"Q2_c_Max_Error_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_c_Max_Error_Adam_Solution_Step_Size_"+to_string(dt)+".csv");

        // num_interval_RK4 = Number of intervals for the RK4 method
        // num_RK4_grid_points = Number of grid points for the RK4 method
        int num_interval_RK4,num_RK4_grid_points;

        num_interval_RK4 = int((T-T0)/dt);

        num_RK4_grid_points = num_interval_RK4+1;

        // t_RK4 = Location of the grid point
        vector<real> t_RK4(num_RK4_grid_points);

        // Function t_interval(address to the vector of the grid points (t_RK4), step size(dt), start of interval(T0), end of interval(T))
        // Writes the location of the grid points to the vector of the grid points (t_RK4)
        t_interval(t_RK4,dt,T0,T);

        // Vector corresponding to RK4 solution
        vector<real> sol_RK4(num_RK4_grid_points);

        // Initial condition
        sol_RK4[0] = y0;

        // Function RK4 solves the ODE using RK4 method
        RK4(sol_RK4,t_RK4, &y_prime);

        // Writes the grid point location in a file
        write_to_file(t_RK4,"Q2_b_iii_Grid_Points_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_b_iii_Grid_Points_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        // Writes the numerical solution using RK4 method at different grid point location in a file
        write_to_file(sol_RK4,"Q2_b_iii_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_b_iii_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        // Vectors corresponding to average error and max error for RK4 method
        vector<real> Average_Error_RK4,Max_Error_RK4;

        // Function avg_error_vector(analytical_sol,numerical_sol)
        // Calculates the average error between the analytical solution and numerical solution
        Average_Error_RK4.push_back(avg_error_vector(analytical_sol_y,sol_RK4));

        // Writes the average error in a file
        write_to_file(Average_Error_RK4,"Q2_c_Average_error_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_c_Average_error_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        // Function max_error_vector(analytical_sol,numerical_sol)
        // Calculates the max error between the analytical solution and numerical solution
        Max_Error_RK4.push_back(max_error_vector(analytical_sol_y,sol_RK4));

        // Writes the max error in a file
        write_to_file(Max_Error_RK4,"Q2_c_Max_Error_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        Output_file_names.push_back("Q2_c_Max_Error_RK4_Solution_Step_Size_"+to_string(dt)+".csv");

        write_to_file(Output_file_names,"Output_file_names.csv");
    }

    return 0;
}
