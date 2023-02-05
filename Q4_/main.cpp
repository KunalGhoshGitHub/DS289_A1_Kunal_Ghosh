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

// Function euler_implicit solves the ODE using implicit euler method
void euler_implicit(vector<real> &sol, vector<real> &x)
{
    real dx;
    for (int i = 1;i < x.size();i++)
    {
        dx = (x[i] - x[i-1]);
        sol[i] = (sol[i-1] + (dx*(199999*exp(-x[i]))))/(1 + (dx*(200000)));
    }
}

int main()
{
    // Interval X0 <= x <= X
    // Interval starts at: X0
    // Interval ends at: X
    // dx = step size
    real dx,X,X0;

    // num_interval_euler_implicit = Number of intervals for the implicit euler method
    // num_euler_implicit_grid_points = Number of grid points for the implicit euler method
    int num_interval_euler_implicit,num_euler_implicit_grid_points;

    // Vector to read from the input files
    vector<real> input;

    // Vector to store the names of all the output files generated
    vector<string> Output_file_names;

    // "Input.txt" has all the parameters
    input = input_parameters("Input.txt");


    X0 = input[0];
    X = input[1];
    dx = input[3];

    num_interval_euler_implicit = int((X-X0)/dx);

    num_euler_implicit_grid_points = num_interval_euler_implicit+1;

    // t_euler = Vector corresponding to location of the grid points
    vector<real> t_euler(num_euler_implicit_grid_points);

    // Function t_interval(address to the vector of the grid points (t_euler), step size(dt), start of interval(X0), end of interval(X))
    // Writes the location of the grid points to the vector of the grid points (t_euler)
    t_interval(t_euler,dx,X0,X);

    // Vector corresponding to implicit euler solution
    vector<real> sol_euler(num_euler_implicit_grid_points);

    // Initial condition
    // Initial conditions: y(x=0)
    sol_euler[0] = input[2];

    // Function euler_implicit solves the ODE using implicit euler method
    euler_implicit(sol_euler,t_euler);

    // Writes the grid point location in a file
    write_to_file(t_euler,"Q4_b_Grid_Points_Euler_Implicit_Solution.csv");

    Output_file_names.push_back("Q4_b_Grid_Points_Euler_Implicit_Solution.csv");

    // Writes the numerical solution using implicit euler method at different grid point location in a file
    write_to_file(sol_euler,"Q4_b_Euler_Implicit_Solution.csv");

    Output_file_names.push_back("Q4_b_Euler_Implicit_Solution.csv");

    write_to_file(Output_file_names,"Output_file_names.csv");

    return 0;
}
