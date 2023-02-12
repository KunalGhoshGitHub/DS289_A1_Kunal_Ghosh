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

// ODE given is mx'' + cx' + kx = 0
// Let's write this equation as two 1st order ODE
// Let x' = y
// So, y' = (-cy/m) + (-kx/m)

// Function corresponding to
// x' = y
real x_prime(real y)
{
    return  y;
}

// Function corresponding to
// y' = (-cy/m) + (-kx/m)
real y_prime(real x, real y, real c, real m, real k)
{
    return (-(c*y/m) - (k*x/m));
}

int main()
{
    // m = mass in kg
    // k = spring constant in N/m
    // Time Interval:
    // Interval starts at: T0
    // Interval ends at: T
    // dt = Time step
    // c = Damping coefficient in Ns/m:
    real dt,T,T0,m,c,k;

    // Vector to read from the input files
    vector<real> input, Values_of_c;

    // Vector to store the names of all the output files generated
    vector<string> Output_file_names;

    // "Input.txt" has all the parameters except the value of c(s)
    input = input_parameters("Input.txt");

    //  "Input_Values_of_c.txt" has the values for different c(s)
    Values_of_c = input_parameters("Input_Values_of_c.txt");

    m = input[0];
    k = input[1];
    T0 = input[2];
    T = input[3];
    dt = input[4];

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
    vector<real> sol1_RK4(num_RK4_grid_points);
    vector<real> sol2_RK4(num_RK4_grid_points);

    // Initial condition

    // Boundary Condition: x(t=0) = 1m
    sol1_RK4[0] = input[5];

    // Boundary Condition: x'(0) = 0 m/s
    sol2_RK4[0] = input[6];

    // Writes the grid point location in a file
    write_to_file(t_RK4,"Q3_Grid_Points_RK4_Solution.csv");

    Output_file_names.push_back("Q3_Grid_Points_RK4_Solution.csv");

    // This for loop will iterate over different values of the c(s)
    for(int i= 0;i<Values_of_c.size();i++)
    {
        // Value of c
        c = Values_of_c[i];

        // Function RK4 solves the ODE using RK4 method for 2 variables
        RK4(sol1_RK4,sol2_RK4,t_RK4, &x_prime,&y_prime,c,m,k);

        // Writes the numerical solutions for 2 variables using RK4 method at different grid point location in a file
        write_to_file(sol1_RK4,"Q3_RK4_1_Solution_c_"+to_string(c)+".csv");
        write_to_file(sol2_RK4,"Q3_RK4_2_Solution_c_"+to_string(c)+".csv");

        Output_file_names.push_back("Q3_RK4_1_Solution_c_"+to_string(c)+".csv");
        Output_file_names.push_back("Q3_RK4_2_Solution_c_"+to_string(c)+".csv");
    }

    real Tol;
    Tol = 0.5;

    // RK4(sol1_RK4,sol2_RK4,t_RK4, &x_prime,&y_prime,c,m,k);
    // This for loop will iterate over different values of the step_size(s)

    c = 5.0;
    vector<real> d_t;

    dt = 1.0;
    real avg_error1, Previous_avg_error;
    for(int i= 0;i<10;i++)
    {
        dt = dt/2.0;

        d_t.push_back(dt);

        cout<<d_t[i]<<endl;
        //cout<<"Average Error: "<<avg_error_vector(sol1_RK4,Previous_Sol1)<<endl;
        num_interval_RK4 = int((T-T0)/d_t[i]);

        num_RK4_grid_points = num_interval_RK4+1;

        // t_RK4 = Location of the grid point
        vector<real> t_RK4(num_RK4_grid_points),sol1_RK4(num_RK4_grid_points),sol2_RK4(num_RK4_grid_points),Previous_Sol1(num_RK4_grid_points);

        // Initial condition

        // Boundary Condition: x(t=0) = 1m
        sol1_RK4[0] = input[5];

        // Boundary Condition: x'(0) = 0 m/s
        sol2_RK4[0] = input[6];

        t_interval(t_RK4,d_t[i],T0,T);
        // Function RK4 solves the ODE using RK4 method for 2 variables
        RK4(sol1_RK4,sol2_RK4,t_RK4, &x_prime,&y_prime,c,m,k);

        // Writes the numerical solutions for 2 variables using RK4 method at different grid point location in a file
        write_to_file(sol1_RK4,"Q3_RK4_1Solution_dt_"+to_string(d_t[i])+".csv");
        write_to_file(sol2_RK4,"Q3_RK4_2Solution_dt_"+to_string(d_t[i])+".csv");

        Output_file_names.push_back("Q3_RK4_1Solution_dt_"+to_string(d_t[i])+".csv");
        Output_file_names.push_back("Q3_RK4_2Solution_dt_"+to_string(d_t[i])+".csv");

        //Previous_Sol1 = sol1_RK4;

        //avg_error1 = avg_error_vector(sol1_RK4,Previous_Sol1);
    }

    write_to_file(d_t,"Values_of_dt.csv");
    Output_file_names.push_back("Values_of_dt.csv");

    write_to_file(Output_file_names,"Output_file_names.csv");

    return 0;
}
