// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

    //Seed positions
    double center[120][3]={{509.59, 248.77, 476.04},{33.21, 338.30, 81.49},{505.05, 306.12, 386.70},{226.61, 344.81, 330.93},
        {470.98, 354.67, 9.09},{7.70, 87.05, 394.12},{59.81, 98.62, 153.67},{357.53, 21.67, 316.08},
        {390.62, 280.59, 509.62},{493.88, 349.31, 94.46},{151.34, 398.41, 184.62},{215.34, 443.33, 412.47},
        {253.00, 102.03, 106.28},{320.61, 175.23, 201.30},{333.23, 176.50, 293.29},{332.83, 382.11, 298.69},
        {206.28, 68.49, 449.61},{424.31, 410.21, 228.35},{445.89, 284.08, 196.18},{340.77, 13.31, 84.00},
        {136.24, 267.18, 51.20},{129.68, 302.53, 272.27},{450.60, 114.83, 203.23},{206.03, 443.10, 94.09},
        {295.40, 194.79, 26.92},{462.75, 324.92, 201.27},{275.99, 274.60, 82.95},{379.44, 259.23, 105.88},
        {98.27, 80.46, 229.12},{442.46, 494.70, 211.65},{213.85, 236.15, 246.96},{248.51, 82.83, 349.01},
        {66.18, 373.04, 7.85},{227.32, 461.56, 298.16},{246.02, 365.39, 258.19},{362.67, 505.84, 427.99},
        {30.04, 58.38, 173.07},{190.55, 288.57, 463.08},{100.37, 342.05, 263.20},{496.79, 86.13, 398.82},
        {137.22, 284.06, 356.56},{198.30, 94.83, 491.90},{251.81, 131.74, 330.13},{201.45, 99.94, 348.80},
        {468.24, 38.98, 433.17},{280.28, 73.26, 280.18},{226.30, 365.98, 305.74},{363.70, 510.94, 221.64},
        {213.80, 434.50, 24.62},{172.52, 323.99, 359.97},{289.72, 25.45, 185.96},{187.06, 258.15, 454.60},
        {33.85, 145.24, 470.89},{397.28, 151.43, 203.69},{328.13, 90.33, 313.89},{302.63, 479.87, 511.41},
        {320.57, 463.92, 213.91},{327.59, 9.28, 332.81},{65.02, 306.37, 207.56},{249.63, 309.61, 452.78},
        {226.91, 396.50, 325.37},{82.46, 227.18, 379.54},{193.53, 261.87, 40.30},{119.73, 161.37, 48.76},
        {463.48, 111.20, 66.75},{262.07, 89.41, 468.33},{400.34, 25.58, 259.35},{162.04, 169.50, 239.76},
        {82.27, 311.04, 481.46},{475.16, 72.37, 63.28},{362.58, 511.76, 494.73},{466.17, 472.55, 332.13},
        {85.41, 423.50, 8.42},{275.90, 472.40, 174.36},{224.58, 175.34, 254.22},{377.78, 383.74, 218.98},
        {268.75, 354.43, 93.47},{25.00, 238.73, 429.98},{60.59, 86.90, 438.88},{1.25, 103.69, 307.21},
        {365.88, 234.14, 102.43},{471.41, 485.72, 134.76},{267.93, 51.01, 324.24},{33.41, 182.77, 114.63},
        {190.82, 211.87, 305.24},{253.32, 157.66, 411.65},{387.51, 310.94, 219.02},{113.66, 135.17, 2.72},
        {249.47, 303.96, 223.74},{175.00, 431.25, 100.14},{13.24, 478.95, 20.52},{261.46, 400.65, 406.22},
        {324.76, 257.73, 312.73},{156.94, 112.84, 277.37},{83.84, 500.20, 270.34},{499.68, 119.76, 200.91},
        {270.74, 158.39, 170.34},{375.81, 437.93, 303.11},{247.75, 226.58, 219.46},{217.48, 148.51, 103.60},
        {379.53, 233.84, 68.53},{489.12, 2.02, 143.02},{128.26, 82.94, 207.24},{145.10, 437.10, 175.48},
        {308.96, 509.69, 426.35},{11.79, 353.99, 193.32},{191.12, 202.13, 414.09},{400.19, 464.48, 485.37},
        {193.33, 21.11, 451.19},{363.82, 146.41, 334.10},{240.53, 411.41, 64.46},{424.69, 419.76, 115.34},
        {139.63, 487.88, 473.73},{241.77, 325.63, 256.00},{427.37, 168.68, 94.95},{81.09, 144.36, 89.22},
        {419.34, 457.06, 184.30},{53.42, 119.98, 206.36},{282.25, 60.10, 482.38},{299.56, 90.00, 308.74}};

    //Seed radii
    double rad[120]={4.539, 4.893, 4.711, 5.636,
        5.203, 5.317, 4.475, 4.874,
        5.286, 4.782, 4.584, 4.982,
        4.901, 5.877, 4.911, 5.406,
        4.761, 4.692, 5.020, 4.803,
        4.706, 4.560, 5.123, 4.791,
        5.168, 4.871, 4.863, 4.848,
        5.314, 5.331, 5.004, 5.201,
        4.457, 4.855, 4.839, 3.966,
        5.508, 4.942, 4.847, 5.080,
        4.221, 5.210, 4.758, 5.585,
        4.661, 4.073, 5.023, 5.071,
        4.463, 3.995, 5.546, 5.486,
        4.948, 5.312, 4.425, 4.546,
        4.716, 5.993, 4.903, 5.187,
        4.761, 4.547, 4.646, 4.088,
        4.903, 4.922, 5.052, 5.038,
        5.183, 5.024, 5.515, 5.403,
        5.474, 5.016, 5.436, 4.577,
        5.474, 3.704, 5.230, 5.374,
        5.012, 5.230, 4.552, 4.389,
        4.813, 4.323, 6.144, 5.038,
        4.886, 5.120, 5.411, 5.096,
        4.423, 6.045, 5.075, 5.498,
        4.953, 4.520, 5.149, 4.929,
        3.763, 4.793, 4.832, 4.588,
        5.205, 4.984, 5.418, 5.295,
        4.522, 4.348, 4.435, 4.584,
        4.902, 4.837, 4.871, 5.243,
        4.657, 4.582, 4.638, 5.875};
    
    double xi=std::sqrt(KnV/W_barrier);
    double interface_coeff=std::sqrt(2.0)*xi;
    double dist;
    scalar_IC = 0.0;
    
    double n0 = 0.0;
    for (unsigned int i=0; i<120; i++){
        dist = 0.0;
        for (unsigned int dir = 0; dir < dim; dir++){
            double mindist=std::min(std::abs(p[dir]-center[i][dir]),userInputs.domain_size[dir]-std::abs(p[dir]-center[i][dir]));
            dist += mindist*mindist;
        }
        dist = std::sqrt(dist);
        n0 += 0.5*(1.0-std::tanh((dist-rad[i])/interface_coeff));
    }
	  // Initial condition for the concentration field
	  if (index == 0){
        scalar_IC = c_avg+n0*(cbtmin-c_avg);
	  }
	  // Initial condition for the structural order parameter field
	  else {
        scalar_IC = n0;
	  }
	  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
