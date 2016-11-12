#    This is a sample input-file for InSpecTOR.

nev = 15        # any integer
pot = FooBar    # a string used later for the name of output-files
#spect =sm      # This option is no more possible.
Energy =0.05    # photon energy!
infinite =true  # whether infinite elements should be allowed
print_quadrature = false # whether quadrature points should be printed or not
pictorious = false # if DO and ESP should be computed explicitly
spherical_analysis= -1 # if a projection onto spherical waves should be performed 
                       # number gives maximum l to be used

cap=false # switches complex absorbing boundaries on and off.
gamma=2   # parameter to scale the cap.

angular=fibonacci # fibonacci, spiral, archimedes, lebedev, womEV, womMD, fliME, womMN, Sdesign #, geodesic4, geodesic6
order = 1
radius= 20.
scheme= son      # son, tm, const, quadr, const_tm, sqrt_tm, tm_300
p=1.0            # if scheme==tm, this option is evaluated.
maxVol= 800.
bending= 20.
circles= 5
r_0= 0.5    # number of neighbours considered in the interpolation. Default: 12
maxiter=700 #default-value.

refine = false  # flag: mesh refinement allowed or not?

mesh_file=Hplus.grid 
mol_file=Hplus.mo
