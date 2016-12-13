#    This is a sample input-file for InSpecTOR.

formulation=symmetric #original,squared, symmetric, root, power
power=0.99      # power of D in the basis functions
nev = 15        # any integer
pot = FooBar    # a string used later for the name of output-files
Energy =0.05    # photon energy (1eV = 0.0367493 )
infinite =true  # whether infinite elements should be allowed
print_quadrature = false # whether quadrature points should be printed or not
pictorious = false # if DO and ESP should be computed explicitly
spherical_analysis= -1 # if a projection onto spherical waves should be performed 
                       # number gives maximum l to be used
cubes=true      # whether or not to make cube files
guessed l=1     # set l and m for the Spherical Harmonic used as initial guess.
guessed m=1     # can be any integers

cap=false # switches complex absorbing boundaries on and off.
gamma=2   # parameter to scale the cap.
offset=1  # distance from the edge to the starting of CAP


angular=fibonacci # fibonacci, spiral, archimedes, lebedev, womEV, womMD, fliME, womMN, Sdesign #, geodesic4, geodesic6
order = 1
radius= 20.
scheme= son      # son, tm, const, quadr, const_tm, sqrt_tm, tm_300, const_son
p=1.0            # if scheme==tm, this option is evaluated.
maxVol= 800.
bending= 20.
circles= 5
r_0= 0.5    # number of neighbours considered in the interpolation. Default: 12
maxiter=700 #default-value.
transform=cayley # cayley, sinv ; else: spectral shift is used.
solver=lapack   # lapack, arnoldi,lanczos ; else: Krylov-Schur is used.

refine = false  # flag: mesh refinement allowed or not?

mesh_file=Hplus.grid 
mol_file=Hplus.mo
