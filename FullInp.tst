#    This is a sample input-file for InSpecTOR.

nev = 15        # any integer
pot = FooBar    # a string used later for the name of output-files
spect =sm       #  sm, lm, (magnitude) sr, lr, (real) si, li, (imaginary) 
Energy =0.05    # photon energy!
infinite =true  # whether infinite elements should be allowed

cap=false # switches complex absorbing boundaries on and off.
gamma=2   # parameter to scale the cap.

angular=fibonacci  # fibonacci, spiral, archimedes, lebedev, womEV, womMD, fliME, womMN, Sdesign
order = 1
radius= 20.
points= 30
maxVol= 800.
bending= 20.
circles= 5
power= 5
maxiter=700 #default-value.

refine = false  # flag: mesh refinement allowed or not?

mesh_file=Hplus.grid 
mol_file=Hplus.mo
