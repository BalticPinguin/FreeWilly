#ifndef SLEPC_CONFIGSOLVER_H
#define SLEPC_CONFIGSOLVER_H

# include "libmesh/solver_configuration.h"
# include "libmesh/slepc_eigen_solver.h"

EXTERN_C_FOR_SLEPC_BEGIN
# include <slepceps.h>
EXTERN_C_FOR_SLEPC_END

/**
 * Defines an \p enum for spectral tronsformations
 * applied before solving the (generalised) eigenproblem
 */
enum SpectralTransform {SHIFT=0,
                        SINVERT,
                        CAYLEY,

                        INVALID_ST
};

class SlepcSolverConfiguration : public libMesh::SolverConfiguration
{
public:

   SlepcSolverConfiguration( libMesh::SlepcEigenSolver<libMesh::Number> & slepc_eigen_solver):
        _slepc_solver(slepc_eigen_solver),
        _st(INVALID_ST)
   {}
   
   ~SlepcSolverConfiguration() {}

   virtual void configure_solver() override;
   
   void SetInterval(libMesh::Real lower, libMesh::Real upper)
   { _upper=upper; _lower=lower;} 

   void SetST(SpectralTransform st)
   { _st=st;}
   
private:
   // The linear solver object that we are configuring
   libMesh::SlepcEigenSolver<libMesh::Number>& _slepc_solver;
   SpectralTransform _st;
   //ST st;
   libMesh::Real _upper=0,_lower=0;

};

#endif // define SLEPC_CONFIGSOLVER_H
