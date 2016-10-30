# include "SlepcConfig.h"

using namespace libMesh;

void SlepcSolverConfiguration::configure_solver()
{
   PetscErrorCode ierr = 0;

   // if a spectral transformation was requested
   if (_st!=INVALID_ST){
    
      // initialise the st with the default values
      //(than, change only the spectral transformation value).
      ST st;
      ierr = EPSGetST(_slepc_solver.eps(), &st);
      libmesh_assert(ierr == 0);
      //STCreate(_slepc_solver.comm().get(), &st);

      // Set it to the desired type of spectral transformation.
      // The value of the respective shift is chosen to be the target
      // specified via \p set_position_of_spectrum().
      switch (_st)
         {
         case SHIFT:
            ierr = STSetType(st, STSHIFT);
            break;
         case SINVERT:
            ierr = STSetType(st, STSINVERT);
            break;
         case CAYLEY:
      #if SLEPC_VERSION_LESS_THAN(2,2,1)
            libmesh_error_msg("SLEPc 2.2.1 is required to call CAYLEY transform.");
            break;
      #else
            ierr = STSetType(st, STCAYLEY);
            break;
      #endif
         default:
            // print a warning but do nothing more.
            break;
         }  //tell the \p EPS object which \p ST to use
      // this is not needed because it is called in the
      // in the \p EPSSetUP() anyway.
      //ierr = EPSSetST(_slepc_solver.eps(), st);

      libmesh_assert(ierr == 0);
   }
}

