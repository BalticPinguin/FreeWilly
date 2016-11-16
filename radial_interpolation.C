// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// C++ includes
#include <iomanip>

// include functions for rbf-interpolation
# include "fsu_soft/rbf_interp_nd.hpp"
# include "fsu_soft/r8lib.hpp"

// Local includes
#include "radial_interpolation.h"
#include "libmesh/point.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"

namespace libMesh
{

//--------------------------------------------------------------------------------
// RBFInterpolation methods
template <unsigned int KDDim>
void RBFInterpolation<KDDim>::construct_kd_tree ()
{
#ifdef LIBMESH_HAVE_NANOFLANN

   LOG_SCOPE ("construct_kd_tree()", "RBFInterpolation<>");
   
   // Initialize underlying KD tree
   if (_kd_tree.get() == libmesh_nullptr)
    _kd_tree.reset (new kd_tree_t (KDDim,
                                   _point_list_adaptor,
                                   nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)));
   
   libmesh_assert (_kd_tree.get() != libmesh_nullptr);
   
   _kd_tree->buildIndex();
#endif
}

template <unsigned int KDDim>
void RBFInterpolation<KDDim>::clear()
{
#ifdef LIBMESH_HAVE_NANOFLANN
   // Delete the KD Tree and start fresh
   if (_kd_tree.get())
      _kd_tree.reset (libmesh_nullptr);
#endif
   
   // Call  base class clear method
   MeshfreeInterpolation::clear();
}


template <unsigned int KDDim>
void RBFInterpolation<KDDim>::interpolate_field_data (const std::vector<std::string> & field_names,
                                                                  const std::vector<Point> & tgt_pts,
                                                                  std::vector<Number> & tgt_vals) const
{
   libmesh_experimental();
      
   // forcibly initialize, if needed
#ifdef LIBMESH_HAVE_NANOFLANN
   if (_kd_tree.get() == libmesh_nullptr)
      const_cast<RBFInterpolation<KDDim> *>(this)->construct_kd_tree();
#endif
   
   LOG_SCOPE ("interpolate_field_data()", "RBFInterpolation<>");
   
   libmesh_assert_equal_to (field_names.size(), this->n_field_variables());
   
   // If we already have field variables, we assume we are appending.
   // that means the names and ordering better be identical!
   if (_names.size() != field_names.size())
      libmesh_error_msg("ERROR:  when adding field data to an existing list the \nvariable list must be the same!");
   
   for (unsigned int v=0; v<_names.size(); v++)
      if (_names[v] != field_names[v])
         libmesh_error_msg("ERROR:  when adding field data to an existing list the \nvariable list must be the same!");
   
   tgt_vals.resize (tgt_pts.size()*this->n_field_variables());
    
#ifdef LIBMESH_HAVE_NANOFLANN
   {
      std::vector<Number>::iterator out_it = tgt_vals.begin();
   
      const size_t num_results = std::min((size_t) _n_interp_pts, _src_pts.size());
   
      std::vector<size_t> ret_index(num_results);
      std::vector<Real>   ret_dist_sqr(num_results);
   
      for (std::vector<Point>::const_iterator tgt_it=tgt_pts.begin();
            tgt_it != tgt_pts.end(); ++tgt_it)
         {
         const Point & tgt(*tgt_it);
         const Real query_pt[] = { tgt(0), tgt(1), tgt(2) };

         _kd_tree->knnSearch(&query_pt[0], num_results, &ret_index[0], &ret_dist_sqr[0]);

         // now check that points are not 2-dimensional only; this can lead to errors.
         double dist;
         for(unsigned int dim=0; dim<KDDim; dim++){
            dist=0;
            for(unsigned int i=0; i<ret_index.size(); i++){
               dist+=(_src_pts[ret_index[i]](dim)-_src_pts[ret_index[0]](dim))*
                     (_src_pts[ret_index[i]](dim)-_src_pts[ret_index[0]](dim));
            }
            if (dist<1e-10){
               //this seems to be a reasonable threshold. Now, search close to nearest point for 
               // its neighbours. This should not give the same problem again.
              const Real query_pt[] = { _src_pts[ret_index[0]](0), 
                                        _src_pts[ret_index[0]](1), 
                                        _src_pts[ret_index[0]](2) };
              _kd_tree->knnSearch(&query_pt[0], num_results, &ret_index[0], &ret_dist_sqr[0]);
              break;
            }
         }
         this->interpolate (tgt, ret_index, ret_dist_sqr, out_it);
        }
   }
#else
 
   libmesh_error_msg("ERROR:  This functionality requires the library to be configured\n" \
                     << "with nanoflann KD-Tree approximate nearest neighbor support!");
 
#endif
}

template <unsigned int KDDim>
void RBFInterpolation<KDDim>::interpolate (const Point               &  pt ,
                                                       const std::vector<size_t> & src_indices,
                                                       const std::vector<Real>   & /*src_dist_sqr*/,
                                                       std::vector<Number>::iterator & out_it) const
{
   // number of variables is restricted to 1 here due to rbf_interpol_nd() at the moment.
   libmesh_assert_equal_to (this->n_field_variables(),1);
   // Compute the interpolation weights & interpolated value
   //const unsigned int n_fv = this->n_field_variables();
   const unsigned int n_fv = 1;
   _vals.resize(n_fv); /**/ std::fill (_vals.begin(), _vals.end(), Number(0.));
  
   // xd= point values 
   const unsigned int n_src=src_indices.size();
   Real* xd= new Real[n_src*KDDim];
   Real* fd= new Real[n_src*KDDim];
   Real max_fd=0;
   Real maxDist=0;
   Real min_fd=_src_vals[src_indices[0]].real();
   // convert (input) vector to pointer-notation for rbf-functions:
   for (unsigned int i=0; i<n_src; i++){
      for (unsigned int j=0; j<KDDim; j++){
          xd[j+KDDim*i]=_src_pts[src_indices[i]](j)-pt(j);
      }
      for (unsigned int j=0; j<n_src; j++){
         if(maxDist<(_src_pts[src_indices[i]]-_src_pts[src_indices[j]]).norm())
            maxDist=(_src_pts[src_indices[i]]-_src_pts[src_indices[j]]).norm();
      }
      fd[i]=_src_vals[src_indices[i]].real();
      if(fd[i]>max_fd)
         max_fd=fd[i];
      if(fd[i]<min_fd)
         min_fd=fd[i];
   }
   
   int closest=-1;
   double minnorm=1e42; // larger than distance to all atoms
   for (unsigned int i=0; i<_geom.size(); i++){
     if ((_geom[i]-pt).norm()<minnorm){
        closest=i;
        minnorm=(_geom[i]-pt).norm();
     }
   }
 
   //Real r0=_power;
   Real r0=1.3*maxDist;  // I hope this is a reasonable number...
   Real inner_range=1.2; // about half of the typical binding length...
   //Real outer_range=4.;  // 'far away' from any nucleus.
         
   Real *w;
   /*
    xi: points where I want to know the values?
   */
   Real* xi= new Real[KDDim];//{0,0,0};  // interpolation point is '0' always
   for(unsigned int i=0; i<KDDim; i++){
      xi[i]=0;
   }
   unsigned int ni = 1;
   Real* fi;

   if (minnorm<inner_range){
      for (unsigned int i=0; i<n_src; i++){
         Real r=(_src_pts[src_indices[i]]-_geom[closest]).norm();
         if (r>1e-7)
            fd[i]+=(Real)_geom[closest].id()/r;
         else
            fd[i]+=(Real)_geom[closest].id()*1e+7;
      }
      //r0*=0.3;
   
      w = rbf_weight (KDDim, n_src, xd, r0, phi1, fd );
      fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi1, w, ni, xi );
      Real r=(pt-_geom[closest]).norm();
      for (unsigned int v=0; v<n_fv; v++){
         if (r>1e-7)
            fi[v]-=(Real)_geom[closest].id()/r;
         else
            fi[v]-=(Real)_geom[closest].id()*1e+7;
      }
   }
   //else if (minnorm > outer_range){
   //   //r0=0.7;
   //   // if the distance to the closest source point is too large: interpolate as
   //   // simple Coulomb-pot. of nearest atom.
   //   if (src_dist_sqr[0]>1.5){
   //      for (unsigned int v=0; v<n_fv; v++, ++out_it){
   //         _vals[v] = fd[0]/(pt-_geom[closest]).norm()*
   //                    (_src_pts[src_indices[0]]-_geom[closest]).norm();
   //         *out_it = _vals[v];
   //      }
   //      return ;
   //   }
   //   w = rbf_weight (KDDim, n_src, xd, r0, phi1, fd );
   //   fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi1, w, ni, xi );
   //}
   else{
      // intermediate case: don't change parameters.

      try{
       w = rbf_weight (KDDim, n_src, xd, r0, phi1, fd );
      }
      // catch the error-case and compute the mean of NN instead...
      catch (int exc){
         // in the case of SVD-failure, use the weighted mean.
         libmesh_assert (exc == 1);
         libmesh_warning("rbf")
         Real fi;
         const unsigned int four=4;
         if (src_dist_sqr[0]>r0){
            // point is outside of the mesh
            for (unsigned int v=0; v<n_fv; v++, ++out_it){
               _vals[v] = fd[0]/(pt-_geom[closest]).norm()*
                        (_src_pts[src_indices[0]]-_geom[closest]).norm();
               *out_it = _vals[v];
            }
            return ;
         }

         // if point is inside the mesh: compute the mean...
         for(unsigned int i=0; i<std::min(n_src, four); i++){
            fi+=fd[i];
         }

         //just in case soemthing went wrong...
         if(fi>=max_fd*n_src || fi<=min_fd*n_src)
            fi=(max_fd+min_fd)/2;
         _vals[0]=fi/n_src;
         *out_it=_vals[0];
         return ;
      }
      fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi1, w, ni, xi );

      // if the interpolation failed: just take average of extremal points.
      // Usually it only failes in the outer part of the mesh where both are 
      // 0 anyways.
      for (unsigned int v=0; v<n_fv; v++)
         {
         if(fi[v]>=max_fd || fi[v]<=min_fd)
            fi[v]=(max_fd+min_fd)/2;
         }
   }
   //w = rbf_weight (KDDim, n_src, xd, r0, phi1, fd );
   //fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi1, w, ni, xi );

   delete xi;
   delete w;
   
      
   for (unsigned int v=0; v<n_fv; v++, ++out_it)
   {
      _vals[v] = fi[v];
      *out_it = _vals[v];
   }
   delete fi;
}
      
// ------------------------------------------------------------
// Explicit Instantiations
template class RBFInterpolation<1>;
template class RBFInterpolation<2>;
template class RBFInterpolation<3>;

} // namespace libMesh
