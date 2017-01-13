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

// Local includes
#include "NN_interpolation.h"
#include "libmesh/point.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"

namespace libMesh
{

//--------------------------------------------------------------------------------
// NeNeInterpolation methods
template <unsigned int KDDim>
void NeNeInterpolation<KDDim>::construct_kd_tree ()
{
#ifdef LIBMESH_HAVE_NANOFLANN

  LOG_SCOPE ("construct_kd_tree()", "NeNeInterpolation<>");

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
void NeNeInterpolation<KDDim>::clear()
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
void NeNeInterpolation<KDDim>::interpolate_field_data (const std::vector<std::string> & field_names,
                                                                  const std::vector<Point> & tgt_pts,
                                                                  std::vector<Number> & tgt_vals) const
{
  libmesh_experimental();

  // forcibly initialize, if needed
#ifdef LIBMESH_HAVE_NANOFLANN
  if (_kd_tree.get() == libmesh_nullptr)
    const_cast<NeNeInterpolation<KDDim> *>(this)->construct_kd_tree();
#endif

  LOG_SCOPE ("interpolate_field_data()", "NeNeInterpolation<>");

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

    //const size_t num_results = (size_t) _src_pts.size();
    const size_t num_results = (size_t) 2;

    std::vector<size_t> ret_index(num_results);
    std::vector<Real>   ret_dist_sqr(num_results);

    for (std::vector<Point>::const_iterator tgt_it=tgt_pts.begin();
         tgt_it != tgt_pts.end(); ++tgt_it)
      {
        const Point & tgt(*tgt_it);
        const Real query_pt[] = { tgt(0), tgt(1), tgt(2) };

        _kd_tree->knnSearch(&query_pt[0], num_results, &ret_index[0], &ret_dist_sqr[0]);

        this->interpolate (tgt, ret_index, ret_dist_sqr, out_it);
        ++out_it;
      }
  }
#else

  libmesh_error_msg("ERROR:  This functionality requires the library to be configured\n" \
                    << "with nanoflann KD-Tree approximate nearest neighbor support!");

#endif
}

// keep it for consistency, actually not needed since there is no real interpolation.
template <unsigned int KDDim>
void NeNeInterpolation<KDDim>::interpolate (const Point               &  pt ,
                                                       const std::vector<size_t> & src_indices,
                                                       const std::vector<Real>   & src_dist_sqr,
                                                       std::vector<Number>::iterator & out_it) const
{
   // number of variables is restricted to 1 here due to rbf_interpol_nd() at the moment.
   libmesh_assert_equal_to (this->n_field_variables(),1);
   // Compute the interpolation weights & interpolated value
   Real threshold=3e-5*3e-5;
   if (src_dist_sqr[0] < threshold){
      *out_it = _src_vals[src_indices[0]];
      return;
   }
   int closestAtom=-1;
   // be always larger than the minimum distance...
   Real minDist=(_geom[0]-pt).norm()+1.;
   for(int atom=0; atom< _geom.size(); atom++){
      if (minDist>(_geom[atom]-pt).norm() ){
          minDist=(_geom[atom]-pt).norm();
          closestAtom=atom;
      }
   }
   // make sure we are close to one nucleus
   //assert((pt-_geom[closestAtom]).norm() < 1e-3);
   //
   if ((pt-_geom[closestAtom]).norm() > 2){
      double pot=0;
      for(int atom=0; atom< _geom.size(); atom++){
         pot += -1./(_geom.size()*(_geom[atom]-pt).norm());
      }
      *out_it=pot;
      return;
   }
   if ((pt-_geom[closestAtom]).norm() > 1e-3)
      err<<"distance to closest point is large "<< src_dist_sqr[0]<<"  "<<(pt-_geom[closestAtom]).norm()<<std::endl;

   // just set it to Z/r:
   *out_it = -_geom[closestAtom].id()/minDist;
   return;

   // fit a+b/r .
   //Real** A;
   //A= new *Real[2];
   //A[0]= new Real[src_dist_sqr.size()];
   //A[1]= new Real[src_dist_sqr.size()];
   
   //Real** pinv;
   //pinv= new *Real[2];
   //pinv[0]= new Real[src_dist_sqr.size()];
   //pinv[1]= new Real[src_dist_sqr.size()];

   //for(unsigned int pt_ind=0; pt_ind<src_dist_sqr.size(); pt_ind++){
   //   A[0][pt_ind]=1;
   //   A[1][pt_ind]=-1/(_src_pts[src_indices[pt_ind]]-_geom[closestAtom]).norm();
   //   pinv[0][pt_ind]=1;
   //   pinv[1][pt_ind]=-1/(_src_pts[src_indices[pt_ind]]-_geom[closestAtom]).norm();
   //}
   //dgesdd(2, src_dist_sqr.size(), pinv)
   //Real x(2)=(pinv.dot(A)*y)

   //*out_it=x[0]+x[1]/min_dist;

   //delete [] A[0];
   //delete [] A[1];
   //delete [] A;
   //delete [] pinv[0];
   //delete [] pinv[1];
   //delete [] pinv;

   //return;

}

// ------------------------------------------------------------
// Explicit Instantiations
template class NeNeInterpolation<1>;
template class NeNeInterpolation<2>;
template class NeNeInterpolation<3>;

} // namespace libMesh
