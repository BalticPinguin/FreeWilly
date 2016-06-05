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
#include "libmesh/point.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "radial_interpolation.h"


namespace libMesh
{

//--------------------------------------------------------------------------------
// RadialInterpolation methods
template <unsigned int KDDim>
void RadialInterpolation<KDDim>::construct_kd_tree ()
{
#ifdef LIBMESH_HAVE_NANOFLANN

  LOG_SCOPE ("construct_kd_tree()", "RadialInterpolation<>");

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
void RadialInterpolation<KDDim>::clear()
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
void RadialInterpolation<KDDim>::interpolate_field_data (const std::vector<std::string> & field_names,
                                                                  const std::vector<Point> & tgt_pts,
                                                                  std::vector<Number> & tgt_vals) const
{
  libmesh_experimental();

  // forcibly initialize, if needed
#ifdef LIBMESH_HAVE_NANOFLANN
  if (_kd_tree.get() == libmesh_nullptr)
    const_cast<RadialInterpolation<KDDim> *>(this)->construct_kd_tree();
#endif

  LOG_SCOPE ("interpolate_field_data()", "RadialInterpolation<>");

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

    // how many points to I want to include in the geometry search?
    const size_t num_results = this.geom.size()+1;

    std::vector<size_t> ret_index(num_results);
    std::vector<Real>   ret_dist_sqr(num_results);
    
    // loop over target points and compute the potential there;
    for (std::vector<Point>::const_iterator tgt_it=tgt_pts.begin();
         tgt_it != tgt_pts.end(); ++tgt_it)
      {
        const Point & tgt(*tgt_it);
        const Real query_pt[] = { tgt(0), tgt(1), tgt(2) };

        _kd_tree->knnSearch(&query_pt[0], num_results, &ret_index[0], &ret_dist_sqr[0]);
        
        //tgt:          point where I seach the result
        //ret_index:    indices of the neighbouring elements
        //ret_dist_sqr: distance squared to these points
        //out_it:       potential at tgt (to be computed).
        this->interpolate (tgt, ret_index, ret_dist_sqr, out_it);

        // libMesh::out << "knnSearch(): num_results=" << num_results << "\n";
        // for (size_t i=0;i<num_results;i++)
        //   libMesh::out << "idx[" << i << "]="
        //       << std::setw(6) << ret_index[i]
        //       << "\t dist["<< i << "]=" << ret_dist_sqr[i]
        //       << "\t val[" << std::setw(6) << ret_index[i] << "]=" << _src_vals[ret_index[i]]
        //       << std::endl;
        // libMesh::out << "\n";

        // libMesh::out << "ival=" << _vals[0] << '\n';
      }
  }
#else

  libmesh_error_msg("ERROR:  This functionality requires the library to be configured\n" \
                    << "with nanoflann KD-Tree approximate nearest neighbor support!");

#endif
}

template <unsigned int KDDim>
void RadialInterpolation<KDDim>::interpolate (const Point               & /* pt */,
                                                       const std::vector<size_t> & src_indices,
                                                       const std::vector<Real>   & src_dist_sqr,
                                                       std::vector<Number>::iterator & out_it) const
{
  // We explicitly assume that the input source points are sorted from closest to
  // farthests.  assert that assumption in DEBUG mode.
#ifdef DEBUG
  if (!src_dist_sqr.empty())
    {
      Real min_dist = src_dist_sqr.front();
      std::vector<Real>::const_iterator it = src_dist_sqr.begin();

      for (++it; it!= src_dist_sqr.end(); ++it)
        {
          if (*it < min_dist)
            libmesh_error_msg(*it << " was less than min_dist = " << min_dist);

          min_dist = *it;
        }
    }
#endif


  libmesh_assert_equal_to (src_dist_sqr.size(), src_indices.size());

  // Compute the interpolation weights & interpolated value
  const unsigned int n_fv = this->n_field_variables();
  _vals.resize(n_fv); /**/ std::fill (_vals.begin(), _vals.end(), Number(0.));

  Real V_0 = 0.;

  std::vector<Real>::const_iterator src_dist_sqr_it=src_dist_sqr.begin();
  std::vector<size_t>::const_iterator src_idx_it=src_indices.begin();

  // Loop over source points
  while ((src_dist_sqr_it != src_dist_sqr.end()) &&
         (src_idx_it      != src_indices.end())){
      for (unsigned int i=0; i<this.geom.size(); i++){
   
         const Real dist_sq = std::max(*src_dist_sqr_it, std::numeric_limits<Real>::epsilon()),
         out_it=[src_indices[i]]=   
   
      }
         ++src_dist_sqr_it;
         ++src_idx_it;
    {

}



// ------------------------------------------------------------
// Explicit Instantiations
template class RadialInterpolation<1>;
template class RadialInterpolation<2>;
template class RadialInterpolation<3>;

} // namespace libMesh
