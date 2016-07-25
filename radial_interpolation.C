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
# include "rbf_interp_nd.hpp"
# include "r8lib.hpp"

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

        //inserted by HUBERT
         //out<<"query point: ";
        //out<<tgt;
        //out<<std::endl;
        //out<<"interpolation_points: ";
        //for(unsigned int aoe=0; aoe<ret_index.size(); aoe++){
        //   out<<_src_pts[ret_index[aoe]];
        //   out<<_src_vals[ret_index[aoe]]<<std::endl;
        //}
        //out<<std::endl;
        //inserted by HUBERT

        this->interpolate (tgt, ret_index, ret_dist_sqr, out_it);

       // libMesh::out << "knnSearch(): num_results=" << num_results << "\n";
       // for (size_t i=0;i<num_results;i++)
       //   libMesh::out << "pt[" << i << "]="
       //       << std::setw(6) <<_src_pts[ret_index[i]]
       //       << "\t dist["<< i << "]=" << ret_dist_sqr[i]
       //       << "\t val[" << std::setw(6) << ret_index[i] << "]=" << _src_vals[ret_index[i]]
       //       << std::endl;
       // libMesh::out << "\n";
       // libMesh::out << "ipt=" << &tgt;
       // libMesh::out << "\t ival=" << _vals[0] << '\n';

      }
        //inserted by HUBERT
        //out<<"value: ";
        //for(unsigned int aoe=0; aoe<tgt_vals.size(); aoe++){
        //   out<<tgt_pts[aoe];
        //   out<<tgt_vals[aoe];
        //   out<<std::endl;
        //}
        //inserted by HUBERT
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
  // convert (input) vector to pointer-notation for rbf-functions:
  //out<<"scr-points:  \n";
  for (unsigned int i=0; i<n_src; i++){
     //out<<_src_pts[src_indices[i]]<<"  ";
     for (unsigned int j=0; j<KDDim; j++){
         xd[j+KDDim*i]=_src_pts[src_indices[i]](j)-pt(j);
         //out<<xd[j+KDDim*i]<<" ";
     }
     fd[i]=_src_vals[src_indices[i]].real();
     //out<<"  ";
     //out<<fd[i]<<std::endl;
  }
  //out<<"done.";
  
  int closest=-1;
  double minnorm=1e42; // larger than distance to all atoms
  for (unsigned int i=0; i<_geom.size(); i++){
     if ((_geom[i]-pt).norm()<minnorm){
        closest=i;
        minnorm=(_geom[i]-pt).norm();
     }
  }
  //out<<"min norm: "<<minnorm;
  //out<<" closest: "<<closest<<std::endl;

  Real r0=_power;
  Real inner_range=1.2; // about half of the typical binding length...
  if (minnorm<inner_range){
      for (unsigned int i=0; i<n_src; i++){
         fd[i]+=(Real)_geom[closest].id()/(_src_pts[src_indices[i]]-_geom[closest]).norm();
         //out<<"distance: "<< (_src_pts[src_indices[i]]-_geom[closest]).norm()<<std::endl;
      }
      r0*=0.2;
  }
  else if (minnorm > 8.){
      r0=0.6;
  }
  else{
      // intermediate case.
  }

  // Compute the interpolation weights & interpolated value
  //---> this needs to be changed accordingly. HUBERT

  /*  m:  dimension (?)
      nd: #points (?)
      xd: points (?)
      fd: function values at points
      r0: parameter.
      phi1: function type (brought from rdb-intperp
      xi: where I want to know point?
  */
  Real *w;
  //w = rbf_weight ( m, nd, xd, r0, phi1, fd );
  if (minnorm<inner_range){
     w = rbf_weight (KDDim, n_src, xd, r0, phi1, fd );
  }
  else if (minnorm > 8.){
     w = rbf_weight (KDDim, n_src, xd, r0, phi4, fd );
  }
  else{
     w = rbf_weight (KDDim, n_src, xd, r0, phi4, fd );
  }
 
  /*
    xi: points where I want to know the values?
  */
  Real* xi= new Real[KDDim];//{0,0,0};  // interpolation point is '0' always
  for(unsigned int i=0; i<KDDim; i++){
     xi[i]=0;
  }
  unsigned int ni = 1;
  Real* fi;
  //fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );
  if (minnorm<inner_range){
     fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi1, w, ni, xi );
     out<<fi[0]<<"  ";
     fi[0]-=(Real)_geom[closest].id()/(pt-_geom[closest]).norm(); //correct for the potential that was removed before.
     //out<<"dist:     "<< (pt-_geom[closest]).norm()<<std::endl;
     out<<fi[0]<<std::endl;
  }
  else if (minnorm > 8.)
     fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi4, w, ni, xi );
  else 
     fi = rbf_interp_nd ( KDDim, n_src, xd, r0, phi4, w, ni, xi );
  delete xi;
  // don't forget setting the output buffer!

  for (unsigned int v=0; v<n_fv; v++, ++out_it)
    {
    _vals[v] = fi[v];
    *out_it = _vals[v];
  }
}

// ------------------------------------------------------------
// Explicit Instantiations
template class RBFInterpolation<1>;
template class RBFInterpolation<2>;
template class RBFInterpolation<3>;

} // namespace libMesh
