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

#ifndef NEARES_NEIGHBOURS_INTERPOLATION_H
#define NEARES_NEIGHBOURS_INTERPOLATION_H

// includes from libmesh
#include "libmesh/point.h"
#include "libmesh/node.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/parallel_object.h"
#ifdef LIBMESH_HAVE_NANOFLANN
#  include "libmesh/nanoflann.hpp"
#endif

// C++ includes
#include <string>
#include <vector>

namespace libMesh
{

/**
 * Radias basis interplation.
 */
template <unsigned int KDDim>
class NeNeInterpolation : public MeshfreeInterpolation
{
protected:

   #ifdef LIBMESH_HAVE_NANOFLANN
   /**
   * This class adapts list of libMesh \p Point types
   * for use in a nanoflann KD-Tree. For more on the
   * basic idea see examples/pointcloud_adaptor_example.cpp
   * in the nanoflann source tree.
   */
   template <unsigned int PLDim>
   class PointListAdaptor
   {
      private:
         const std::vector<Point> & _pts;
   
      public:
         PointListAdaptor (const std::vector<Point> & pts) :
                  _pts(pts){}
   
         /**
         * libMesh \p Point coordinate type
         */
         typedef Real coord_t;
   
         /**
         * Must return the number of data points
         */
         inline size_t kdtree_get_point_count() const { return _pts.size(); }
      
         /**
         * Returns the distance between the vector "p1[0:size-1]"
         * and the data point with index "idx_p2" stored in the class
         */
         inline coord_t kdtree_distance(const coord_t * p1, const size_t idx_p2, size_t size) const
         {
            libmesh_assert_equal_to (size, PLDim);
            libmesh_assert_less (idx_p2, _pts.size());
   
            const Point & p2(_pts[idx_p2]);
   
            switch (size)
            {
            case 3:
               {
                  const coord_t d0=p1[0] - p2(0);
                  const coord_t d1=p1[1] - p2(1);
                  const coord_t d2=p1[2] - p2(2);
      
                  return d0*d0 + d1*d1 + d2*d2;
               }
      
            case 2:
               {
                  const coord_t d0=p1[0] - p2(0);
                  const coord_t d1=p1[1] - p2(1);
      
                  return d0*d0 + d1*d1;
               }
      
            case 1:
               {
                  const coord_t d0=p1[0] - p2(0);
      
                  return d0*d0;
               }
      
            default:
               libmesh_error_msg("ERROR: unknown size " << size);
            }
   
         return -1.;
      }

    /**
     * Returns the dim'th component of the idx'th point in the class:
     * Since this is inlined and the "dim" argument is typically an immediate value, the
     *  "if's" are actually solved at compile time.
     */
    inline coord_t kdtree_get_pt(const size_t idx, int dim) const
    {
      libmesh_assert_less (dim, (int) PLDim);
      libmesh_assert_less (idx, _pts.size());
      libmesh_assert_less (dim, 3);

      const Point & p(_pts[idx]);

      if (dim==0) return p(0);
      if (dim==1) return p(1);
      return p(2);
    }

    /**
     * Optional bounding-box computation: return false to default to a standard bbox computation loop.
     * Return true if the BBOX was already computed by the class and returned in "bb" so it can be
     * avoided to redo it again. Look at bb.size() to find out the expected dimensionality
     * (e.g. 2 or 3 for point clouds)
     */
    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }
   };
   
   PointListAdaptor<KDDim> _point_list_adaptor;
   
   // template <int KDDIM>
   // class KDTree : public KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<num_t, PointListAdaptor >,
   //  PointListAdaptor,
   //  KDDIM>
   // {
   // };
      
   typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<Real, PointListAdaptor<KDDim> >,
                                                   PointListAdaptor<KDDim>, KDDim> kd_tree_t;
      
   mutable UniquePtr<kd_tree_t> _kd_tree;
   std::vector<Node> _geom;
   
   #endif // LIBMESH_HAVE_NANOFLANN
   
   /**
   * Build & initialize the KD tree, if needed.
   */
   virtual void construct_kd_tree ();
      
   /**
   * Performs inverse distance interpolation at the input point from
   * the specified points.
   */
   virtual void interpolate (const Point               & pt,
                            const std::vector<size_t> & src_indices,
                            const std::vector<Real>   & src_dist_sqr,
                            std::vector<Number>::iterator & out_it) const;
   
   const unsigned int _n_interp_pts;
   
   /**
   * Temporary work array.  Object level scope to avoid cache thrashing.
   */
   mutable std::vector<Number> _vals;

public:
   
   /**
   * Constructor. Takes the inverse distance power,
   * which defaults to 2.
   */
   NeNeInterpolation (const libMesh::Parallel::Communicator & comm_in,
                                 const unsigned int n_interp_pts ,
                                 const Real  /*power    */           ,
                                 const std::vector<Node> geometry) :
   #if LIBMESH_HAVE_NANOFLANN
       MeshfreeInterpolation(comm_in),
      _point_list_adaptor(_src_pts),
      _geom(geometry),
      _n_interp_pts(n_interp_pts)
   #else 
       MeshfreeInterpolation(comm_in),
      _geom(geometry)
   #endif
   {}
   
   /**
   * Clears all internal data structures and restores to a
   * pristine state.
   */
   virtual void clear() libmesh_override;
   
   /**
   * Interpolate source data at target points.
   * Pure virtual, must be overriden in derived classes.
   */
   virtual void interpolate_field_data (const std::vector<std::string> & field_names,
                                       const std::vector<Point>  & tgt_pts,
                                       std::vector<Number> & tgt_vals) const libmesh_override;
   };
   
} // namespace libMesh


#endif // #define NEARES_NEIGHBOURS_INTERPOLATION_H
