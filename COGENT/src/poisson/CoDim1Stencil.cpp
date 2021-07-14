#include "CoDim1Stencil.H"
#include "EllipticOpBC.H"

#include "NamespaceHeader.H"



CoDim1Stencil::CoDim1Stencil( const CoDim1Stencil& a_bc_stencil_coupling )
   : m_value_type(a_bc_stencil_coupling.m_value_type),
     m_dir(a_bc_stencil_coupling.m_dir),
     m_side(a_bc_stencil_coupling.m_side),
     m_order(a_bc_stencil_coupling.m_order),
     m_box(a_bc_stencil_coupling.m_box),
     m_dx(a_bc_stencil_coupling.m_dx)
{
   if (a_bc_stencil_coupling.m_boundary_values != NULL) {
      m_boundary_values = new FArrayBox(a_bc_stencil_coupling.m_boundary_values->box(),
                                        a_bc_stencil_coupling.m_boundary_values->nComp());
      m_boundary_values->copy(*a_bc_stencil_coupling.m_boundary_values);
   }
   else {
      m_boundary_values = NULL;
   }

   m_defined = true;
}



CoDim1Stencil::~CoDim1Stencil()
{
   if (m_boundary_values) delete m_boundary_values;
}



void
CoDim1Stencil::define( const FArrayBox&      a_boundary_values,
                       const int&            a_value_type,
                       const RealVect&       a_dx,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const int             a_order )
{
   CH_assert(!m_defined);
   m_boundary_values = new FArrayBox(a_boundary_values.box(), a_boundary_values.nComp());
   m_boundary_values->copy(a_boundary_values);
   m_dx = a_dx;
   m_value_type = a_value_type;
   m_dir = a_dir;
   m_side = a_side;
   m_order = a_order;

   double width;

   switch( m_order )
      {
      case 2:
         width = 1;
         break;
      case 4:
         width = 2;
         break;
      }

   if (a_side == Side::LoHiSide::Lo) {
      m_box = m_boundary_values->box();
      m_box.shiftHalf(m_dir,-1);
      m_box.growLo(m_dir,width-1);

      // Neumann boundary conditions are specified relative to the outward
      // normal, so we need to flip the sign on low side boundaries to
      // convert to non-directional derivatives.
      if (a_value_type == EllipticOpBC::NEUMANN) {
         m_boundary_values->negate();
      }
   }
   else {
      m_box = m_boundary_values->box();
      m_box.shiftHalf(m_dir,1);
      m_box.growHi(m_dir,width-1);
   }

   m_defined = true;
}



void
CoDim1Stencil::getStencil( const IntVect&   a_iv,
                           vector<IntVect>& a_points,
                           vector<double>&  a_weights,
                           double&          a_bv_contrib ) const
{
   a_points.resize(0);
   a_weights.resize(0);
   a_bv_contrib = 0.;

   if ( m_box.contains(a_iv) ) {

      if ( m_value_type == EllipticOpBC::DIRICHLET ) {

         if ( m_order == 2 ) {

            int num_points = 1;
            double weights[1] = {-1.};

            IntVect point = a_iv;
            if (m_side == Side::LoHiSide::Lo) {
               for (int i=0; i<num_points; ++i) {
                  point[m_dir]++;
                  a_points.push_back(point);
                  a_weights.push_back(weights[i]);
               }

               IntVect bv_point = a_iv;
               bv_point[m_dir]++;

               a_bv_contrib = 2. * (*m_boundary_values)(bv_point,0);
            }
            else {
               for (int i=0; i<num_points; ++i) {
                  point[m_dir]--;
                  a_points.push_back(point);
                  a_weights.push_back(weights[i]);
               }

               IntVect bv_point = a_iv;

               a_bv_contrib = 2. * (*m_boundary_values)(bv_point,0);
            }
         }
         else if ( m_order == 4 ) {

            int num_ghosts = 2;
            int num_points = 3;
            double weights[2][3] = {{-3., 1., -1./5.}, {-18., 8., -9./5.}};
            double bv_weight[2] = {16./5., 64./5.};

            IntVect point = a_iv;
            if (m_side == Side::LoHiSide::Lo) {
               int dist = m_box.bigEnd(m_dir) - a_iv[m_dir];
               CH_assert(dist >= 0 && dist < num_ghosts);

               for (int i=0; i<dist; ++i) {
                  point[m_dir]++;
               }

               for (int i=0; i<num_points; ++i) {
                  point[m_dir]++;
                  a_points.push_back(point);
                  a_weights.push_back(weights[dist][i]);
               }

               IntVect bv_point = a_iv;
               bv_point[m_dir] += (dist + 1);

               a_bv_contrib = bv_weight[dist] * (*m_boundary_values)(bv_point,0);
            }
            else {
               int dist = a_iv[m_dir] - m_box.smallEnd(m_dir);
               CH_assert(dist >= 0 && dist < num_ghosts);

               for (int i=0; i<dist; ++i) {
                  point[m_dir]--;
               }

               for (int i=0; i<num_points; ++i) {
                  point[m_dir]--;
                  a_points.push_back(point);
                  a_weights.push_back(weights[dist][i]);
               }

               IntVect bv_point = a_iv;
               bv_point[m_dir] -= dist;

               a_bv_contrib = bv_weight[dist] * (*m_boundary_values)(bv_point,0);
            }

         }
         else {
            MayDay::Error("CoDim1Stencil::getStencil(): unknown order");
         }
      }
      else if (m_value_type == EllipticOpBC::NEUMANN) {

         if ( m_order == 2 ) {

            int num_points = 1;
            double weights[] = {1.};
            double bv_weight = -1.;

            IntVect point = a_iv;
            if (m_side == Side::LoHiSide::Lo) {
               for (int i=0; i<num_points; ++i) {
                  point[m_dir]++;
                  a_points.push_back(point);
                  a_weights.push_back(weights[i]);
               }

               IntVect bv_point = a_iv;
               bv_point[m_dir]++;

               a_bv_contrib = m_dx[m_dir] * bv_weight * (*m_boundary_values)(bv_point,0);
            }
            else {
               for (int i=0; i<num_points; ++i) {
                  point[m_dir]--;
                  a_points.push_back(point);
                  a_weights.push_back(weights[i]);
               }

               IntVect bv_point = a_iv;

               a_bv_contrib = -m_dx[m_dir] * bv_weight * (*m_boundary_values)(bv_point,0);
            }
         }
         else if (m_order == 4) {

            int num_points = 3;
            double weights[2][3] = {{21./23., 3./23., -1./23.}, {-54./23., 104./23., -27./23.}};
            double bv_weights[2] = {-24./23., -96./23.};

            IntVect point = a_iv;
            if (m_side == Side::LoHiSide::Lo) {
               int dist = m_box.bigEnd(m_dir) - a_iv[m_dir];
               CH_assert(dist >= 0 && dist < 2);

               for (int i=0; i<dist; ++i) {
                  point[m_dir]++;
               }

               for (int i=0; i<num_points; ++i) {
                  point[m_dir]++;
                  a_points.push_back(point);
                  a_weights.push_back(weights[dist][i]);
               }

               IntVect bv_point = a_iv;
               bv_point[m_dir] += (dist + 1);

               a_bv_contrib = m_dx[m_dir] * bv_weights[dist] * (*m_boundary_values)(bv_point,0);
            }
            else {
               int dist = a_iv[m_dir] - m_box.smallEnd(m_dir);
               CH_assert(dist >= 0 && dist < 2);

               for (int i=0; i<dist; ++i) {
                  point[m_dir]--;
               }

               for (int i=0; i<num_points; ++i) {
                  point[m_dir]--;
                  a_points.push_back(point);
                  a_weights.push_back(weights[dist][i]);
               }

               IntVect bv_point = a_iv;
               bv_point[m_dir] -= dist;

               a_bv_contrib = -m_dx[m_dir] * bv_weights[dist] * (*m_boundary_values)(bv_point,0);
            }

         }
         else {
            MayDay::Error("CoDim1Stencil::getStencil(): unknown order");
         }
      }
      else if (m_value_type == EllipticOpBC::EXTRAPOLATED) {
         
         if ( m_order == 2 ) {
            int num_points = 2;
            double weights[] = {2., -1.};
         
            IntVect point = a_iv;
            if (m_side == Side::LoHiSide::Lo) {
               for (int i=0; i<num_points; ++i) {
                  point[m_dir]++;
                  a_points.push_back(point);
                  a_weights.push_back(weights[i]);
               }
            }
            else {
               for (int i=0; i<num_points; ++i) {
                  point[m_dir]--;
                  a_points.push_back(point);
                  a_weights.push_back(weights[i]);
               }
            }
         }
         else if ( m_order == 4 ) {
            MayDay::Error("CoDim1Stencil::getStencil(): fourth-order extrapolation option not implemented");
         }
         else {
            MayDay::Error("CoDim1Stencil::getStencil(): unknown order");
         }
      }
      else {
         MayDay::Error("CoDim1Stencil::getStencil(): unknown value type");
      }
   }
}


#include "NamespaceFooter.H"
