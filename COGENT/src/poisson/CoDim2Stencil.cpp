#include "CoDim2Stencil.H"
#include "EllipticOpBC.H"

#include "NamespaceHeader.H"



CoDim2Stencil::CoDim2Stencil( const CoDim2Stencil& a_bc_stencil_coupling )
   : m_dir(a_bc_stencil_coupling.m_dir),
     m_transverse_dir(a_bc_stencil_coupling.m_transverse_dir),
     m_side(a_bc_stencil_coupling.m_side),
     m_transverse_side(a_bc_stencil_coupling.m_transverse_side),
     m_order(a_bc_stencil_coupling.m_order),
     m_transverse_boundary(a_bc_stencil_coupling.m_transverse_boundary),
     m_box(a_bc_stencil_coupling.m_box),
     m_dx(a_bc_stencil_coupling.m_dx),
     m_value_type(a_bc_stencil_coupling.m_value_type),
     m_defined(true)
{
}



void
CoDim2Stencil::define( const Box&            a_codim1_box,
                       const int             a_value_type,
                       const RealVect&       a_dx,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const int&            a_transverse_dir,
                       const Side::LoHiSide& a_transverse_side,
                       const int             a_order,
                       const bool            a_transverse_boundary )
{
   CH_assert(!m_defined);

   int nghosts;
   if (a_order == 4) {
      nghosts = 2;
   }
   else {
      nghosts = 1;
   }

   m_box = adjCellBox(a_codim1_box, a_transverse_dir, a_transverse_side, nghosts);
   m_value_type = a_value_type;
   m_dx = a_dx;
   m_dir = a_dir;
   m_side = a_side;
   m_transverse_dir = a_transverse_dir;
   m_transverse_side = a_transverse_side;
   m_order = a_order;
   m_transverse_boundary = a_transverse_boundary;
   m_defined = true;
}



void
CoDim2Stencil::getStencil( const IntVect&   a_iv,
                           vector<IntVect>& a_points,
                           vector<double>&  a_weights ) const
{
   CH_assert(m_box.contains(a_iv));

   a_points.resize(0);
   a_weights.resize(0);

   double fac = m_transverse_boundary? 0.5: 1.0;

   if (m_order == 2) {

      int num_weights = 2;
      double weights[] = {2.,-1.};

      IntVect point = a_iv;

      for (int i=0; i<num_weights; ++i) {

         if ( m_value_type == EllipticOpBC::EXTRAPOLATED ) {
            if (m_side == Side::LoHiSide::Lo) {
               point[m_dir]++;
            }
            else if (m_side == Side::LoHiSide::Hi) {
               point[m_dir]--;
            }
         }
         else {
            point[m_transverse_dir] += (1 - 2*m_transverse_side);
         }
         a_points.push_back(point);
         a_weights.push_back(weights[i]*fac);
      }
   }
   else if (m_order == 4) {

      int num_ghosts = 2;
      int num_points = 3;
      double weights[2][3] = {{3., -3., 1.}, {6., -8., 3.}};

      IntVect point = a_iv;
      if (m_transverse_side == Side::LoHiSide::Lo) {
         int dist = m_box.bigEnd(m_transverse_dir) - a_iv[m_transverse_dir];
         CH_assert(dist >= 0 && dist < num_ghosts);

         for (int i=0; i<dist; ++i) {
            point[m_transverse_dir]++;
         }

         for (int i=0; i<num_points; ++i) {
            point[m_transverse_dir]++;
            a_points.push_back(point);
            a_weights.push_back(weights[dist][i]*fac);
         }
      }
      else {
         int dist = a_iv[m_transverse_dir] - m_box.smallEnd(m_transverse_dir);
         CH_assert(dist >= 0 && dist < num_ghosts);
         
         for (int i=0; i<dist; ++i) {
            point[m_transverse_dir]--;
         }

         for (int i=0; i<num_points; ++i) {
            point[m_transverse_dir]--;
            a_points.push_back(point);
            a_weights.push_back(weights[dist][i]*fac);
         }
      }
   }
   else {
      MayDay::Error("CoDim2Stencil::getCoDim2Stencil(): Unimplemented order");
   }
}



#include "NamespaceFooter.H"
