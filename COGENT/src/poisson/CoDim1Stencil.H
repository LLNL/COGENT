#ifndef _CODIM1STENCIL_
#define _CODIM1STENCIL_

#include "FArrayBox.H"
#include "RealVect.H"

#include "NamespaceHeader.H"



class CoDim1Stencil
{
public:

   CoDim1Stencil()
      : m_boundary_values(NULL),
        m_defined(false) {};
   
   CoDim1Stencil( const CoDim1Stencil& a_bc_stencil_coupling );

   virtual ~CoDim1Stencil();
   
   virtual void define( const FArrayBox&      boundary_values,
                        const int&            value_type,
                        const RealVect&       dx,
                        const int&            dir,
                        const Side::LoHiSide& side,
                        const int             order );

   virtual void getStencil( const IntVect&   iv,
                            vector<IntVect>& points,
                            vector<double>&  weights,
                            double&          bv_contrib ) const;

   inline Box box() const {return m_box;}
   inline int dir() const {return m_dir;}
   inline int side() const {return m_side;}

private:
   
   FArrayBox*         m_boundary_values;
   int                m_value_type;
   int                m_dir;
   Side::LoHiSide     m_side;
   int                m_order;
   Box                m_box;
   RealVect           m_dx;
   bool               m_defined;
};


#include "NamespaceFooter.H"

#endif


