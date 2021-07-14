#ifndef _CODIM2STENCIL_H_
#define _CODIM2STENCIL_H_

#include "Box.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

class CoDim2Stencil
{
public:

   CoDim2Stencil()
      : m_defined(false) {};
   
   CoDim2Stencil( const CoDim2Stencil& a_bc_stencil_coupling );

   virtual ~CoDim2Stencil() {};
   
   virtual void define( const Box&            box,
                        const int             value_type,
                        const RealVect&       dx,
                        const int&            dir,
                        const Side::LoHiSide& side,
                        const int&            transverse_dir,
                        const Side::LoHiSide& transverse_side,
                        const int             order,
                        const bool            transverse_boundary );

   virtual void getStencil( const IntVect&   iv,
                            vector<IntVect>& points,
                            vector<double>&  weights ) const;

   inline Box box() const {return m_box;}
   inline int dir() const {return m_dir;}
   inline int transverse_dir() const {return m_transverse_dir;}
   inline int side() const {return m_side;}
   inline int transverse_side() const {return m_transverse_side;}
   inline bool transverse_boundary() const {return m_transverse_boundary;}

private:
   
   int                m_dir;
   int                m_transverse_dir;
   Side::LoHiSide     m_side;
   Side::LoHiSide     m_transverse_side;
   int                m_order;
   bool               m_transverse_boundary;
   Box                m_box;
   RealVect           m_dx;
   int                m_value_type;
   bool               m_defined;
};


#include "NamespaceFooter.H"

#endif


