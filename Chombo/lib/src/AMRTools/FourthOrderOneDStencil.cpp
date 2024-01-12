#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

#include "FourthOrderOneDStencil.H"
#include "BoxIterator.H"


#include "NamespaceHeader.H"
//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderOneDStencil::FourthOrderOneDStencil()
{
  m_defined = false;
}




void FourthOrderOneDStencil::define(const int& a_refRatio,
                                    const int& a_order   )
{
  m_refRatio = a_refRatio;
  m_order    = a_order;
  int numStencilWeights = m_order + 1;
 
  m_weights = new Real* [m_refRatio];
  for (int iWeight = 0; iWeight <numStencilWeights ; ++ iWeight)
    {
      m_weights[iWeight] = new Real[5];
    }
  
  //first row
  m_weights[0][0] = -  35.0/1024;
  m_weights[0][1] =   294.0/1024;
  m_weights[0][2] =   960.0/1024;
  m_weights[0][3] =  -230.0/1024;
  m_weights[0][4] =    35.0/1024;

  //second row
  m_weights[1][0] = -  13.0/1024;
  m_weights[1][1] =    58.0/1024;
  m_weights[1][2] =  1088.0/1024;
  m_weights[1][3] = - 122.0/1024;
  m_weights[1][4] =    13.0/1024;

  //third row
  m_weights[2][0] =    13.0/1024;
  m_weights[2][1] = - 122.0/1024;
  m_weights[2][2] =  1088.0/1024;
  m_weights[2][3] =    58.0/1024;
  m_weights[2][4] = -  13.0/1024;

  // fourth row
  m_weights[3][0] =    35.0/1024;
  m_weights[3][1] = - 230.0/1024;
  m_weights[3][2] =   960.0/1024;
  m_weights[3][3] =   294.0/1024;
  m_weights[3][4] = -  35.0/1024;
}

void FourthOrderOneDStencil::calculateStencilWeights(LAPACKMatrix       & a_weights  ,
                                                             const LAPACKMatrix & a_matrix   ,
                                                             const LAPACKMatrix & a_coarseData)
{
  LAPACKMatrix inverseMatrix(a_matrix);
  inverseMatrix.invert();
  multiply(a_weights,inverseMatrix,a_coarseData);
}
  
void FourthOrderOneDStencil::makeMomentMatrix(LAPACKMatrix & a_matrix,
                                                        const Real   & a_dx    ,
                                                        const Real   & a_vol   )
{
  //calculate the integral over a given control volume of m(x) = (x[idir] - a_offset[idir])^a_moment[idir]
  //1/dx*dy int (x-xBar)^mi  dx  =
  // (1/m0+1)(1/m1+1)(1/m2+1)(x-xBar)(m0+1) (y-yBar)(m1+1) (z-zBar)(m1+1)
    
  Vector <Real> momentMap(9);
  calculateMoments (momentMap,a_dx,a_vol);
    
  for (int irow = 0; irow < 3; ++ irow)
    for (int icol = 0; icol < 3; ++icol)
      {
        a_matrix(irow,icol) = momentMap[irow + icol];
      }
}
    
    
void FourthOrderOneDStencil::calculateMoments(Vector <Real> & a_oneDMomentMap,
                                                      const Real    & a_dx           ,
                                                      const Real    & a_vol          )
{
  for (int hilo = -1; hilo <= 1; hilo += 2)
    {
      for (int ioffset = 0; ioffset <=2; ++ ioffset)
        {
          // ioffset and hilo together indicate which one of the five cells were using: two to left, one to the left, center cell, one to the right, two to the right
          Real upperLimit = hilo*(ioffset + 0.5) * a_dx;
          //Real lowerLimit = hilo*(ioffset - 0.5) * a_dx;
          for (int imom = 0; imom <= 3; ++ imom)
            {
              //anti-derivative of (x-a_xBar)^imom = (x-xBar)^(imom + 1)/(imom + 1)
              Real moment = pow(upperLimit,imom + 1)/(imom + 1);
              a_oneDMomentMap.push_back(moment);
            }
        }
    }
}

 void FourthOrderOneDStencil::apply(FArrayBox      & a_fineFab    ,
                                    const FArrayBox& a_coarseFab  ,
                                    const IntVect  & a_coarseIv   ,
                                    const int      & a_interpDir  ,
                                    const IntVect  & a_changeVecDx)
 {
   // each row comprises the stencil weights for one fine cell out of 4
   
   int ncomp = a_fineFab.nComp();

   Box fineBox(a_coarseIv,a_coarseIv);
   fineBox.refine(a_changeVecDx);
  
   
   for (int icomp = 0; icomp<ncomp; ++icomp)
     {
       for (int iFine = 0;iFine < m_refRatio; ++iFine)
         {
           IntVect fineIv = fineBox.smallEnd();
           fineIv += iFine*BASISV(a_interpDir);
           // pout()<<"fineIv = "<<fineIv<<endl;
           
           if (a_fineFab.box().contains(fineIv) == false)
             {
               MayDay::Abort("a_fineFab.box().contains(fineIv) = false");
             }

           int iWeightIndex;
           int jWeightIndex;
           
           iWeightIndex = iFine;
           for (int hiLo = -2; hiLo <=2; ++hiLo)
             {
               jWeightIndex = hiLo + 2;
               Real weight =  m_weights[iWeightIndex][jWeightIndex];
               IntVect offset =  hiLo*BASISV(a_interpDir);

               if (a_coarseFab(a_coarseIv + offset,icomp) == 0.0)
                 {
                   cout<<"this is the stop point"      <<endl;

                   cout <<"offset             = "  <<offset             <<endl;
                   cout<<"a_coarseIv + offset = "  <<a_coarseIv + offset<<endl;
                   cout<<"iFine               = "  <<iFine              <<endl;
                   cout<<"hiLo                = "  <<hiLo               <<endl;
                   cout<<"weight              = "  <<weight             <<endl;
                   cout<<"fineIv              = "  <<fineIv             <<endl;
                   cout<<"a_coarseIv          = "  <<a_coarseIv         <<endl;
                   cout<<"a_coarseFab(coarseIv + offset,icomp) = " <<a_coarseFab(a_coarseIv + offset,icomp)<<endl;
                   cout<<"Before:a_fineFab(fineIv,icomp) = " <<a_fineFab(fineIv,icomp)<<endl;
                  
                   MayDay::Abort("bad Ghost values");
                 }
               a_fineFab(fineIv,icomp) += weight*a_coarseFab(a_coarseIv + offset,icomp);
               
#if 0
               cout<<"this is the stop point"      <<endl;

               cout <<"offset             = "  <<offset             <<endl;
               cout<<"a_coarseIv + offset = "  <<a_coarseIv + offset<<endl;
               cout<<"iFine               = "  <<iFine              <<endl;
               cout<<"hiLo                = "  <<hiLo               <<endl;
               cout<<"weight              = "  <<weight             <<endl;
               cout<<"fineIv              = "  <<fineIv             <<endl;
               cout<<"a_coarseIv          = "  <<a_coarseIv         <<endl;
               cout<<"a_coarseFab(coarseIv + offset,icomp) = " <<a_coarseFab(a_coarseIv + offset,icomp)<<endl;
               cout<<"Before:a_fineFab(fineIv,icomp) = " <<a_fineFab(fineIv,icomp)<<endl;
               cout<<"After:a_fineFab(fineIv,icomp)  = " <<a_fineFab(fineIv,icomp)<<endl;
#endif
             }
         }
     }
 }

#include "NamespaceFooter.H"
