#include "inspect.H"
#include "inspectF_F.H"

#include "NamespaceHeader.H"



void
inspect(const FArrayBox& data)
{
   FORT_INSPECT(CHF_BOX(data.box()),
                CHF_CONST_FRA(data));
}



void
inspect(const FluxBox& data)
{
   for (int dir=0; dir<SpaceDim; ++dir) {
      const FArrayBox& this_fab = data[dir];
      FORT_INSPECT(CHF_BOX(this_fab.box()),
                   CHF_CONST_FRA(this_fab));
   }
}



void
inspect(const LevelData<FArrayBox>& data)
{
   DataIterator dit = data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      const FArrayBox& this_fab = data[dit];
      FORT_INSPECT(CHF_BOX(this_fab.box()),
                   CHF_CONST_FRA(this_fab));
   }
}



void
inspect(const LevelData<FluxBox>& data)
{
   DataIterator dit = data.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      for (int dir=0; dir<SpaceDim; ++dir) {
         const FArrayBox& this_fab = data[dit][dir];
         FORT_INSPECT(CHF_BOX(this_fab.box()),
                      CHF_CONST_FRA(this_fab));
      }
   }
}



#include "NamespaceFooter.H"
