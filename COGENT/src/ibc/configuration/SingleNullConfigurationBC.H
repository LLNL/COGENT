#ifndef  _SINGLENULLCONFIGURATIONBC_H_
#define  _SINGLENULLCONFIGURATIONBC_H_

#include "FluidVarBC.H"

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for fluid quantities in single null geometries
 *
*/
class SingleNullConfigurationBC
   : public FluidVarBC
{
   public:

      // Constructor
      SingleNullConfigurationBC( const std::string& species_name,
                                 const std::string& variable_name,
                                 const int&         poloidal_blocks_per_sector,
                                 const int&         verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~SingleNullConfigurationBC();

   protected:
   
      virtual std::string getBCType(const int block_number,
                                    const int dir,
                                    const int side );
      
      virtual RefCountedPtr<GridFunction> getBCFunction(const int block_number,
                                                        const int dir,
                                                        const int side );
   
      virtual void setNames();
   
      
   private:
   
      /// Print object parameters.
      /**
       */
      void printParameters() const;


      /// Parse the input database for parameters.
      /**
       */
      void parseParameters(ParmParse& pp);

      // prohibit copying
      SingleNullConfigurationBC( const SingleNullConfigurationBC& );
      SingleNullConfigurationBC& operator=( const SingleNullConfigurationBC& );
 

      int m_verbosity;
      int m_poloidal_blocks_per_sector;

   enum {INVALID=-1,
         RADIAL_CORE, RADIAL_SOL, RADIAL_PF,
         POLOIDAL_INNER_DIV, POLOIDAL_OUTER_DIV,
#if CFG_DIM==3
         TOROIDAL_CORE, TOROIDAL_SOL, TOROIDAL_PF,
         TOROIDAL_INNER_DIV, TOROIDAL_OUTER_DIV,
#endif
         NUM_BOUNDARIES};
  
};

#include "NamespaceFooter.H"

#endif

