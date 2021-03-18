#ifndef  _LOGRECTCONFIGURATIONBC_H_
#define  _LOGRECTCONFIGURATIONBC_H_

#include "FluidVarBC.H"
#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for fluid quantities in single null geometries
 *
*/
class LogRectConfigurationBC
   : public FluidVarBC
{
   public:

      // Constructor
      LogRectConfigurationBC( const std::string&  species_name,
                              const std::string&  variable_name,
                              const int&          verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~LogRectConfigurationBC();


   
   protected:
      
      virtual RefCountedPtr<GridFunction> getBCFunction(const int block_number,
                                                        const int dir,
                                                        const int side );

      virtual std::string getBCType(const int block_number,
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
      LogRectConfigurationBC( const LogRectConfigurationBC& );
      LogRectConfigurationBC& operator=( const LogRectConfigurationBC& );

      int m_verbosity;

  enum {INVALID=-1,
	RADIAL_LOWER, RADIAL_UPPER,
	POLOIDAL_LOWER, POLOIDAL_UPPER,
#if CFG_DIM==3
	TOROIDAL_LOWER, TOROIDAL_UPPER,
#endif
	NUM_BOUNDARIES};

};

#include "NamespaceFooter.H"

#endif

