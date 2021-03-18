#ifndef  _SNCORECONFIGURATIONBC_H_
#define  _SNCORECONFIGURATIONBC_H_


#include "GridFunction.H"
#include "FluidVarBC.H"
#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

/**
 * Boundary condition interface class for fluid quantities in single null geometries
 *
*/
class SNCoreConfigurationBC
   : public FluidVarBC
{
   public:

      // Constructor
      SNCoreConfigurationBC( const std::string&  species_name,
                             const std::string&  variable_name,
                             const int&          verbosity = 0 );

      /// Destructor.
      /**
       */
      virtual ~SNCoreConfigurationBC();

   
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
      SNCoreConfigurationBC( const SNCoreConfigurationBC& );
      SNCoreConfigurationBC& operator=( const SNCoreConfigurationBC& );


      int m_verbosity;

      enum {INVALID=-1,
            RADIAL_INNER,
            RADIAL_OUTER,
            NUM_BOUNDARIES};

};

#include "NamespaceFooter.H"

#endif

