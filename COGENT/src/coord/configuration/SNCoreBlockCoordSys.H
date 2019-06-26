#ifndef _SNCOREBLOCKCOORDSYS_H_
#define _SNCOREBLOCKCOORDSYS_H_

#include "SingleNullBlockCoordSys.H"

#include "NamespaceHeader.H"
namespace POL = POL_NAMESPACE;

class SNCoreBlockCoordSys
   : public SingleNullBlockCoordSys
{
   public:

      /// Constructor
      /**
       * Constructor with initialization.  Input data is read from
       * the parm_parse object, the cellSpacing is computed and
       * the define() member is called.
       *
       * @param[in] pp               the ParmParse database.
       * @param[in] domain           the configuration space domain.
       * @param[in] dx               mapped grid cell size
       * @param[in] poloidal_index   poloidal block index
       * @param[in] toroidal_index   toroidal block index
       */
      SNCoreBlockCoordSys( ParmParse&            pp,
                           const ProblemDomain&  domain,
                           const RealVect&       dx,
                           const int             poloidal_index,
                           const int             toroidal_index)
         : SingleNullBlockCoordSys(pp, domain, dx, poloidal_index, toroidal_index) {};

      /// Destructor
      /**
       */
      virtual ~SNCoreBlockCoordSys() {};

      virtual void printInit() const;

      virtual bool blockNameIsValid( const string& block_name ) const;

      virtual void definePoints( const ParmParse&  pp,
                                 int&              block_poloidal,
                                 const int         block_full_poloidal,
                                 const IntVect&    mapping_block_size,
                                 const int&        n_poloidal_extend,
                                 double &          dtheta,
                                 double *          theta_pts ) const;

      bool isValid(const RealVect& xi, const bool omit_toroidal) const;

      virtual string geometryType() const {return "SNCore";}
   
      static const std::string pp_name;

      enum SNCoreBlockType {LCORE,RCORE,MCORE,NUM_SNCORE_BLOCKS};

};

#include "NamespaceFooter.H"

#endif
