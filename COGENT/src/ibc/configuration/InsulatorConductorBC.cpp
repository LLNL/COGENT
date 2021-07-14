#include "InsulatorConductorBC.H"
#include "TimeFunctionFactory.H"
#include "InsulatorConductorBCUtils.H"

#include "NamespaceHeader.H"

InsulatorConductorBC::InsulatorConductorBC( const string&          a_species_name,
                                            const string&          a_variable_name,
                                            const Vector<string>&  a_bdry_name )
: m_species_name(a_species_name),
  m_variable_name(a_variable_name)
{

   m_bdry_is_conductor.resize( a_bdry_name.size(), 0 );
   m_bdry_is_insulator_conductor.resize( a_bdry_name.size(), 0 );
   
   parseParameters(a_bdry_name);

}

void InsulatorConductorBC::defineInsulatorConductorBinary( const LevelData<FArrayBox>&      a_data,
                                                           const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                                           const MagGeom&                   a_geometry ) const
{
   CH_TIME("InsulatorConductorBC::defineInsulatorConductorBinary() for FArrayBox");
   if( !m_ICbinary.isDefined() || (m_ICbinary.ghostVect() < a_data.ghostVect()) ) {
 
      const DisjointBoxLayout& grids( a_geometry.grids() );
      m_ICbinary.define(grids, 1, a_data.ghostVect());
      for (DataIterator dit( grids ); dit.ok(); ++dit) {
         m_ICbinary[dit].setVal(0.0); // initialize to conductor everywhere
      }
 
      for(int b(0); b<a_all_bdry_layouts.size(); b++) {
         if(m_bdry_is_insulator_conductor[b]) {
            // WARNING!!! THIS NEEDS TO BE MODIFIED IF INSULATOR 
            //            EXISTS ON MORE THAN ONE BOUNDARY
            //            ALSO NEED TO SPECIFY WHICH BOUNDARIES ARE CONDUCTORS
            const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
            InsulatorConductorBCUtils::defineInsulatorConductorBinary( m_ICbinary,
                                                                       bdry_layout,
                                                                       m_Xmin_insulator,
                                                                       m_Xmax_insulator,
                                                                       a_geometry );
         }
      }
   }
}

void InsulatorConductorBC::defineInsulatorConductorBinary( const LevelData<EdgeDataBox>&    a_data,
                                                           const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                                           const MagGeom&                   a_geometry ) const
{
   CH_TIME("InsulatorConductorBC::defineInsulatorConductorBinary() for EdgeDataBox");
   if( !m_ICbinary_ce.isDefined() || (m_ICbinary_ce.ghostVect() < a_data.ghostVect()) ) {
 
      const DisjointBoxLayout& grids( a_geometry.grids() );
      m_ICbinary_ce.define(grids, 1, a_data.ghostVect());
      for (DataIterator dit( grids ); dit.ok(); ++dit) {
         for (int dir=0; dir<SpaceDim; dir++) {
            m_ICbinary_ce[dit].setVal(0.0); // initialize to conductor everywhere
         }
      }
 
      for(int b(0); b<a_all_bdry_layouts.size(); b++) {
         if(m_bdry_is_insulator_conductor[b]) {
            // WARNING!!! THIS NEEDS TO BE MODIFIED IF INSULATOR 
            //            EXISTS ON MORE THAN ONE BOUNDARY
            //            ALSO NEED TO SPECIFY WHICH BOUNDARIES ARE CONDUCTORS
            const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
            InsulatorConductorBCUtils::defineInsulatorConductorBinary( m_ICbinary_ce,
                                                                       bdry_layout,
                                                                       m_Xmin_insulator,
                                                                       m_Xmax_insulator,
                                                                       a_geometry );
         }
      }
   }
}

void InsulatorConductorBC::applyBC( LevelData<FArrayBox>&      a_data,
                              const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                              const Vector<string>&            a_all_bc_type,
                              const MagGeom&                   a_geometry,   
                              const Real                       a_time ) const
{
   CH_TIME("InsulatorConductorBC::applyBC()");
      
   defineInsulatorConductorBinary( a_data, a_all_bdry_layouts, a_geometry );
   
   for(int b(0); b<a_all_bdry_layouts.size(); b++) {
      if(m_bdry_is_insulator_conductor[b] || m_bdry_is_conductor[b]) {
      
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         CH_assert(this_bc_type=="insulatorConductor_Bv"); 
         Real this_ft = 0.0;
         if(m_bdry_is_insulator_conductor[b]) m_timeFunction->getValue(this_ft,a_time);
         InsulatorConductorBCUtils::setInsulatorConductorBC( a_data,
                                                             m_ICbinary,
                                                             bdry_layout,
                                                             this_bc_type,
                                                             this_ft,
                                                             m_Xmin_insulator,
                                                             a_geometry );
      }
   }

}

void InsulatorConductorBC::setInsulatorBC( LevelData<FArrayBox>&      a_dst,
                                     const LevelData<FArrayBox>&      a_src,
                                     const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                     const Vector<string>&            a_all_bc_type,
                                     const MagGeom&                   a_geometry,   
                                     const Real                       a_time ) const
{
   CH_TIME("InsulatorConductorBC::applyBC()");
      
   defineInsulatorConductorBinary( a_dst, a_all_bdry_layouts, a_geometry );
   
   for(int b(0); b<a_all_bdry_layouts.size(); b++) {
      if(m_bdry_is_insulator_conductor[b]) {
         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         InsulatorConductorBCUtils::setInsulatorBC( a_dst,
                                                    a_src,
                                                    m_ICbinary,
                                                    bdry_layout );
      }
   }

}

void InsulatorConductorBC::applyEdgeBC( LevelData<EdgeDataBox>&    a_data,
                                  const BoundaryBoxLayoutPtrVect&  a_all_bdry_layouts,
                                  const Vector<string>&            a_all_bc_type,
                                  const MagGeom&                   a_geometry ) const
{
   CH_TIME("InsulatorConductorBC::applyEdgeBC()");
   
   defineInsulatorConductorBinary( a_data, a_all_bdry_layouts, a_geometry );

   // Loop over boundaries and find those that contain insulator
   for(int b(0); b<a_all_bdry_layouts.size(); b++) {
      
      if(m_bdry_is_insulator_conductor[b]) {

         const BoundaryBoxLayout& bdry_layout( *(a_all_bdry_layouts[b]) );
         const std::string this_bc_type (a_all_bc_type[b]);
         if(this_bc_type=="insulatorConductor_currentDensity") {
            InsulatorConductorBCUtils::setInsulatorConductorBC( a_data,
                                                                m_ICbinary_ce,
                                                                bdry_layout,
                                                                this_bc_type );
         }
      }

   }
}

void
InsulatorConductorBC::parseParameters( const Vector<string>&  a_bdry_name )
{
   for (int i(0); i<a_bdry_name.size(); i++) {
      string prefix0 = "BC.";
      prefix0 += a_bdry_name[i];
      ParmParse fpp0( prefix0.c_str() );
      fpp0.query("is_conductor", m_bdry_is_conductor[i]);
      fpp0.query("is_insulator_conductor", m_bdry_is_insulator_conductor[i]);
      if( m_bdry_is_conductor[i] ) {
         if(!procID()) cout << "- species: "   << m_species_name << endl;
         if(!procID()) cout << "- variable: "  << m_variable_name << endl;
         if(!procID()) cout << "  " << a_bdry_name[i] << " is conductor " << endl;
      }
      if( m_bdry_is_insulator_conductor[i] ) {
        
         fpp0.get( "rmin_insulator", m_Xmin_insulator[0] );
         fpp0.get( "rmax_insulator", m_Xmax_insulator[0] );
         fpp0.get( "zmin_insulator", m_Xmin_insulator[1] );
         fpp0.get( "zmax_insulator", m_Xmax_insulator[1] );
        
         // create time-function for insulator BC
         std::string prefixtf( fpp0.prefix() );
         prefixtf += ".time_function";
         ParmParse tfpp( prefixtf.c_str() );
         TimeFunctionFactory  timeFactory;
         m_timeFunction = timeFactory.create(tfpp,1);
         
         if(!procID()) cout << "- species: "   << m_species_name << endl;
         if(!procID()) cout << "- variable: "  << m_variable_name << endl;
         if(!procID()) cout << "- knows about insulator on boundary " << endl;
         if(!procID()) cout << "  " << a_bdry_name[i] << ": " << endl;
         if(!procID()) cout << "m_rmin_insulator = " << m_Xmin_insulator[0] << endl;
         if(!procID()) cout << "m_rmax_insulator = " << m_Xmax_insulator[0] << endl;
         if(!procID()) cout << "m_zmin_insulator = " << m_Xmin_insulator[1] << endl;
         if(!procID()) cout << "m_zmax_insulator = " << m_Xmax_insulator[1] << endl;
      }
   }
   
   // make sure at least one boundary defined
   // as conductor or insulator_conductor (but not as both)
   int num_conductor_bdrys = 0;
   int num_insulator_conductor_bdrys = 0;
   int num_overlap_bdrys = 0;
   for (int i(0); i<m_bdry_is_conductor.size(); i++) {
      num_conductor_bdrys += m_bdry_is_conductor[i];
      num_insulator_conductor_bdrys += m_bdry_is_insulator_conductor[i];
      num_overlap_bdrys += m_bdry_is_conductor[i]*m_bdry_is_insulator_conductor[i];
   }
   CH_assert(num_conductor_bdrys > 0 || num_insulator_conductor_bdrys > 0);
   CH_assert(num_overlap_bdrys == 0);

}

#include "NamespaceFooter.H"

