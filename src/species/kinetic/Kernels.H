#ifndef _KERNELS_H_
#define _KERNELS_H_

#include "KineticSpecies.H"

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

class KineticSpecies;

/// Kernel interface.
/**
 * This is a function object that is used to compute the moment integrand.
 */
class Kernel
{
   public:
  
      /// virtual destructor to prevent compiler complaints (DFM 2/4/09)
      virtual ~Kernel() {;}

      /// Computes the kernel.
      /**
       * Computes the kernel of a moment calculation.  
       *
       * If kcomp is the number of kernel components, and the distribution function
       * has more than one component, the first kcomp components of the result
       * will contain the moment of the first distribution function component,
       * the next kcomp components will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first kcomp components contain the first
       * component of the distribution function, the next kcomp components
       * contain the second component of the distribution function, etc.
       *
       * @param[in]  result Cell-averaged distribution function
       * @param[out] result Cell-averaged integrand
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const = 0;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const = 0;

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual int nComponents() const = 0;
    
      ///Computes cell-centered psi and theta velocity projections

      virtual void computeVelCell(LevelData<FArrayBox>& velCell,
                                  const PhaseGeom& phase_geom,
                                  const LevelData<FluxBox>& field) const;
};

/// Density kernel.
/**
 */
class DensityKernel : public Kernel
{
   public:

      /// destructor to prevent compiler from complaining (DFM 2/4/09)
      virtual ~DensityKernel() {;}

      /// Computes the kernel.
      /**
       * Computes the integrand, which is just one times the distribution
       * function.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
      */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is one.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.0; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
       virtual int nComponents() const { return 1; }
};

/// Charge density kernel.
/**
 */
class ChargeDensityKernel : public Kernel
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~ChargeDensityKernel() {;}

      /// Computes the kernel.
      /**
       * Computes the integrand and returns the charge as the scale.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is the unit charge.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return kinetic_species.charge(); }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
       virtual int nComponents() const { return 1; }
};

/// Mass density kernel.
/**
 */
class MassDensityKernel : public Kernel
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~MassDensityKernel() {;}

      /// Computes the kernel.
      /**
       * Computes the integrand and returns the charge as the scale.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is the unit charge.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return kinetic_species.mass(); }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
       virtual int nComponents() const { return 1; }
};

/// Momentum kernel.
/**
 */
class MomentumDensityKernel : public Kernel
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~MomentumDensityKernel() {;}

      /// Computes the integrand for a momentum density moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the velocity and species mass and returns the species mass as the scale.
       *
       * The number of kernel components is VEL_DIM.  If the distribution function
       * has more than one component, the first VELDIM components of the result
       * will contain the moment of the first distribution function component,
       * the next VEL_DIM components will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first VEL_DIM components contain the first
       * component of the distribution function, the next VEL_DIM components
       * contain the second component of the distribution function, etc.
       *
       * @param[in] result Cell-averaged distribution function.
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is the species mass.
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return kinetic_species.mass(); }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return VEL_DIM; }
};

/// 4th Moment kernel (used for transport models).
/**
 */
class FourthMomentKernel : public Kernel
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~FourthMomentKernel() {;}

      /// Computes the integrand for a fourth moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the scalar velocity to the 4th powerte.
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[in] result Cell-averaged distribution function.
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }
};

/// Pressure kernel.
/**
 */
class PressureKernel : public Kernel
{
   public:

      PressureKernel(const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift);

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~PressureKernel() {;}

      /// Computes the integrand for a pressure moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the scalar velocity squared.
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }

    private:

       const CFG::LevelData<CFG::FArrayBox>& m_v_parallel_shift;

};

/// perp energy kernel.
/**
 */
class PerpEnergyKernel : public Kernel
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~PerpEnergyKernel() {;}

      /// Computes the integrand for the perp pressure moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the scalar perp velocity squared.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }
};

/// Parallel momenutm kernel.
/**
 */
class ParallelMomKernel : public Kernel
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~ParallelMomKernel() {;}

      /// Computes the integrand for a ParallelMom moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the parallel velocity component.
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[out] result Cell-averaged integrand.
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }
};

/// Particle flux kernel.
/**
 */
class ParticleFluxKernel : public Kernel
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       ParticleFluxKernel( const LevelData<FluxBox>& field );

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~ParticleFluxKernel() {;}

      /// Computes the integrand for a ParticleFlux moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the radial component of the drift velocity.
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }

    private:

        const LevelData<FluxBox>& m_field;
};

/// Heat flux kernel.
/**
 */
class HeatFluxKernel : public Kernel
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       HeatFluxKernel( const LevelData<FluxBox>& field,
                       const LevelData<FArrayBox>& phi);

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~HeatFluxKernel() {;}

      /// Computes the integrand for a HeatFlux moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the scalar velocity squared times the radial component of the drift velocity.
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[out] result Cell-averaged integrand.
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }

    private:

        const LevelData<FluxBox>& m_field;
        const LevelData<FArrayBox>& m_phi;

};

/// Perpendicular current density kernel.
/**
 */
class PerpCurrentDensityKernel : public Kernel
{
public:
   /// Constructor with initialization.
   /**
    *  Initialized with potentail, which is used to get phase-velocity
    */
   PerpCurrentDensityKernel( const LevelData<FluxBox>& field );
   
   /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
   virtual ~PerpCurrentDensityKernel() {;}
   
   /// Computes the integrand for a PerpCurrent moment calculation.
   /**
    * Computes the integrand by multiplying the input distribution function
    * by the perpendicular component of the drift velocity.
    *
    * The number of kernel components is CFG_DIM.  If the distribution function
    * has more than one component, the first component of the result
    * will contain the moment of the first distribution function component,
    * the next component will contain the moment of the second
    * distribution function component, etc.  Upon entry, it is assumed
    * that result has already been initialized with the distribution function
    * correspondingly, ie., the first component contains the first
    * component of the distribution function, the next component
    * contain the second component of the distribution function, etc.
    *
    * @param[out] result Cell-averaged integrand.
    * @param[in] kinetic_species Kinetic species object.
    */
   virtual void eval( LevelData<FArrayBox>& result,
                     const KineticSpecies& kinetic_species ) const;
   
   /// Returns the kernel scale.
   /**
    * Returns the kernel scale, which is 1..
    *
    * @param[in] kinetic_species Kinetic species object.
    * @return real-valued scale to be applied in configuration space.
    */
   virtual Real scale( const KineticSpecies& kinetic_species ) const
   { return kinetic_species.charge(); }
   
   /// Returns the number of kernel components
   /**
    * Returns the number of kernel components.
    *
    * @param[in] kinetic_species Kinetic species object.
    */
   virtual int nComponents() const { return CFG_DIM; }
   
private:
   
   const LevelData<FluxBox>& m_field;
};


/// Guding Center Poloidal Momentum moment.
/**
 */
class GuidingCenterPoloidalMomKernel : public Kernel
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       GuidingCenterPoloidalMomKernel( const LevelData<FluxBox>& field );

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~GuidingCenterPoloidalMomKernel() {;}

      /// Computes the integrand for a PoloidalMom moment calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the scalar velocity squared times the radial component of the drift velocity.
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[out] result Cell-averaged integrand.
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }

    private:

        const LevelData<FluxBox>& m_field;
};

/// Magnetization kernel.
/**
 */
class MagnetizationKernel : public Kernel
{
public:
   
   /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
   virtual ~MagnetizationKernel() {;}
   
   /// Computes the integrand for a Magnetization moment calculation.
   /**
    * Computes the integrand by multiplying the input distribution function
    * by the magneitc moment.
    *
    * The number of kernel components is one.  If the distribution function
    * has more than one component, the first component of the result
    * will contain the moment of the first distribution function component,
    * the next component will contain the moment of the second
    * distribution function component, etc.  Upon entry, it is assumed
    * that result has already been initialized with the distribution function
    * correspondingly, ie., the first component contains the first
    * component of the distribution function, the next component
    * contain the second component of the distribution function, etc.
    *
    * @param[out] result Cell-averaged integrand.
    * @param[in] kinetic_species Kinetic species object.
    */
   virtual void eval( LevelData<FArrayBox>& result,
                      const KineticSpecies& kinetic_species ) const;
   
   /// Returns the kernel scale.
   /**
    * Returns the kernel scale, which is 1..
    *
    * @param[in] kinetic_species Kinetic species object.
    * @return real-valued scale to be applied in configuration space.
    */
   virtual Real scale( const KineticSpecies& kinetic_species ) const
   { return 1.; }
   
   /// Returns the number of kernel components
   /**
    * Returns the number of kernel components.
    *
    * @param[in] kinetic_species Kinetic species object.
    */
   virtual int nComponents() const { return 1; }
};

/// Delta F kernel.
/**
 */
class DeltaFKernel : public Kernel
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       DeltaFKernel(const CFG::LevelData<CFG::FArrayBox>& density,
                    const CFG::LevelData<CFG::FArrayBox>& temeprature,
                    const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift);


      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~DeltaFKernel() {;}

      /// Computes the cell-averaged difference the current dfn and Maxwellian function with the same first three moments
      /**
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[in] result Cell-averaged distribution function.
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }

    private:

        const CFG::LevelData<CFG::FArrayBox>& m_density;
        const CFG::LevelData<CFG::FArrayBox>& m_temperature;
        const CFG::LevelData<CFG::FArrayBox>& m_v_parallel_shift;

};

/// Maxwellian kernel.
/**
 */
class MaxwellianKernel : public Kernel
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       MaxwellianKernel(const CFG::LevelData<CFG::FArrayBox>& density,
                        const CFG::LevelData<CFG::FArrayBox>& temeprature,
                        const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift);


      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~MaxwellianKernel() {;}

      /// Computes the cell-averaged Maxwellian function with the same first three moments
      /**
       *
       * The number of kernel components is one.  If the distribution function
       * has more than one component, the first component of the result
       * will contain the moment of the first distribution function component,
       * the next component will contain the moment of the second
       * distribution function component, etc.  Upon entry, it is assumed
       * that result has already been initialized with the distribution function
       * correspondingly, ie., the first component contains the first
       * component of the distribution function, the next component
       * contain the second component of the distribution function, etc.
       *
       * @param[in] result Cell-averaged distribution function.
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<FArrayBox>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is 1..
       *
       * @param[in] kinetic_species Kinetic species object.
       * @return real-valued scale to be applied in configuration space.
       */
      virtual Real scale( const KineticSpecies& kinetic_species ) const
      { return 1.; }

      /// Returns the number of kernel components
      /**
       * Returns the number of kernel components.
       *
       * @param[in] kinetic_species Kinetic species object.
       */
       virtual int nComponents() const { return 1; }

    private:

        const CFG::LevelData<CFG::FArrayBox>& m_density;
        const CFG::LevelData<CFG::FArrayBox>& m_temperature;
        const CFG::LevelData<CFG::FArrayBox>& m_v_parallel_shift;

};



#include "NamespaceFooter.H"

#endif
