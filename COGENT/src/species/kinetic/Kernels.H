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

template <class T>
class Kernel
{
   public:
  
      /// virtual destructor to prevent compiler complaints (DFM 2/4/09)
      virtual ~Kernel<T>() {;}

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
      virtual void eval( LevelData<T>& result,
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
    
      ///Computes physical velocity components (e.g., V_R, V_Z) in configuration directions

      virtual void computeVelCfgComp(LevelData<T>&             vel_cfg_comp,
                                     const PhaseGeom&          phase_geom,
                                     const LevelData<FluxBox>& field,
                                     const int                 option) const;

      ///Computes projections of the physical velocity components onto normals in configuration directions

      virtual void computeVelCfgCompNormals(LevelData<T>&             vel_cfg_comp,
                                            const PhaseGeom&          phase_geom,
                                            const LevelData<FluxBox>& field,
                                            const int                 option) const;

};

/// Density kernel.
/**
 */
template <class T>
class DensityKernel : public Kernel<T>
{
   public:

      /// destructor to prevent compiler from complaining (DFM 2/4/09)
      virtual ~DensityKernel<T>() {;}

      /// Computes the kernel.
      /**
       * Computes the integrand, which is just one times the distribution
       * function.
       *
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
      */
      virtual void eval( LevelData<T>& result,
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
template <class T>
class ChargeDensityKernel : public Kernel<T>
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
      virtual void eval( LevelData<T>& result,
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
template <class T>
class MassDensityKernel : public Kernel<T>
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
      virtual void eval( LevelData<T>& result,
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
template <class T>
class MomentumDensityKernel : public Kernel<T>
{
   public:

      /// Constructor with initialization.
      /**
      *  Initialized with field, which is used to get phase-velocity
      */

      MomentumDensityKernel( const LevelData<FluxBox>& field )
         : m_field(field) {;}
      
      /// virtual destructor added to silence compiler complaints (MAD 1/5/18)
      virtual ~MomentumDensityKernel() {;}

      /// Computes the integrand for a guiding center momentum density moment (i.e., m*n*Vgc) calculation.
      /**
       * Computes the integrand by multiplying the input distribution function
       * by the physical guiding center velocity components and species mass
       * and returns the species mass as the scale.
       *
       * @param[in] result Cell-averaged distribution function.
       * @param[out] result Cell-averaged integrand. 
       * @param[in] kinetic_species Kinetic species object.
       */
      virtual void eval( LevelData<T>& result,
                         const KineticSpecies& kinetic_species ) const;

      /// Returns the kernel scale.
      /**
       * Returns the kernel scale, which is the species mass.
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
       */
      virtual int nComponents() const { return  CFG_DIM; }

private:
   
      const LevelData<FluxBox>& m_field;

};

/// 4th Moment kernel (used for transport models).
/**
 */
template <class T>
class FourthMomentKernel : public Kernel<T>
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
      virtual void eval( LevelData<T>& result,
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
template <class T>
class PressureKernel : public Kernel<T>
{
   public:

      PressureKernel(const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
         : m_v_parallel_shift(v_parallel_shift) {;}

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
      virtual void eval( LevelData<T>& result,
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

template <class T>
class ParallelPressureKernel : public Kernel<T>
{
   public:

      ParallelPressureKernel(const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
         : m_v_parallel_shift(v_parallel_shift) {;}

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~ParallelPressureKernel() {;}

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
      virtual void eval( LevelData<T>& result,
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

template <class T>
class PerpPressureKernel : public Kernel<T>
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~PerpPressureKernel() {;}

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
      virtual void eval( LevelData<T>& result,
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

/// Parallel heat flux kernel.
/**
 */
template <class T>
class ParallelHeatFluxKernel : public Kernel<T>
{
public:
   
   ParallelHeatFluxKernel(const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
      : m_v_parallel_shift(v_parallel_shift) {;}
   
   /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
   virtual ~ParallelHeatFluxKernel() {;}
   
   /// Computes the integrand for a ParallelHeatFlux moment calculation.
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
    * @param[in] result Cell-averaged distribution function.
    * @param[out] result Cell-averaged integrand.
    * @param[in] kinetic_species Kinetic species object.
    */
   virtual void eval( LevelData<T>& result,
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

/// Energy kernel.
/**
 */
template <class T>
class EnergyKernel : public Kernel<T>
{
   public:

      /// Constructor with initialization.
      /**
        *  Initialized with potentail, which is used to get phase-velocity
      */
      EnergyKernel( const CFG::LevelData<CFG::FArrayBox>& phi )
        : m_phi(phi) {;}
  
      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~EnergyKernel() {;}

      /// Computes the integrand for a energy moment calculation.
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
      virtual void eval( LevelData<T>& result,
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

      const CFG::LevelData<CFG::FArrayBox>& m_phi;

};

/// Kinetic energy kernel.
/**
 */
template <class T>
class KineticEnergyKernel : public Kernel<T>
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~KineticEnergyKernel() {;}

      /// Computes the integrand for a energy moment calculation.
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
      virtual void eval( LevelData<T>& result,
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

/// perp energy kernel.
/**
 */
template <class T>
class PerpEnergyKernel : public Kernel<T>
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
      virtual void eval( LevelData<T>& result,
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

/// Parallel Energy kernel.
/**
 */
template <class T>
class ParallelEnergyKernel : public Kernel<T>
{
   public:

      /// virtual destructor added to silence compiler complaints (DFM 2/4/09)
      virtual ~ParallelEnergyKernel() {;}

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
      virtual void eval( LevelData<T>& result,
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

/// Parallel momentum kernel.
/**
 */
template <class T>
class ParallelMomKernel : public Kernel<T>
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
      virtual void eval( LevelData<T>& result,
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
template <class T>
class ParticleFluxKernel : public Kernel<T>
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       ParticleFluxKernel(const LevelData<FluxBox>& field,
                          const int velocity_option)
          : m_field(field),
            m_velocity_option(velocity_option) {;}

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
      virtual void eval( LevelData<T>& result,
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
       virtual int nComponents() const { return 2; }

    private:

        const LevelData<FluxBox>& m_field;
        const int m_velocity_option;
};

/// Heat flux kernel.
/**
 */
template <class T>
class HeatFluxKernel : public Kernel<T>
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       HeatFluxKernel(const LevelData<FluxBox>& field,
                      const LevelData<FArrayBox>& phi,
                      const int velocity_option )
          : m_field(field),
            m_phi(phi),
            m_velocity_option(velocity_option) {;}

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
      virtual void eval( LevelData<T>& result,
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
       virtual int nComponents() const { return 2; }

    private:

        const LevelData<FluxBox>& m_field;
        const LevelData<FArrayBox>& m_phi;
        const int m_velocity_option;

};

/// Guding Center Poloidal Momentum moment.
/**
 */
template <class T>
class GuidingCenterPoloidalMomKernel : public Kernel<T>
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       GuidingCenterPoloidalMomKernel( const LevelData<FluxBox>& field )
          : m_field(field) {;}

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
      virtual void eval( LevelData<T>& result,
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
template <class T>
class MagnetizationKernel : public Kernel<T>
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
   virtual void eval( LevelData<T>& result,
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
template <class T>
class DeltaFKernel : public Kernel<T>
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       DeltaFKernel(const CFG::LevelData<CFG::FArrayBox>& density,
                    const CFG::LevelData<CFG::FArrayBox>& temperature,
                    const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
          : m_density(density),
            m_temperature(temperature),
            m_v_parallel_shift(v_parallel_shift) {;}

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
      virtual void eval( LevelData<T>& result,
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
template <class T>
class MaxwellianKernel : public Kernel<T>
{
   public:
       /// Constructor with initialization.
       /**
        *  Initialized with potentail, which is used to get phase-velocity
        */
       MaxwellianKernel(const CFG::LevelData<CFG::FArrayBox>& density,
                        const CFG::LevelData<CFG::FArrayBox>& temperature,
                        const CFG::LevelData<CFG::FArrayBox>& v_parallel_shift)
          : m_density(density),
            m_temperature(temperature),
            m_v_parallel_shift(v_parallel_shift) {;}

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
      virtual void eval( LevelData<T>& result,
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
