#ifndef PARTICLES_H
#define PARTICLES_H

#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>

#include "Tools.h"
#include "TimeSelection.h"

class Particle;

class Params;
class Patch;



//----------------------------------------------------------------------------------------------------------------------
//! Particle class: holds the basic properties of a particle
//----------------------------------------------------------------------------------------------------------------------
class Particles
{
public:

    //! Constructor for Particle
    Particles();

    //! Destructor for Particle
    virtual ~Particles();

    //! Create nParticles null particles of nDim size
    void initialize( unsigned int nParticles, unsigned int nDim, bool keep_position_old );

    //! Create nParticles null particles of nDim size
    void initialize( unsigned int nParticles, Particles &part );

    //! Set capacity of Particles vectors
    void reserve( unsigned int n_part_max, unsigned int nDim );

    //! Initialize like another particle, but only reserve space
    void initializeReserve( unsigned int n_part_max, Particles &part );

    //! Resize Particles vectors
    void resize( unsigned int nParticles, unsigned int nDim, bool keep_position_old );

    //! Resize Particles vectors
    void resize( unsigned int nParticles);

    //! Resize the cell_keys vector
    void resizeCellKeys(unsigned int nParticles);

    //! Remove extra capacity of Particles vectors
    void shrinkToFit();

    //! Reset Particles vectors
    void clear();

    //! Get number of particules
    inline unsigned int size() const
    {
        return Weight.size();
    }

    //! Get number of particules
    inline unsigned int capacity() const
    {
        return Weight.capacity();
    }

    //! Get dimension of particules
    inline unsigned int dimension() const
    {
        return Position.size();
    }

    //! Copy particle iPart at the end of dest_parts
    void copyParticle( unsigned int iPart, Particles &dest_parts );
    //! Copy particle iPart at the end of the current array
    void copyParticle( unsigned int iPart );
    //! Insert particle iPart at dest_id in dest_parts
    void copyParticle( unsigned int ipart, Particles &dest_parts, int dest_id );

    //! Insert nPart particles starting at ipart to dest_id in dest_parts
    void copyParticles( unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id );

    //! Make a new particle at the position of another
    void makeParticleAt( Particles &source_particles, unsigned int ipart, double w, short q=0., double px=0., double py=0., double pz=0. );

    //! Suppress particle iPart
    void eraseParticle( unsigned int iPart );
    //! Suppress nPart particles from iPart
    void eraseParticle( unsigned int iPart, unsigned int nPart );

    //! Suppress all particles from iPart to the end of particle array
    void eraseParticleTrail( unsigned int iPart );

    //! Print parameters of particle iPart
    void print( unsigned int iPart );

    friend std::ostream &operator << ( std::ostream &, const Particles &particle );

    //! Exchange particles part1 & part2 memory location
    void swapParticle( unsigned int part1, unsigned int part2 );
    void swapParticles( std::vector<unsigned int> parts );
    void translateParticles( std::vector<unsigned int> parts );
    void swapParticle3( unsigned int part1, unsigned int part2, unsigned int part3 );
    void swapParticle4( unsigned int part1, unsigned int part2, unsigned int part3, unsigned int part4 );

    //! Exchange particles part1 & part2 memory location
    void swapParticle( unsigned int part1, unsigned int part2, unsigned int N );

    //! Overwrite particle part1 into part2 memory location. Erasing part2
    void overwriteParticle( unsigned int part1, unsigned int part2 );

    //! Overwrite particle part1->part1+N into part2->part2+N memory location. Erasing part2->part2+N
    void overwriteParticle( unsigned int part1, unsigned int part2, unsigned int N );

    //! Overwrite particle part1->part1+N into part2->part2+N of dest_parts memory location. Erasing part2->part2+N
    void overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2, unsigned int N );

    //! Overwrite particle part1 into part2 of dest_parts memory location. Erasing part2
    void overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2 );

    //! Move iPart at the end of vectors
    void pushToEnd( unsigned int iPart );

    //! Create new particle
    void createParticle();

    //! Create n_additional_particles new particles
    void createParticles( int n_additional_particles );

    //! Create n_additional_particles new particles at position pstart in the particles data structure
    void createParticles( int n_additional_particles, int pstart );

    //! Move ipart at new_pos in the particles data structure
    void moveParticles( int iPart, int new_pos );

    //! Compress the particles vectors according to the provided mask
    //! between istart and iend
    void eraseParticlesWithMask( int istart, int iend, std::vector <int> & mask );

    //! Compress the particles vectors using cell_keys as a mask
    //! between istart and iend
    void eraseParticlesWithMask( int istart, int iend);

    //! This method erases particles according to the provided mask
    //! between istart and iend
    // void eraseParticlesWithMask( int istart, int iend, vector <bool> & to_be_erased);

    //! Test if ipart is in the local patch
    bool isParticleInDomain( unsigned int ipart, Patch *patch );

    //! Method used to get the Particle position
    inline double  position( unsigned int idim, unsigned int ipart ) const
    {
        return Position[idim][ipart];
    }
    //! Method used to set a new value to the Particle former position
    inline double &position( unsigned int idim, unsigned int ipart )
    {
        return Position[idim][ipart];
    }

    //! Method used to get the Particle position
    inline double distance2ToAxis( unsigned int ipart ) const
    {
        return Position[1][ipart] * Position[1][ipart] + Position[2][ipart] * Position[2][ipart];
    }

    //! Method used to get the Particle position
    inline double  position_old( unsigned int idim, unsigned int ipart ) const
    {
        return Position_old[idim][ipart];
    }
    //! Method used to set a new value to the Particle former position
    inline double &position_old( unsigned int idim, unsigned int ipart )
    {
        return Position_old[idim][ipart];
    }

    //! Method used to get the list of Particle position
    inline std::vector<double>  position( unsigned int idim ) const
    {
        return Position[idim];
    }

    //! Method used to get the Particle momentum
    inline double  momentum( unsigned int idim, unsigned int ipart ) const
    {
        return Momentum[idim][ipart];
    }
    //! Method used to set a new value to the Particle momentum
    inline double &momentum( unsigned int idim, unsigned int ipart )
    {
        return Momentum[idim][ipart];
    }
    //! Method used to get the Particle momentum
    inline std::vector<double>  momentum( unsigned int idim ) const
    {
        return Momentum[idim];
    }

    //! Method used to get the Particle weight
    inline double  weight( unsigned int ipart ) const
    {
        return Weight[ipart];
    }
    //! Method used to set a new value to the Particle weight
    inline double &weight( unsigned int ipart )
    {
        return Weight[ipart];
    }
    //! Method used to get the Particle weight
    inline std::vector<double>  weight() const
    {
        return Weight;
    }
    //! Method used to get the Acc_path
    inline double  acc_path( unsigned int ipart ) const
    {
        return Acc_path[ipart];
    }
    //! Method used to set a new value to the Acc_path
    inline double &acc_path( unsigned int ipart )
    {
        return Acc_path[ipart];
    }
    //! Method used to get the Acc_path
    inline std::vector<double>  acc_path() const
    {
        return Acc_path;
    }

    //! Method used to get the Particle charge
    inline short  charge( unsigned int ipart ) const
    {
        return Charge[ipart];
    }
    //! Method used to set a new value to the Particle charge
    inline short &charge( unsigned int ipart )
    {
        return Charge[ipart];
    }
    //! Method used to get the list of Particle charges
    inline std::vector<short>  charge() const
    {
        return Charge;
    }


    //! Method used to get the Particle Lorentz factor
    inline  double LorentzFactor( unsigned int ipart )
    {
        return sqrt( 1.+pow( momentum( 0, ipart ), 2 )+pow( momentum( 1, ipart ), 2 )+pow( momentum( 2, ipart ), 2 ) );
    }

    //! Method used to get the inverse Particle Lorentz factor
    inline double inverseLorentzFactor( unsigned int ipart )
    {
        return 1./sqrt( 1.+ momentum( 0, ipart )*momentum( 0, ipart ) + momentum( 1, ipart )*momentum( 1, ipart ) + momentum( 2, ipart )*momentum( 2, ipart ) );
    }

    //! Method used to get the momentum norm which is also the normalized photon energy
    inline double momentumNorm( unsigned int ipart )
    {
        return sqrt( pow( momentum( 0, ipart ), 2 )+pow( momentum( 1, ipart ), 2 )+pow( momentum( 2, ipart ), 2 ) );
    }

    //! Partiles properties, respect type order : all double, all short, all unsigned int

    //! array containing the particle position
    std::vector< std::vector<double> > Position;

    //! array containing the particle former (old) positions
    std::vector< std::vector<double> >Position_old;

    //! array containing the particle moments
    std::vector< std::vector<double> >  Momentum;
    //! array containing the particle accumulated path
    std::vector<double> Acc_path;

    //! containing the particle weight: equivalent to a charge density
    std::vector<double> Weight;

    //! containing the particle quantum parameter
    std::vector<double> Chi;

    //! Incremental optical depth for the Monte-Carlo process
    std::vector<double> Tau;

    //! charge state of the particle (multiples of e>0)
    std::vector<short> Charge;

    //! Id of the particle
    std::vector<uint64_t> Id;

    //! cell_keys of the particle
    std::vector<int> cell_keys;

    // TEST PARTICLE PARAMETERS
    bool is_test;

    //! True if tracking the particles
    bool tracked;

    void resetIds()
    {
        unsigned int s = Id.size();
        for( unsigned int iPart=0; iPart<s; iPart++ ) {
            Id[iPart] = 0;
        }
    }

    //! Method used to get the Particle Id
    inline uint64_t id( unsigned int ipart ) const
    {
        DEBUG( ipart << " of " << Id.size() );
        return Id[ipart];
    }
    //! Method used to set the Particle Id
    inline uint64_t &id( unsigned int ipart )
    {
        return Id[ipart];
    }
    //! Method used to get the Particle Ids
    inline std::vector<uint64_t> id() const
    {
        return Id;
    }
    void sortById();

    //! Quantum parameter for particles that are submitted
    //! to a radiation reaction force (CED or QED)
    bool isQuantumParameter;

    //! Parameters for particles that are submitted to a
    //! Monte-Carlo process such as:
    //! - discontinuous radiation reaction force
    bool isMonteCarlo;

    //! Method used to get the Particle chi factor
    inline double  chi( unsigned int ipart ) const
    {
        return Chi[ipart];
    }
    //! Method used to set a new value to the Particle chi factor
    inline double &chi( unsigned int ipart )
    {
        return Chi[ipart];
    }
    //! Method used to get the Particle chi factor
    inline std::vector<double>  chi() const
    {
        return Chi;
    }

    //! Method used to get the Particle optical depth
    inline double  tau( unsigned int ipart ) const
    {
        return Tau[ipart];
    }
    //! Method used to set a new value to
    //! the Particle optical depth
    inline double &tau( unsigned int ipart )
    {
        return Tau[ipart];
    }
    //! Method used to get the Particle optical depth
    inline std::vector<double>  tau() const
    {
        return Tau;
    }

    //! Method to keep the positions for the next timesteps
    void savePositions();

    std::vector< std::vector<double  >*> double_prop_;
    std::vector< std::vector<short   >*> short_prop_;
    std::vector< std::vector<uint64_t>*> uint64_prop_;

    //! Specific pointers
    double * __restrict__ position_x;
    double * __restrict__ position_y;
    double * __restrict__ position_z;

    double * __restrict__ momentum_x;
    double * __restrict__ momentum_y;
    double * __restrict__ momentum_z;

#ifdef __DEBUG
    bool testMove( int iPartStart, int iPartEnd, Params &params );

    inline double dist2( unsigned int iPart )
    {
        double dist( 0. );
        for( unsigned int iDim = 0 ; iDim < Position.size() ; iDim++ ) {
            double delta = position( iDim, iPart )-position_old( iDim, iPart );
            dist += delta*delta;
        }
        return dist;
    }
    inline double dist( unsigned int iPart, unsigned int iDim )
    {
        double delta = std::abs( position( iDim, iPart )-position_old( iDim, iPart ) );
        return delta;
    }
#endif

    Particle operator()( unsigned int iPart );

    //! Methods to obtain any property, given its index in the arrays double_prop_, uint64_prop_, or short_prop_
    void getProperty( unsigned int iprop, std::vector<uint64_t> *&prop )
    {
        prop = uint64_prop_[iprop];
    }
    void getProperty( unsigned int iprop, std::vector<short> *&prop )
    {
        prop = short_prop_[iprop];
    }
    void getProperty( unsigned int iprop, std::vector<double> *&prop )
    {
        prop = double_prop_[iprop];
    }

    //! Indices of first and last particles in each bin/cell
    std::vector<int> first_index, last_index;

    virtual void initGPU() { std::cout << "Should not came here" << std::endl; };
    virtual void syncGPU() { std::cout << "Should not came here" << std::endl; };
    virtual void syncCPU() { std::cout << "Should not came here" << std::endl; };
    virtual double* getPtrPosition( int idim ) {
        return &(Position[idim][0]);
    };
    virtual double* getPtrPosition_old( int idim ) {
        return &(Position_old[idim][0]);
    };
    virtual double* getPtrMomentum( int idim ) {
        return &(Momentum[idim][0]);
    };
    virtual double* getPtrWeight() {
        return &(Weight[0]);
    };
    virtual double* getPtrAccPath() {
        return &(Acc_path[0]);
    };
    virtual double* getPtrChi() {
        return &(Chi[0]);
    };
    virtual short* getPtrCharge() {
        return &(Charge[0]);
    };
    virtual uint64_t* getPtrId() {
        return &(Id[0]);
    };
    virtual double* getPtrTau() {
        return &(Tau[0]);
    };
    virtual int* getPtrCellKeys() {
        return &(cell_keys[0]);
    };


private:

};

#endif
