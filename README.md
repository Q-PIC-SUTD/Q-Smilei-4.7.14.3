<img src="doc/Sphinx/_static/smileiLogo.svg" width=300 />

# Particle-In-Cell code for plasma simulation

Open-source, user-friendly and designed for high performances on super-computers, it is applied to a wide range of physics studies: from relativistic laser-plasma interaction to astrophysics.

## Documentation [![Build Status](https://travis-ci.com/SmileiPIC/Smilei.svg?branch=master)](https://travis-ci.com/SmileiPIC/Smilei)

* [General documentation](https://smileipic.github.io/Smilei) 
* [Tutorials](https://smileipic.github.io/tutorials)

## Feedback

* [Issues](https://github.com/SmileiPIC/Smilei/issues): bug reports, feature requests, reporting on bad doc or unexpected behaviour.
* [Chat room](https://app.element.io/#/room/!LQrdVpOJEohPSWMlmf:matrix.org): general discussions, suggestions, remarks, sharing results & papers or just to say hi!
* [Contribute](https://smileipic.github.io/Smilei/contribute.html).

Copied from Smilei-4.7.3.7
Customized:
- epsilon
- mfp/mean free time
- In PusherBoris_mfp and PusherBoris_mfp_2D, particles are deflected before pushed by fields
- Const external field not involved in Maxwell's solver
- In external field, do not use it to apply B fields
- particle pusher from prescribed oscillation amplitude and frequency. So far only oscillation along x at 0 phase (cos(wt)) is available

- particle pusher for massless electrons ("massless", must set mass>0 still, else it will be assumed to be photon with no charge. Calculations involving mass are dummy (*m then /m), provide a constant velocity through V_const). So far only available for quasistatic case (no rotation in magnetic field)
- Quasistatic solver (without Faraday equation) is available with keyword "quasistatic" in Main 
- nSurf field (defined similar as Rho) added
- Interpolator for nSurf (0 order) available only for 2D simulation without vectorization (momentum conserving)
- Boris Pusher with nSurf available
- particle pusher from prescribed oscillation amplitude and frequency, with time-dependent amplitude (prescribed_profile). So far only built-in time profile is available, user-define python function result in segmentation fault.
- Add nsurf Interpolator for 3D (no vectorization)
- surface reflection is included in nSurf interpolator for both 2D and 3D. Reflect on both momentum and position. The position of surface is assumed to be the mid-point of position and position_old. Field interpolator is done for the newly-reflected position.

- quasistatic solver, laser applied as prescribed field. In PrescribedField, use keyword 'is_it_uniform' to save time applying prescribed field as this is done on 1 processor
- Dirac electron in graphene pusher needs to be re-configured with a separate momentum and velocity. In original code, momentum is treated as velocity to deposit in the current. Need to add another momentum parameter (momentum_raw) that participates in the particle pusher. Momentum = V_const*momentum_raw/(abs(momentum_raw)). mass cannot set =0 as the code treat mass=0 for photons. 
- Add Massless pusher with mfp

- in nSurf interpolator add diffuse scattering: Bechmarked for Interpolator2D0OrderWnSurf_diffuse.cpp, Pending benchmarking for Interpolator3D0OrderWnSurf_diffuse.cpp
Changes in this version:
- Diffuse scattering has a cosine distribution around nsurf vector
- Diffuse scattering with specularity = probability of specular scattering

