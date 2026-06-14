# GNC2: Self-Consistent Solution of Stellar Dynamics in Nuclear Star Clusters with a Central Massive Black hole

GNC2 is an updated version of the open-source **Monte Carlo code** [GNC](https://github.com/zhangfupeng-gzhu/GNC). It obtains self-consistent solutions for the stellar dynamics of nuclear star cluster (NSC) around a central Massive black hole (MBH). Developed by **Fupeng Zhang** and **Pau Amaro Seoane**. Compared to the previous version, GNC2 now includes the self-consistent stellar potential, adiabatic invariant theory, stellar evolution, and gravitational wave radiation acting on stellar orbits. 

**Note that currently it is in beta. The stable version will be released soon.**

-----

## Key Features

  * **Two-Dimensional Fokker-Planck Method:** Tracks the evolution of particles in energy and angular momentum space, which is highly efficient for simulating long-term relaxation.
  * **Multi-Mass Components:** Natively handles clusters containing different types of stars and compact objects (e.g., stellar-mass black holes, neutron stars).
  * **Core Physics Included:** Accurately models **two-body relaxation** and the effects of the **loss cone**, where objects are consumed by the central black hole.
  * **High Accuracy for Rare Particles:** Implements a weighting method to improve the statistical accuracy for rare but dynamically important objects.
  * **Flexible and Extensible:** The code is designed to be easily extended to include more complex physical processes, such as relativistic effects.
  * **Self-consistent Stellar Potential:** The stellar potential of stars is solved after each iteration of simulation according to Poisson Equation.
  * **Integrate Adiabatic Invariant Theory:** The energy of stellar objects is automatically adjusted according to adiabatic invariant theory, enabling self-consistent solutions of the stellar density profile with a slowly varying potential (including that of the central MBH). 
  * **Massive Black Hole Mass Growth:** Includes the mass accretion of MBH from loss cone accretion (e.g., tidal disruption of stars, direct swallowing of stellar objects), stellar evolution and from extreme mass ratio inspirals (EMRIs). 
  * **Gravitational Wave Radiation:** The energy and angular momentum of stellar objects decay due to gravitational wave radiation, for both bound and unbound orbits with respect to the central MBH.
  * **Integration of MOBSE Stellar Evolution Code:** The mass, stellar type, and physical radius of the stellar objects in the simulation are automatically updated according to MOBSE, an upgraded version of BSE code ( Hurley, J. R., Pols, O. R., & Tout, C. A. 2000, MNRAS, 315, 543; Hurley, J. R., Tout, C. A., & Pols, O. R. 2002, MNRAS, 329, 897 ): 

  	``MOBSE``: Giacobbo, N., Mapelli, M., & Spera, M. 2018, MNRAS, 474, 2959; Giacobbo, N., & Mapelli, M. 2018, MNRAS, 480, 2011.
  
-----

## Getting Started

### Prerequisites

The code is written in **Fortran 2003**. To compile and run it, you'll need: 

* `gfortran` (GNU fortran 13.3.0) with `Open MPI`

* Open MPI can installed via: apt install openmpi-bin libopenmpi-dev 

* HDF5 library is required: apt install libhdf5-dev 

### Installation and Usage
1. In your ``~/.bashrc``, added the following lines:
```shell
loadGNC() {
    local script_dir="$installed_gnc_source_path"   
    if [ -f "$script_dir/set_env.sh" ]; then
        pushd "$script_dir" > /dev/null
        source set_env.sh
		export HDF5LIBDIR=$installed_hdf5lib_path
		export HDF5INCDIR=$installed_hdf5inc_path
        popd > /dev/null
    else
        echo "error：can not find $script_dir/set_env.sh"
    fi
}


```
Replace:
* ``$installed_gnc_source_path`` with the location of the source files you downloaded.
* ``$installed_hdf5lib_path`` and ``$installed_hdf5inc_path`` with the library and include paths of HDF5, which can be obtained by running ``h5fc -show''.

2. Change to the source directory and run ``make". Wait until compilation is complete. In Linux systems where stack memory execution is prohibited, you need to set ``usestack=1`` in ``main/makefile`` and install ``execstack`` to mark the program's stack as executable.

3. *Start a new shell*, run the following commands in shell :
```shell
> loadGNC
```

4. You are now ready to run examples to verify your installation. Change to any subfolder of the ``examples`` folder, then run:
```shell
> mpirun -np number_of_threads ini
```
for model initialization, followed by
```shell
> mpirun -np number_of_threads main
```
for the main simulation run. Here ``number_of_threads`` can be 4, 8, 12, 24, or 48, depending on your machine.

For models that include gravitational wave radiation, you can run 
```shell
>event_analysis
```
after the simulation finishes to obtain EMRI data. The output is located in the folder : output/ecev/dms/event_analysis/

-----

## Examples

* ``M1`` is a cluster consists of $1$ solar mass stars. The density follow Dehnen profile, with a total mass of $4\times10^7$ solar mass, $r_a=2.17$ pc and $\gamma=1$. The spin effect of MBH is ignored. The mass of MBH $=10^7$ solar mass and always remain fixed. In this example, there is **No** stellar evolution and mass accretion.
The simulation usually finishes within 1 hours when running with 4 threads. After that, you can plot the evolution of the density profile, tidal disruption event rates, effective radius of the cluster, and the influence radius of the MBH in ``plot'' by running:
```
> loadGNC
> python3 plot.py
```

* ``M2_sp`` is a cluster consists of two mass components (stars + a number fraction of 0.001 stellar mass black holes). The density follow Dehnen profile, with a total mass of $4\times10^7$ solar mass, $r_a=2.17$ pc and $\gamma=1$. The spin of MBH =1, the mass of MBH $=4\times10^6$ solar mass and always remain fixed. In this example, there is **No** stellar evolution and mass accretion.
The simulation usually takes about 3 hours when running with 24 threads. After that, you can plot the evolution of the tidal disruption and EMRI event rates in ``plot'' by running:
```
> loadGNC
> python3 rates.py
```

* ``MKG82_sp_stacc_0.06`` is a cluster consists of 12 mass components ranges from 0.01 to 150 solar masses. Initially all particles are zero-age main-sequence stars, follow the Kroupa initial mass function. The density follow Dehnen profile, with a total mass of $7\times10^7$ solar mass, $r_a=1.5$ pc and $\gamma=1$. The spin of MBH =1, the initial mass of MBH is $=10^4$ solar mass. Stellar evolution and mass accretion of MBH are **included**. The simulation usually takes about 18 hours when running with 24 threads. After that, you can plot the evolution of the MBH mass growth, tidal disruption and EMRI event rates in ``plot'' by running:
```
> loadGNC
> python3 evl.py
```

-----

## Setting Parameters of the Simulation

When run ``ini`` or ``main``, you need to have a parameter file ``model.in`` in the folder, which is the custom parameters of the simulation; There is a list of all default parameters in ``source/main/model_default.in``. Parameters set in ``model.in`` will overwrite the default parameters. A simple example of ``model.in`` is shown in ``examples/M2_sp``, as following: 
```
##================================================================
#Specify the path of model_default.in respect to the root of GNC folder
model_default_path= /source/main/model_default.in
#lines starts with ``#`` is considered as comments
#Set the unit of mass to 4*10^6msun
mass unit =4e6
#Consider a maximally spinning MBH, so that the ISO orbit will be modified according to the orbital inclination
mbh spin = 1
#Includes the gravitational wave orbital radiation.
GW outer orbit  =1
--EMRI criteria  =0.001
--output EMRI data bin =1
#Set the maximum value of energy. Usually, 10^5 or 10^6 is enough.
set emax=1d6
#Set the inner and outer boundary of distance to the MBH. Usually, it can be ``1d-5 1000`` or ``1d-6  1000``
rmin rmax = 1d-6 1000
#Set the chattery of simulation output. ``0`` has the minimum output, ``1`` has more details, values larger than ``2`` are for debug use.
chattery=0
```
For the complete list of all available parameters, see ``source/main/model_default.in``.

## Setting Initial Parameters of the Cluster

In the same folder, ``mfrac.in`` sets the initial value of parameters of cluster's mass components and density profiles. An example of a ``Denhen`` cluster with two discrete mass component (stars+ stellar mass black holes) are :
```
#number of mass bin 
2  GIVEN  
###########
sg data mode = given
model mode = SAME
#first line: m1 mc m2 parameter_of_weighting_of_particles clone_factor
#second line: star_fraction  sbh_fraction  ns_fraction wd_fraction bd_fraction, which is the fraction of components within each mass bin
#``star'': MS stars
#``sbh'': stellar black holes
#``ns``: neutron stars
#``wd``: white dwarfs
#``bd``: brown dwarfs
----------------------------------------------
   1.00   1.00    1.00    100000 30
   1.00   0.00    0.00    0.0  0.0  0.0   
#----------------------------------------------
   10.00   10.00    10.00   100 8
   0.00   1.00    0.00    0.0  0.0  0.0   
#----------------------------------------------
# the following line sets a 10^-3 fraction of stellar black holes respective to stars
1    0.001   
##########################################
density model = Dehnen
unit= MsunPc
#mtot, in unit of msun
mtot=4e7
#ra, in unit of pc
ra=2.17
gamma=1.0
----------------------------------------------
```

**Setting `particle_weighting_factor`**
When particle cloning is ignored, ``particle_weighting_factor`` is the real number of stellar objects represented by each simulation particle. In ``GNC``, simulations should typically use the cloning scheme by setting ``clone_scheme=1`` in ``model.in``. In such cases, the weight of a particle can be further split to a number of ``clone_factor`` clone particles, depending on its energy. Thus, `particle_weighting_factor` and `clone_factor` together determine the total number of Monte-Carlo particles initially generated for a given mass bin. The larger `particle_weighting_factor`, the fewer particles in simulation, because each particle now represents more real objects.

Note that ``particle_weighting_factor`` values across different mass bins must all be integer multiples of the smallest ``particle_weighting_factor`` among them.

**Setting `clone_factor`** clone factor sets the amplification of clone particles at higher energies. Due to the mass segregation, the following values are recommended:
for stars or smaller mass particles, `clone_factor=20~40`;
for stellar mass black holes or more massive particles, `clone_factor = 4 ~ 10`;

If you find that your simulation initially generates too many particles in a mass bin, make ``clone_factor`` smaller or ``particle_weighting_factor`` larger, then try again.

If the density profile initially generated does not resolve the inner regions of the cluster well, make ``clone_factor`` larger.

**``particle_weighting_factor`` and ``clone_factor`` should be set, so that in each mass bin, the total number of particles per thread is around 1000~100000**. Depending on your machine, the total number of particles across all mass bins should be smaller than $2\times10^6$, unless you have a more than 32 GB of memory. For example, when running M2_sp, you will get lines of output like the following:
```
              NStar           NSbh            NNs            NWd            NBd            
TOT=           9847           2331              0              0              0             
   1           9847              0              0              0              0              
   2              0           2331              0              0              0             
```
It means that in the first mass bin, the number of stars is ``9847`` and in the second mass bin the number of SBHs is ``2331``.

The second example of cluster is a cluster consists of 12 mass bins, following Kroupa Initial mass function:
```
#number of mass bin   
12  KROUPA
#Use MOBSE stellar evolution code
sg data mode  = mobse
model mode    =SAME
----------------------------------------------
metalicity = 0.02
#line for each mass bin: 
# bin_number minimum_mass maximum_mass parameters_of_particle_weighting clone_factor
1   0.01   0.1  50000  20
----------------------------------------------
2   0.1    0.5  50000  20
----------------------------------------------
3   0.5    1.0  30000  16
----------------------------------------------
4   1      2   5000  16
----------------------------------------------
5   2      4   5000  16
----------------------------------------------
6   4      8   2000   8
----------------------------------------------
7   8      12  1000   8
----------------------------------------------
8   12     24  400   6
----------------------------------------------
9   24     32  200   6
----------------------------------------------
10   32    64  200   6
----------------------------------------------
11   64    128  100   6
----------------------------------------------
12   128   150  10   4
----------------------------------------------
###########
density model = Dehnen
unit =MsunPc
#mtot, in unit of msun
mtot=7d7
#ra, in unit of pc
ra=1.5
gamma=1.0
----------------------------------------------
```
Particles in 12 mass bins all follow the Dehnen density profile (with different normalization).

Similarly, for different Dehnen profiles, ``parameters_of_particle_weighting`` and ``clone_factor`` should be adjusted so that the number of particles in each mass bin does not have too many or too few samples.

-----

## Citing GNC2

If you have used GNC2 in your work, we kindly ask you to please cite the original reference paper:

* Zhang, F., & Amaro-Seoane, P. 2026, The Astrophysical Journal, 999, 224
* Zhang, F., & Amaro-Seoane, P. 2025, The Astrophysical Journal, 980, 210
* Zhang, F., & Amaro-Seoane, P. 2024, The Astrophysical Journal, 961, 232

-----

## The other source codes used in GNC2

We have used the following codes from other people, 

* DOPRI5, Dormand, J. R., & Prince, P. J. 1980, JCoAM, 6, 19; 
* The Gaussian integrator has used the ``SLATEC Common Mathematical Library`` by Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro and Lee Walton;
* MOBSE: Giacobbo, N., Mapelli, M., & Spera, M. 2018, MNRAS, 474, 2959; Giacobbo, N., & Mapelli, M. 2018, MNRAS, 480, 2011.

-----

## Questions and Support

If you have any questions, find a bug, or have suggestions for improvements, please feel free to contact the authors:

  * **Fupeng Zhang:** `zhangfupeng@gzhu.edu.cn`
  * **Pau Amaro Seoane:** `amaro@upv.es`

  


