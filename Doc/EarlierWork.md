## Early models:
Those models focused on a particular phenomenon and animated the smokes density directly without modelling its velocity.
The implementation involves density modelling, animating density, and adding detail.

Density modelling is used to determine the distribution of smoke in space, which could be achieved by
defining a three-dimensional grid and storing density values at each grid point. 
Initially, some regions might have higher density values to simulate the initial shape of the smoke. 
After modelling the density, the models then animate the density to simulate its movement, which is implemented 
by changing the density values at each grid point. For example, density values could be increased
in the direction of smoke movement to simulate the effect of smoke moving in a certain direction. 
To increase the realism of the simulation, additional details such as textures could be added. 
These details could include textures, changes in colour, or other forms of noise to make the smoke appear more realistic. 

Those models are relatively simple and have their limitations. The biggest problem is that those models
only focused on changes in smoke density without considering its velocity. This leads to simulation results
lacking dynamicity and failing to accurately capture the movement characteristics of smoke. 
They cant accurately simulate the curling and swirling motions typical of smoke since only density was considered.


## Models with Velocity Fields:
Based on early models, the model with velocity fields marks a significant advancement in smoke modelling. 
It is implemented by adding three-dimensional arrays where each cell contains a velocity vector representing
the speed and direction of motion at that point in space to represent the motion of smoke. 

Initially, these velocity fields are populated with random values to create a turbulent and chaotic initial motion,
resembling the behaviour of real-world smoke. They are then used to transport the smoke density. 
This means that the smoke density is moved along the velocity field, causing it to spread and deform over time. 

However, even this model has its limitations; the model lacks dynamical feedback even though it captures
the general motion characteristics of smoke. Realistic smoke behaviour often involves interactions between density and velocity,
such as buoyancy effects and turbulence generation, which are not fully captured by this approach. 
Additionally, implementing and simulating random velocity fields can be computationally expensive, 
especially when dealing with high-resolution grids and complex scenes. Generating and advecting these fields requires significant 
computational resources, limiting the scale and detail of simulations that can be achieved. 
Furthermore, since the initial velocity fields are random, fine-tuning and controlling the behaviour of the simulated smoke
can be challenging. Achieving desired effects may require manual adjustment of parameters, 
making it less intuitive compared to physically based simulation methods.

## Direct Simulation of Fluid Dynamics (CFD):
CFD is a more natural approach which solves the Navier-Stokes equations directly,
which govern the behaviour of fluids to simulate the motion of smoke. The Navier-Stokes equations describe 
the motion of fluid substances. These equations, coupled with appropriate boundary conditions, 
govern the evolution of fluid velocity and pressure fields over time. For smoke simulation, 
these equations are solved numerically using computational methods. 

The simulation domain is discretized into a three-dimensional grid or mesh. The Navier-Stokes equations are then
solved at discrete points or cells within this grid, typically using finite difference, finite volume, 
or finite element methods. Additionally, time-stepping methods, such as the Euler method or implicit methods like 
the Runge-Kutta method, are used to integrate the equations forward in time. And the time step size is chosen based on 
stability and accuracy considerations. Even some boundary conditions are imposed to model interactions between the fluid 
and its surroundings. These conditions specify the behaviour of the fluid at the boundaries of the simulation domain, 
such as inflow, outflow, or solid walls. It also employed turbulence models to simulate the effects of small-scale turbulent
motion that cannot be resolved directly by the computational grid. These models add additional terms to the Navier-Stokes 
equations to capture turbulent effects. Finally, efficient numerical solvers are used to solve the discretized equations, 
considering the complex interplay between pressure, velocity, and viscosity terms. 

Comparing with earlier described models, despite this gives a better simulation result, it is expensive to do high-resolution
simulations or complex geometries. The need to solve the Navier-Stokes equations at each grid point over time requires
significant computational resources. In addition, generating and managing the computational mesh can be challenging, 
especially for irregular or complex geometries. Mesh quality and resolution directly impact the accuracy and stability 
of the simulation. Also, the turbulence modelling introduces additional complexity and uncertainty.
One more factor that may lead to inaccuracies in the simulation results is that there are several assumptions about 
fluid behaviour with Navier-Stokes equations such as incompressibility and Newtonian viscosity. 
Which may not always hold true for real-world fluids like smoke. 

## Further development by Foster and Metaxas:
Later, Foster and Metaxas made significant improvement in smoke simulation by simulating smoke motions in three dimensions
using relatively coarse grids. Instead of needing an extremely fine grid to capture all the details, 
they could use a more manageable grid resolution while still producing visually convincing results. 
Their simulations captured the characteristic swirling motions of smoke, adding realism to the simulations. 
Which was achieved through careful modelling and computation of the fluid dynamics equations.

The limitations of this approach lies partly on limited detail due to the fact that coarse grids inherently limit the
level of detail that can be captured in simulations. While they were able to produce convincing swirling motions, 
finer details, such as small-scale vortices typical of smoke, may have been lost due to the coarse resolution.

Additionally, their explicit integration scheme required small time steps to maintain stability, particularly when
fluid velocities were high. This led to slower simulations, as larger time steps couldn't be used without risking 
instability even though their simulations produced visually appealing results. Furthermore, like many explicit 
integration schemes, their simulations suffered from numerical dissipation. This means that small-scale details in
the flow, such as fine vortices, tended to dissipate too quickly, resulting in a loss of realism in the simulation.

## Stam's Model:
So far, all the models have stability problems due to different factors. Stam's model made achievement on unconditional stability.
Unlike previous approaches that required small time steps for stability, Stam's model could maintain stability even with
larger time steps. This greatly improved the efficiency of simulations, allowing them to run faster without sacrificing
stability. Stam employed a semi-Lagrangian advection scheme, which involves tracing back the flow of the fluid to advect
quantities like density and velocity. This method is computationally efficient and helps reduce numerical dissipation 
compared to some other advection schemes. He even used implicit solvers to solve the equations governing fluid motion. 
Which are more stable and can handle larger time steps compared to explicit solvers, contributing to the model's stability.

However, it still had limitations, despite its stability, Stam's model suffered from numerical dissipation. 
This means that small-scale features in the flow, such as fine vortices and details in the smoke, could be overly 
smoothed out, leading to a loss of realism in the simulation. Stam's model used a first-order integration scheme, 
which limited its ability to accurately capture the finer details of smoke motion. Higher-order integration schemes 
could provide better accuracy but might introduce additional complexity and computational cost. Like other named 
simulation models,  Stam's model required careful tuning of parameters to achieve the desired results. 
Finding the right balance between stability, accuracy, and computational efficiency could be challenging 
and time-consuming. For implementing Stam's model, particularly the combination of semi-Lagrangian advection 
and implicit solvers, could be complex and require expertise in numerical methods and computational fluid dynamics.


## Compressible flow equations VS Incompressible flow equations
compressible flow equations and incompressible flow equations are two fundamental frameworks for modeling fluid dynamics,
each suitable for different scenarios depending on the magnitude of density variations within the fluid.
Yngve et al. proposed using compressible equations to model explosions, but these introduced strict time step restrictions.
Most practitioners prefer using incompressible equations due to their flexibility and avoidance of strict time step 
restrictions associated with compressible equations. The key distinction between compressible and incompressible flow is
the treatment of density variation. Compressible flow equations account for density changes, while incompressible flow 
equations assume constant density. Compressible flow equations require consideration of the speed of sound, 
whereas incompressible flow equations typically ignore its effects. Compressible flow equations are used in scenarios
where density changes are significant, such as high-speed aerodynamics and combustion. 
Incompressible flow equations are used in scenarios where density changes are negligible, 
such as most everyday fluid flow applications.

## Lattice Gas Solvers (LGS)
Another alternative Lattice Gas Solvers(LGS), LGS are based on the concept of cellular automata, where the fluid is
represented as discrete particles moving on a lattice grid. Each cell in the grid represents a small volume of fluid,
and particles move between neighboring cells according to predefined rules. Unlike continuum-based methods like finite
difference or finite volume, LGS simulate fluid dynamics at a microscopic level. They track the motion and interactions
of individual particles to simulate the behavior of the entire fluid.

The simulation proceeds in two main steps: collision and streaming. During collision, particles interact within
each cell according to collision rules. Then, during streaming, particles move to neighboring cells based on streaming
rules. Despite the simple rules governing particle motion, complex fluid behavior emerges from the collective
interactions of the particles. This allows LGS to capture phenomena like turbulence, shock waves, and vortices.
LGS can be highly parallelized, making them well-suited for implementation on parallel computing architectures like
GPUs. This enables efficient simulation of large-scale fluid flows.

while Lattice Gas Solvers offer a computationally efficient and conceptually simple approach to simulating fluid dynamics,
they also have limitations. LGS are inherently less accurate compared to continuum-based methods for fluid dynamics, particularly in modeling complex fluid behaviors or capturing fine-scale details. This is because they rely on simplified rules and discrete particle motion rather than solving the full continuum equations.
The discrete nature of the lattice grid can introduce artifacts in the simulation results, such as grid-induced
turbulence or numerical dissipation. These artifacts can affect the accuracy and realism of the simulation, 
especially in regions of high fluid flow gradients. While LGS can be parallelized, simulating large-scale fluid flows
with high resolution can still be computationally expensive. The need to track individual particle motion and
interactions imposes computational overhead, limiting the size and resolution of simulations that can be performed. 
Despite their conceptual simplicity, implementing and tuning LGS can be complex. Designing appropriate collision and
streaming rules to accurately represent fluid behavior requires careful consideration, and achieving stable and
realistic simulations may require extensive parameter tuning. LGS are well-suited for certain types of fluid flow
problems, such as simulating simple flows or studying emergent phenomena. However, they may not be suitable for all
types of fluid flows, particularly those involving complex geometries or boundary conditions.

