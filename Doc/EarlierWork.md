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

