Calculates an approximation of the movement of three charged particles due to Coulomb's Law during some time interval,
restricted to the x-y plane due to limitations in graphing library. The file particles.py includes force, acceleration, and velocity calculations and plotting in two dimensions, and the file 3d_particles.py includes the same but without plotting in three dimensions.

To run: enter the particles' initial positions, velocity in the X and Y (and Z if in three-dimensional mode) directions, mass, and charge on line 12 (either version), saves, then runs "python particles.py" or "python 3d_particles.py". Change time interval on line 195 (particles) or 204 (3d_particles); it is the argument to the function simulate(). Change increments through the time interval on line 145 (particles) or 152 (3d_particles); it is the upper bound of the xrange() call.

Initial conditions:

particles:
{'loc':[1,1, 1], 'vx':0, 'vy':0, 'vz':0, 'mass':1, 'charge':1}
{'loc':[0,0, 0], 'vx':0, 'vy':0, 'vz':0, 'mass':1, 'charge':1}
{'loc':[.1,0, 0], 'vx':0, 'vy':0, 'vz':0, 'mass':1, 'charge':1}

t: 10

increments: 100
