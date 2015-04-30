import sys
import matplotlib.pyplot as plotter
import math

# Inputs to the program are the initial positions, velocities, masses and charges of the particles.

# Output is a graphical representation of the trajectories for a given time interval.


class Particles():

	def __init__(self):
		self.K = 8.99 * math.pow(10, 9)
		self.particles = [{'loc':[.1,.1], 'vx':0, 'vy':0, 'mass':1, 'charge':.01}, {'loc':[0,0], 'vx':0, 'vy':0, 'mass':1, 'charge':.01}, {'loc':[.1,0], 'vx':0, 'vy':0, 'mass':1, 'charge':.01}]
		self.forces = [self.coulomb(self.particles[0], self.particles[1]), self.coulomb(self.particles[1], self.particles[2]), self.coulomb(self.particles[0], self.particles[2])]


	# calculate distance between two particles using pythagorean theorem
	def distance(self, loc1, loc2):
		x_dist = loc1[0] - loc2[0]
		y_dist = loc1[1] - loc2[1]
		
		return math.sqrt(((math.pow(x_dist, 2)) + (math.pow(y_dist, 2))))


	# define interaction of any two particles: Kq1q2/r^2
	def coulomb(self, a, b):
		numerator = self.K * a['charge'] * b['charge']
		denominator = math.pow(self.distance(a['loc'], b['loc']), 2)
		force = (numerator/denominator)

		if ((a['charge'] > 0.0) & (b['charge'] > 0.0)) | ((a['charge'] < 0.0) & (b['charge'] < 0.0)): # attractive or repulsive forces
			force *= -1
		return force


	# calculates the total force on a particle in either the x or y coordinate
	def combine(self, f2, f3, curr_loc, loc2, loc3, coord):
		curr = curr_loc[coord]
		p1 = loc2[coord]
		p2 = loc3[coord]
		r1 = self.distance(curr_loc, loc2)
		r2 = self.distance(curr_loc, loc3)

		f1 = f2 * ((curr - p1)/r1) # f1cos(theta)
		f2 = f3 * ((curr - p2)/r2) # f2cos(theta)

		f_tot = f1 + f2
		return f_tot # pythagorean theorem


	# calculate vfinal from given information
	def calc_vfinal(self, vi, a, t):	
		vf = (math.pow(vi, 2) + (2*a*t))
		if vf < 0: # guard against negative velocity -- can be negative, but can't calculate this way
			vf *= -1
			vf = math.sqrt(vf)
			vf *= -1
			return vf # ewww this style is gross but whatever
		else:
			return math.sqrt((math.pow(vi, 2) + (2*a*t))) # kinematics


	# determine acceleration from given information
	# (simulate based on constant-acceleration even though that's not what really happens)
	def calc_acceleration(self, index, coord):
		curr = self.particles[index]
		f1 = self.forces[index]
		if index + 1 < len(self.particles):
			loc1 = self.particles[index + 1]['loc']
		else:
			loc1 = self.particles[0]['loc']
		f2 = self.forces[index - 1]
		loc2 = self.particles[index - 1]['loc']
		f_tot = self.combine(f1, f2, curr['loc'], loc1, loc2, coord)
		a = f_tot/curr['mass']
		return a


	# calculate x position after some time period
	def calc_pos_final(self, t, index, coord):
		p = self.particles[index]
		a = self.calc_acceleration(index, coord)

		# x = x + vt + 1/2at^2
		if coord == 0:
			final = p['loc'][coord] + p['vx'] + (1/2)*a*(math.pow(t, 2))
			p['vx'] = self.calc_vfinal(p['vx'], a, t)

		else:
			final = p['loc'][coord] + p['vy'] + (1/2)*a*(math.pow(t, 2))
			p['vy'] = self.calc_vfinal(p['vy'], a, t)

		# update vx and loc for next iteration
		p['loc'][coord] = final
		return final


# -----------------------------------------------------
	# takes arrays of x and y positions and graphs them on the same axes
	def display(self, x_pos, y_pos):
		plotter.plot(x_pos[0], y_pos[0], 'r-', x_pos[1], y_pos[1], 'b-', x_pos[2], y_pos[2], 'g-')
		plotter.xlabel('X Positions Over Time')
		plotter.ylabel('Y Positions Over Time')
		plotter.show()


	# increments through the time interval and calculates the position of the particles after that time
	def generate_pos(self, t, coord):
		wrapper = []
		count = 0
		for particle in self.particles: # iterate through indexes of particles
			locations = []
			locations.append(particle['loc'][coord])
			for increment in range(1, 100): # incremental progression through the time period
				if coord == 0:
					position = self.calc_pos_final(t/increment, count, coord)
				else:
					position = self.calc_pos_final(t/increment, count, coord)
				locations.append(position)
			wrapper.append(locations)
			count += 1
		return wrapper


	# wrapper function for particle interactions
	def simulate(self, t):
		x_pos = self.generate_pos(t, 0)
		y_pos = self.generate_pos(t, 1)

		self.display(x_pos, y_pos)


particles = Particles()
particles.simulate(10)