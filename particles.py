import sys
import matplotlib.pyplot as plotter
import math

# Inputs to the program are the initial positions, velocities, masses and charges of the particles.

# Output is a graphical representation of the trajectories for a given time interval.


class Particles():

	def __init__(self):
		self.K = 8.99 * math.pow(10, 9)
		self.particles = [{'loc':[0,.1], 'vx':0, 'vy':0, 'mass':1, 'charge':.01}, {'loc':[0,0], 'vx':0, 'vy':0, 'mass':1, 'charge':.01}, {'loc':[.1,0], 'vx':0, 'vy':0, 'mass':1, 'charge':.01}]
		self.forces = self.calc_forces()


	# calculate distance between two particles using pythagorean theorem
	def distance(self, loc1, loc2):
		x_dist = loc1[0] - loc2[0]
		y_dist = loc1[1] - loc2[1]
		
		return math.sqrt(((math.pow(x_dist, 2)) + (math.pow(y_dist, 2))))


	def calc_forces(self):
		return [self.coulomb(0, 1), self.coulomb(0, 2), self.coulomb(1, 0), self.coulomb(1, 2), self.coulomb(2, 0), self.coulomb(2, 1)]


	# define interaction of any two particles: Kq1q2/r^2
	def coulomb(self, a, b):
		p1 = self.particles[a]
		p2 = self.particles[b]
		numerator = self.K * p1['charge'] * p2['charge']
		denominator = math.pow(self.distance(p1['loc'], p2['loc']), 2)
		f = (numerator/denominator)
		dist = self.distance(p1['loc'], p2['loc'])
		fx = f * ((p1['loc'][0] - p2['loc'][0])/dist) # f1cos(theta) --> opp/hyp
		fy = f * ((p1['loc'][1] - p2['loc'][1])/dist) # f2cos(theta) --> opp/hyp

		if ((p1['charge'] > 0.0) & (p2['charge'] > 0.0)) | ((p1['charge'] < 0.0) & (p2['charge'] < 0.0)): # repulsive forces
			if p1['loc'][0] < p2['loc'][0]:
				fx *= -1 # left particle pushed left
			if p1['loc'][1] < p2['loc'][1]:
				fy *= -1 # lower particle pushed down
		else: # attractive forces
			if p1['loc'][0] > p2['loc'][0]:
				fx *= -1 # right particle pulled left
			if p1['loc'][1] > p2['loc'][1]:
				fy *= -1 # lower particle pulled up

		force = [fx, fy]
		return force


	# calculate vfinal from given information -- x or y direction
	def calc_vfinal(self, vi, a, t):	
		vf = (math.pow(vi, 2) + (2*a*t))
		return math.sqrt(vf) # kinematics


	# combine two forces acting on a particle in x or y direction
	def combine(self, f1, f2, coord):
		return f1[coord] + f2[coord]


	# determine acceleration in x or y direction
	# (simulate based on constant-acceleration even though that's not what really happens)
	def calc_acceleration(self, index, coord):
		curr = self.particles[index]
		i = index*2 # eg index 0 in particles array needs indices 0 and 1 in forces array
		f1 = self.forces[i]
		f2 = self.forces[i + 1]
		if (index + 1) < len(self.particles):
			loc1 = self.particles[index + 1]['loc']
		else:
			loc1 = self.particles[0]['loc'] # roll over to initial force
		loc2 = self.particles[index - 1]['loc'] # python automatically wraps around

		# total force in x or y direction
		f_tot = self.combine(f1, f2, coord)

		a = f_tot/curr['mass']
		return a


	# calculate x position after some time period
	def calc_pos_final(self, t, index, coord):
		p = self.particles[index]
		a = self.calc_acceleration(index, coord)

		# xf = xi + vt + 1/2at^2
		if coord == 0:
			final = p['loc'][coord] + p['vx'] + (1/2)*a*(math.pow(t, 2))
			p['vx'] = self.calc_vfinal(p['vx'], a, t) # update for next iteration

		else:
			final = p['loc'][coord] + p['vy'] + (1/2)*a*(math.pow(t, 2))
			p['vy'] = self.calc_vfinal(p['vy'], a, t) # update for next iteration

		p['loc'][coord] = final # update loc for next iteration
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
		# SOMETHING IS WRONG WITH THE ITERATION/PLOTTING
		for increment in range(1, 100): # incremental progression through the time period
			locations = []
			count = 0
			for particle in self.particles: # update all particles for each time increment
				locations.append(particle['loc'][coord])
				position = self.calc_pos_final(t/increment, count, coord)
				locations.append(position)
				count += 1 # indicates index of particle in self.particles array
			self.forces = self.calc_forces()
			wrapper.append(locations)
		return wrapper


	# wrapper function for particle interactions
	def simulate(self, t):
		x_pos = self.generate_pos(t, 0) # 0 indicates X in X,Y pair
		y_pos = self.generate_pos(t, 1) # 1 indicates Y in X,Y pair
		for elem in x_pos:
			print elem
			print '\n'
		for elem in y_pos:
			print elem
			print '\n'

		self.display(x_pos, y_pos)


particles = Particles()
particles.simulate(.1)