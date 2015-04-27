import sys
import matplotlib.pyplot as plotter
import math

# Inputs to the program are the initial positions, velocities, masses and charges of the particles.

# Output is a graphical representation of the trajectories for a given time interval.


class Particles():

	def __init__(self):
		self.K = 8.99 * math.pow(10, 9)
		self.particles = [{'loc':[1,1.1], 'vx':0, 'vy':0, 'mass':1, 'charge':10}, {'loc':[1,1], 'vx':0, 'vy':0, 'mass':1, 'charge':10}, {'loc':[1,0.9], 'vx':0, 'vy':0, 'mass':1, 'charge':10}]
		self.forces = [self.coulomb(self.particles[0], self.particles[1]), self.coulomb(self.particles[1], self.particles[2]), self.coulomb(self.particles[0], self.particles[2])]

	# calculate distance between two particles using pythagorean theorem
	def distance(self, loc1, loc2):
		x_dist = loc1[0] - loc2[0]
		y_dist = loc1[1] - loc2[1]
		
		return math.sqrt(((math.pow(x_dist, 2)) + (math.pow(y_dist, 2))))


	# define interaction of any two particles
	def coulomb(self, a, b):
		# Kq1q2/r^2
		numerator = self.K * a['charge'] * b['charge']
		denominator = math.pow(self.distance(a['loc'], b['loc']), 2)
		force = (numerator/denominator)

		if ((a['charge'] > 0) & (b['charge'] > 0)) | (a['charge'] < 0 & b['charge'] < 0): # attractive or repulsive forces
			force *= -1
		return force


	# calculates the total force on a particle
	def combine(self, f2, f3, curr, loc2, loc3):
		curr_x = curr[0]
		x2 = loc2[0]
		x3 = loc3[0]
		curr_y = curr[1]
		y2 = loc2[1]
		y3 = loc3[0]
		r2 = self.distance(curr, loc2)
		r3 = self.distance(curr, loc3)

		f1x = f2 * ((curr_x - x2)/r2) # f1cos(theta)
		f2x = f3 * ((curr_x - x3)/r3) # f2cos(theta)
		f1y = f2 * ((curr_y - y2)/r2) # f1sin(theta)
		f2y = f3 * ((curr_y - y3)/r3) # f2sin(theta)

		f_totx = f1x + f2x
		f_toty = f1y + f2y
		return math.sqrt(math.pow(f_totx, 2) + math.pow(f_toty, 2)) # pythagorean theorem

	# calculate vfinal from given information
	def calc_vfinal(self, vi, f1, f2, curr, loc2, loc3, t, p):
		m = p['mass']
		f_tot = self.combine(f1, f2, curr, loc2, loc3)
		a = f_tot/m # acceleration
		
		return math.sqrt((math.pow(vi, 2) + (2*a*t))) # kinmatics


	# calculate x position after some time period
	def calc_xfinal(self, t, index):
		p = self.particles[index]
		f1 = self.forces[index]
		if index+1 < len(self.particles):
			p1 = self.particles[index + 1]
		else:
			p1 = self.particles[0]
		f2 = self.forces[index - 1]
		p2 = self.particles[index - 1]

		v_final = self.calc_vfinal(p['vx'], f1, f2, p['loc'], p1['loc'], p2['loc'], t, p)
		# vi and vf should account for negative or positive displacement
		displacement = ((p['vx'] + v_final)/2) * t # kinematics
		
		return (p['loc'][0] + displacement)


	# calculate y position after some time period
	def calc_yfinal(self, t, index):
		p = self.particles[index]
		f1 = self.forces[index]
		if index+1 < len(self.particles):
			p1 = self.particles[index + 1]
		else:
			p1 = self.particles[0]
		f2 = self.forces[index - 1]
		p2 = self.particles[index - 1]

		v_final = self.calc_vfinal(p['vy'], f1, f2, p['loc'], p1['loc'], p2['loc'], t, p)		
		# vi and vf should account for negative or positive displacement
		displacement = ((p['vy'] + v_final)/2) * t # kinematics
		
		return (p['loc'][1] + displacement)


	# finds lowest entry in position vector
	def find_min(self, pos):
		min_val = float('Inf')
		for arr in pos:
			for val in arr:
				if val < min_val:
					min_val = val

		return min_val


	# finds highest entry in position vector
	def find_max(self, pos):
		max_val = float('-Inf')
		for arr in pos:
			for val in arr:
				if val > max_val:
					max_val = val

		return max_val


	# takes arrays of x and y positions and graphs them on the same axes
	def display(self, x_pos, y_pos):
		xmin = self.find_min(x_pos)
		xmax = self.find_max(x_pos)
		ymin = self.find_min(y_pos)
		ymax = self.find_max(y_pos)


		plotter.plot(x_pos[0], y_pos[0], 'r-')
		#plotter.axis(xmin, xmax, ymin, ymax)
		plotter.xlabel('X Positions Over Time')
		plotter.ylabel('Y Positions Over Time')
		plotter.show()



	# wrapper function for particle interactions
	def simulate(self, particles, t):
		# initial position array, halfway point array, final position array
		# easily modified to show more points along trajectory 
		self.x_pos = [[self.particles[0]['loc'][0], self.calc_xfinal((t/2), 0), self.calc_xfinal(t, 0)], [self.particles[1]['loc'][0], self.calc_xfinal((t/2), 1), self.calc_xfinal(t, 1)], [self.particles[2]['loc'][0], self.calc_xfinal((t/2), 1), self.calc_xfinal(t, 2)]]
		self.y_pos = [[self.particles[0]['loc'][1], self.calc_yfinal((t/2), 0), self.calc_yfinal(t, 0)], [self.particles[1]['loc'][1], self.calc_yfinal((t/2), 1), self.calc_yfinal(t, 1)], [self.particles[2]['loc'][1], self.calc_yfinal((t/2), 1), self.calc_yfinal(t, 2)]]
		
		self.display(self.x_pos, self.y_pos)


particles = Particles()
particles.simulate(particles, 10)



