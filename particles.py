import sys
import matplotlib.pyplot as plotter
import math

# Inputs to the program are the initial positions, velocities, masses and charges of the particles.

# Output is a graphical representation of the trajectories for a given time interval.


class Particles():

	def __init__(self):
		# CHECKED
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


	# calculates the total force on a particle
	def combine(self, f2, f3, curr, loc2, loc3):
		#print '>>>>>>>>>>>>>>>>>>>>>  ' + str(curr)
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
		#print '>>>>>>>>>>>>>>>>>>>>>  x force: ' + str(f1x) + ' ' + str(f2x) + ' ' + str(f_totx)
		#print '>>>>>>>>>>>>>>>>>>>>>  y force: ' + str(f1y) + ' ' + str(f2y) + ' ' + str(f_toty)
		total = math.sqrt(math.pow(f_totx, 2) + math.pow(f_toty, 2))
		#print '>>>>>>>>>>>>>>>>>>>>>  net force: ' + str(total)
		return total # pythagorean theorem

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
		if index + 1 < len(self.particles):
			p1 = self.particles[index + 1]
		else:
			p1 = self.particles[0]
		f2 = self.forces[index - 1]
		p2 = self.particles[index - 1] # this wraps around so no checks are needed (arr[-1] is well defined)

		v_final = self.calc_vfinal(p['vx'], f1, f2, p['loc'], p1['loc'], p2['loc'], t, p)
		# vi and vf should account for negative or positive displacement
		displacement = ((p['vx'] + v_final)/2) * t # kinematics
		p['loc'][0] += displacement
		return p['loc'][0]


	# calculate y position after some time period
	def calc_yfinal(self, t, index):
		p = self.particles[index]
		f1 = self.forces[index]
		if index + 1 < len(self.particles):
			p1 = self.particles[index + 1]
		else:
			p1 = self.particles[0]
		f2 = self.forces[index - 1]
		p2 = self.particles[index - 1]

		v_final = self.calc_vfinal(p['vy'], f1, f2, p['loc'], p1['loc'], p2['loc'], t, p)	

		# vi and vf should account for negative or positive displacement
		displacement = ((p['vy'] + v_final)/2) * t # kinematics)
		p['loc'][1] += displacement # store to use for new distance calculation in next iteration
		return p['loc'][1]


# -----------------------------------------------------


	# takes arrays of x and y positions and graphs them on the same axes
	def display(self, x_pos, y_pos):
		print 'XPOSITIONS'
		for array in x_pos:
			print array
			print '\n'
		print 'YPOSITIONS'
		for array in y_pos:
			print array
			print '\n'
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
					position = self.calc_xfinal(t/increment, count)
				else:
					position = self.calc_yfinal(t/increment, count)
				locations.append(position)
			wrapper.append(locations)
			count += 1
		return wrapper


	# wrapper function for particle interactions
	def simulate(self, t):
		# initial position array, halfway point array, final position array
		# easily modified to show more points along trajectory 
		#print 'x1 i ' + str(self.particles[0]['loc'][0])
		#print 'x1 t/2 ' + str(self.calc_xfinal((t/2), 0))
		#print 'x1 t ' + str(self.calc_xfinal(t, 0))
		#rint 'x2 i ' + str(self.particles[1]['loc'][0])
		#print 'x2 t/2 ' + str(self.calc_xfinal((t/2), 1))
		#print 'x2 t ' + str(self.calc_xfinal(t, 1))
		#print 'x3 i ' + str(self.particles[2]['loc'][0])
		#print 'x3 t/2 ' + str(self.calc_xfinal((t/2), 1))
		#print 'x3 t ' + str(self.calc_xfinal(t, 2)) + '\n'

		#print 'y1 i ' + str(self.particles[0]['loc'][1])
		#print 'y1 t/2 ' + str(self.calc_yfinal((t/2), 0))
		#print 'y1 t ' + str(self.calc_yfinal(t, 0))
		#print 'y2 i ' + str(self.particles[1]['loc'][1])
		#print 'y2 t/2 ' + str(self.calc_yfinal((t/2), 1))
		#print 'y2 t ' + str(self.calc_yfinal(t, 1))
		#print 'y3 i ' + str(self.particles[2]['loc'][1])
		#print 'y3 t/2 ' + str(self.calc_yfinal((t/2), 1))
		#print 'y3 t ' + str(self.calc_yfinal(t, 2)) + '\n'
		#print self.particles[0]
		#print self.particles[1]
		#print self.particles[2]
		#print '\n'

		#print '1 and 2: ' + str(self.forces[0])
		#print '2 and 3: ' + str(self.forces[1])
		#print '1 and 3: ' + str(self.forces[2])
		x_pos = self.generate_pos(t, 0)
		y_pos = self.generate_pos(t, 1)

		#self.x_pos = [[self.particles[0]['loc'][0], self.calc_xfinal((t/2), 0), self.calc_xfinal(t, 0)], [self.particles[1]['loc'][0], self.calc_xfinal((t/2), 1), self.calc_xfinal(t, 1)], [self.particles[2]['loc'][0], self.calc_xfinal((t/2), 1), self.calc_xfinal(t, 2)]]
		#self.y_pos = [[self.particles[0]['loc'][1], self.calc_yfinal((t/2), 0), self.calc_yfinal(t, 0)], [self.particles[1]['loc'][1], self.calc_yfinal((t/2), 1), self.calc_yfinal(t, 1)], [self.particles[2]['loc'][1], self.calc_yfinal((t/2), 1), self.calc_yfinal(t, 2)]]
		
		self.display(x_pos, y_pos)


particles = Particles()
particles.simulate(.01)



