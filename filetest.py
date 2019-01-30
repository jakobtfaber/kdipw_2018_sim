import sys
import os

def main():

	params = [] # List of parameters inputted


	# Check if an input file was provided
	if os.isatty(0):
		print("Usage: python ./filetest.py < params.txt")
		exit(2)

	# Read inputfile and store params in a list
	for line in sys.stdin:
		for c in line:
			# Skip lines that start with #
			if c == '#':
				break
			# Read number from line
			elif c.isdigit():
				params.append(eval(line))
				break

	if len(params) != 4:
		print("Invalid parameter file")
		exit(3)

	# Assign params to appropriate values
	vel = params[0] # Meters per second
	nscreen = params[1]
	trials = params[2] # Number of raypaths
	sigma0 = params[3] # Overall angular scale (milliarcseconds = mas)

	
	print(vel, nscreen, trials, sigma0)



if __name__ == "__main__":
	main()