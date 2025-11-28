function output = getDesignfromWindspeed(B, radiusLimits, lambda, chordProfile, AoAProfile, noElements, V, kinematicViscosity, convergenceTolerance, iterationLimit, plot)

B = %Number of Blades
radiusLimits = %Start and End radius for blade in meters
lambda = %Tip Speed ratio
chordProfile = %Chord Lengths at root to tip
AoAProfile = %Angle of Attacks
noElements = %Number of Blade Elements
V = %Velocity of Free Stream
kinematicViscosity = %Kinematic Velocity of the Air
convergenceTolerance = %the maximum allowable fluctuation in a and aPrime to conclude convergence
iterationLimit = %The iteration Limit
plot = %0 or 1 plot or not


CP = BEMT_Demo(B, radiusLimits, lambda, chordProfile, AoAProfile, noElements, V, kinematicViscosity, convergenceTolerance, iterationLimit, plot)

output = CP