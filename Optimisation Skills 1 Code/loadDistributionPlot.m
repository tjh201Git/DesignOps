% Load distribution plot

N = 100;
R = 1.5;

r = linspace(0,R,N);

load = 50*r.*(1-(1./exp(-50*(r-R))));

normalised_radial_position = r / R;


% Plot the load distribution
fig = figure;
theme(fig, 'light')
plot(normalised_radial_position, load);
xlabel('Normalised Radial Position, r/R');
ylabel('Load Distribution, N/m');
% title('Load Distribution vs Radius');
% grid on;