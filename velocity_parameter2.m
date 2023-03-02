function [ integral ] = velocity_parameter2( v0, T, int)
%This calculates the integral over the velocity distribution function for a
%given velocity v0.  So if v0 is very high, the answer will be close to 1
%because almost all velocities in the distribution are < v0.  This is for a
%molecule ballistically moving across a planet's surface

%Stefan-Boltzmann constant, m^2 kg s^-2 K^-1
k = 1.38e-23;


%Mass of molecule, kg
m = 44.01*.001*(1/6.02e23);  % For CO2
%m = 28*.001*(1/6.02e23);   %For N2
%m = 18.02*.001*(1/6.02e23);   %for H20
%m = 2.66e-26;                  %for CH4


a = m/(2*k*T);
c = 4*pi*((m/(2*pi*k*T))^1.5);

integral = c * ((.25 * ((pi/(a^3))^0.5) * erf(v0 * a^.5))  -  (v0/(2*a))*exp(-a*v0^2)) - int;



end