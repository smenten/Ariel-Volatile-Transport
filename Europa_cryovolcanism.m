


%% Friction equation (temperature change) **WORKING**

k = 2.3  ;              % thermal conductivity of water ice 
rho_i = 917 ;           % density of ice in kg/m^3
C = 2050  ;             % specific heat capacity of water ice
u = 3.171e-9;              % slab velocity (m/s)
phi = 0;                % angle between direction of velocity and slab subduction
theta = 4.4 ;           % angle of subduction 
miu = 0.55 ;            % coefficient of friction (ice sliding on ice)
d_v = 6900 ;            % depth to melt
g = 1.315 ;             % gravity of Europa in m/s^2 
kappa = k/(rho_i*C) ;       % thermal diffusivity of water ice


tau = miu*rho_i*g*d_v     ;                 % shear stress

a = ((2*tau)/k) ;
b = (kappa*d_v*u) ;
c = (pi*cosd(phi)*sind(theta)) ;

d_T = a*sqrt(b/c) ;
