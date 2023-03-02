function [ Tsurf ] = thermal_ariel_f(latitude)

%% SEMI-IMPLICIT 1D THERMAL CONDUCTION CODE
% Features: DIFFERENT LAYER PROPERTIES, CO2 FROST, SCATTERED RADIATION, SLOPES AND HORIZON EFFECTS (5/16/2016)
% Author: Ali Bramson



%% Input parameters
a = 19.8;          % Orbital semi-major axis in AU
ecc = 0.0012;        % Orbital eccentricity
obl = 98.0;             % Obliquity in degrees
Lsp = 96.99;            % Ls of perihelion in degrees
rot = 217728;           % Length of a solar (not sidereal) day in seconds
Q = .003;                % Heat flux
A = 0.23;                 % Bond Albedo
Afrost = 0.7;            % Frost albedo
emis = 0.90;              % Emissivity
emisFrost = 0.7;         % Emissivity of frost
lat = latitude;                % Latitude
k = 0.00017;              % Thermal conductivity (W m^-1 K^-1)
c = 867;             % Specific heat J/(kg K)
rho = 1615;               % Density of subsurface (kg/m3)
Tfrost = 0;            % Frost point in Kelvins
L_CO2 = 5.899e5;         % Latent heat of CO2 frost (J/kg)
slope = 0;               % Slope of terrain in degrees
slope_aspect = 0;        % Direction slope faces
saveFlat = 0;            % Whether to save output of radiation scattered off of 'surroundings' to be used for sloped case
downwellingPerc = 0;  % Percentage of noon-time IR scattered
scatteredVisPerc = 0; % Percentage of scattered visible radiation

saveTsurf = 1;                     % If this is 1, the script saves surface temps to the file given below; set to 0 if you don't want to save output
TsurfSaved = '90LatSurfaceTemp.mat';    % File to store surface temperatures and Ls at the end of the year



%% Numerical properties
dt = 1e4;              % timestep, in seconds
totalRunTime = 5;      % number of years to run model
f = .5;                 % explicit (0) vs implicit (1); 0.5 = Crank-Nicholson

%% Constants
au  = 1.4959787061e11;           % Size of one astronomical unit in meters
gc  = 6.67259e-11;               % Gravitational constant
sm  = 1.9891e30;                 % Solar mass
d2r = pi/180;                    % Convert degrees to radians
r2d = 180/pi;                    % Convert radians to degrees
yr  = 2*pi*sqrt(au^3/(gc*sm));   % Length of one Earth year in seconds
f1au = 1367;                     % Solar Flux at 1AU (W/m^2)
sigmasb = 5.67e-8;               % Stefan Boltzmann constant W m^-2 K^-4

year = 2*pi/sqrt(gc*sm/(a*au)^3);

%% Properties

Kappa = k ./ (rho .* c);
diurnalSkinDepths = sqrt(Kappa*rot/pi); % meters
%  fprintf ('Diurnal thermal skin depth of top layer = %.8f m\n', diurnalSkinDepths);
annualSkinDepths = sqrt(Kappa*year/pi); % meters
%  fprintf ('Annual thermal skin depth of bottom layer = %.8f m\n', annualSkinDepths);

layerGrowth = 1.03;
dailyLayers = 15;
annualLayers = 6;

firstLayerThickness = diurnalSkinDepths / ((1-layerGrowth^dailyLayers)/(1-layerGrowth)); % Thickness of first layer based on diurnal skin depth of surface material
nLayers = ceil(log(1-(1-layerGrowth)*(annualLayers*annualSkinDepths/firstLayerThickness))/log(layerGrowth) ); % Number of subsurface layers based on annual skin depth of deepest layer

dz = (firstLayerThickness * layerGrowth.^[0:nLayers-1])'; % transpose to make column vector

% Change here if want properties to change at a certain depth
k = zeros(nLayers,1) + k;
%k(6:end) = 2;
c = zeros(nLayers,1) + c;
%c(6:end) = 3000;
rho = zeros(nLayers,1) + rho;
%rho(6:end) = 2000;
%dz = (0.005 * 1.03.^[0:nLayers-1])'; % changing dz's

depthsAtMiddleOfLayers = cumsum(dz) - dz/2;
depthsAtLayerBoundaries = cumsum(dz);
depthBottom = sum(dz);

%% Describe orbit of a body
% Get solar distance and subsolar latitude for given time taking into
% account eccentricity, semi-major axis, obliquity and longitude of
% perihelion. Then convert to solar flux.

% Low-res calculation of times at uniformly spaced true anomalies
year = 2*pi/sqrt(gc*sm/(a*au)^3);
n=1e4;
TA = ((0:n)'+0.5)*2*pi/n;        % create evenly spaced true anomalies
TA=TA(1:end-1);   % do because matlab indexes start at 1
EA = acos( (ecc+cos(TA))./(1+ecc*cos(TA)) );   % Calculate eccentric anomalies
MA = EA - ecc*sin(EA);                   % Calculate mean anomalies
t  = MA/sqrt(gc*sm/(a*au)^3);     %Time along Mars orbital path (will be irregularly spaced)
t(find(TA > pi)) = year - t(find(TA > pi));  % correction of 2nd half of the orbit

% High-Res interpolation of true anomalies at uniformly spaced times
% Note: evenly spaced timesteps but doesn't go to the end of a solar
% day- ends at about 1 pm on the last day, rather than midnight
t2 = ((0:(year/dt + 1))' + .5)/(year/dt + 1) * year; % Now use these uniformly spaced times
t2=t2(1:end-1);  % different in matlab than idl- matlab starts at 1 whereas idl indexes start at 0

% Now use these evenly spaced times to calculate the rest
TA2 = interp1(t,TA,t2, 'pchip'); %Interpolate the True anomalies at these times
EA2 = acos( (ecc+cos(TA2))./(1+ecc*cos(TA2)) );  % Calculate the eccentric anomalies at these times
nStepsInYear = size(t2,1);

% Solar distance/declination/longitude and local time calculations
sol_dist = a*(1 - ecc*cos(EA2));               % Solar Distance
lsrad = mod((TA2 + Lsp*d2r), (2*pi));  % in radians for taking sin of
ls = rad2deg(TA2 + Lsp*d2r); % so it wont wrap around in plots
lsWrapped = rad2deg(mod((TA2 + Lsp*d2r), (2*pi)));   % Solar longitude in degrees - wraps around after 360 degrees to make plots better
sin_dec = sin(obl*d2r)*sin(lsrad);               % Sine of solar declination

cos_dec = sqrt(1 - sin_dec.^2);
hr = mod((t2/rot * 2*pi), (2*pi)) - pi;  % Local mean true solar time in hour angle (-pi to pi)- cant be used to compare to ephemeris time because fixed rotation...calculated as if was a mean solar
ltst = hr*r2d/15 + 12;  % between 0 and 24

% Incidence angles and fluxes
cosi = ( sin(lat*d2r)*sin_dec + cos(lat*d2r)*cos_dec.*cos(hr) );

% Slopes
sini   = sqrt(1-cosi.^2);
az  = acos( (sin_dec - sin(lat*d2r).*cosi) ./ (cos(lat*d2r).*sini) );
az(find(ltst > 12)) = 2*pi - az(find(ltst > 12)); % flip around since arccos degenerate

cosi_slope = cosi.*cos(slope*d2r) + sini.*sin(slope*d2r).*cos(slope_aspect*d2r - az);
cosi_slope(find(cosi_slope < 0)) = 0;  % Set any value less than 0 to 0...shadowed
cosi_slope(find(cosi < 0)) = 0; 

% Solar Flux
sf   = f1au./(sol_dist.^2) .* cosi_slope;

if slope == 0
    cosi(find(cosi < 0)) = 0; 
    sf   = f1au./(sol_dist.^2) .* cosi;
end

annual_sf = sum(sf*dt);   % Total annual energy


% fprintf('Minimum Solar Flux: %8.4f [W/m^2]\n', min(sf));
% fprintf('Maximum Solar Flux: %8.4f [W/m^2]\n', max(sf));
% fprintf('Mean Solar Flux: %8.4f [W/m^2]\n', mean(sf));
% 
% fprintf ('Annual Solar Flux = %.6e [W/m^2] \n', annual_sf);

% Daily noontime flux and the value of this to be used for downwelling IR
% radiation
% Will later need to change when including effects of sloped terrain
sf_noon = f1au./(sol_dist.^2) .* ( sin(lat*d2r)*sin_dec + cos(lat*d2r)*cos_dec);
IRdown = downwellingPerc .* sf_noon;

% Scattered Visible Light at each timestep (will get 0 scattered vis light
% at night)- this assumes light scatters isotropically even though in
% actuality more light is scattered near the disk of the sun in the sky
% Will later need to change when including effects of sloped terrain
visScattered = scatteredVisPerc .* sf;

if slope ~= 0
    load flatSaved.mat;
    sky = cos(slope*d2r/2)^2;
else
    flatVis = zeros(nStepsInYear,1);
    flatIR = zeros(nStepsInYear,1);
    sky = 1;
end

%% Get ready to run model
% First column of these arrays will be "duds" to kick off the year
Temps = zeros(nLayers, nStepsInYear);
Tsurf = zeros(1,nStepsInYear);

oldTemps = zeros(nLayers, 1);

frostMass = 0;
frostMasses = zeros(1,nStepsInYear);

Tref = 55;
oldTemps(:) = Tref;
%fprintf ('Initial Surface Temperature = Tref = %2.2f K\n', Tref);

alpha_u = (2*k.*circshift(k,[1 0])./(k.*circshift(dz,[1 0]) + circshift(k,[1 0]).*dz)).*(dt./(rho.*c.*dz));
alpha_u(1) = 0;
alpha_d = (2*k.*circshift(k,[-1 0])./(k.*circshift(dz,[-1 0]) + circshift(k,[-1 0]).*dz)).*(dt./(rho.*c.*dz));
alpha_d(end) = 0;

dia_e = zeros(nLayers, 1);
dia_e = 1 - (1-f)*alpha_u - (1-f)*alpha_d;
dia_e(nLayers) = 1 - (1-f)*alpha_u(end);
boundary = zeros(nLayers, 1); % column vector
boundary(end) = dt*Q/(rho(end)*c(end)*dz(end));
k1_e = (1-f)*alpha_d;
k3_e = (1-f)*alpha_u;

dia_i = zeros(nLayers, 1);
dia_i = 1 + f*alpha_u + f*alpha_d;
dia_i(nLayers) = 1+f*alpha_u(end);
Amatrix_i = zeros(nLayers, nLayers);
Amatrix_i = diag(dia_i) + diag(-f*alpha_u(2:end),-1) + diag(-f*alpha_d(1:end-1),1);
Amatrix_i = sparse(Amatrix_i);

beta = k(1)*dt/(rho(1)*c(1)*dz(1)*dz(1));
        
sf_i = [sf; sf(1)];
IRdown_i = [IRdown; IRdown(1)];
visScattered_i = [visScattered; visScattered(1)];
T_e = zeros(nLayers,1);

Fin = (sf + visScattered.*sky + flatVis.*(1-sky)).*(1-A) + (IRdown.*sky + flatIR.*(1-sky)).*emis;
Fin_frost = (sf + visScattered.*sky + flatVis.*(1-sky)).*(1-Afrost) + (IRdown.*sky + flatIR.*(1-sky)).*emisFrost;
Fin_i = (circshift(sf, [-1 0]) + circshift(visScattered, [-1 0]).*sky + circshift(flatVis, [-1 0]).*(1-sky)).*(1-A) + (circshift(IRdown, [-1 0]).*sky + circshift(flatIR, [-1 0]).*(1-sky)).*emis;
%Fin_frost_e = (circshift(sf, [-1 0]) + circshift(visScattered, [-1 0]) + circshift(flatVis, [-1 0]).*(1-sky)).*(1-Afrost) + (circshift(IRdown, [-1 0]) + circshift(flatIR, [-1 0]).*(1-sky)).*emisFrost;

% Calculate a and b's for surface temperature calculation
aa = (dz(1)/(2*k(1)))*(Fin(1) + 3*emis*sigmasb*Tref^4)/(1+(4*emis*sigmasb*Tref^3*dz(1)/(2*k(1))));
b = 1/(1+(4*emis*sigmasb*Tref^3*dz(1)/(2*k(1))));
Tsurf(1) = aa+b*Tref;

% For frost mass calculations
gamma_frost = -(1/L_CO2)*(2*k(1)*dt/dz(1));
theta_frost = (dt/L_CO2)*(2*k(1)*Tfrost/dz(1) - Fin_frost + emisFrost*sigmasb*Tfrost^4);


for yr=1:totalRunTime
    
    for n = 1:nStepsInYear
        
        if frostMass == 0
            
            b = 1/(1+(4*emis*sigmasb*Tref^3*dz(1)/(2*k(1))));
            a_e = (dz(1)/(2*k(1)))*(Fin(n) + 3*emis*sigmasb*Tref^4)/(1+(4*emis*sigmasb*Tref^3*dz(1)/(2*k(1))));
            a_i = (dz(1)/(2*k(1)))*(Fin_i(n) + 3*emis*sigmasb*Tref^4)/(1+(4*emis*sigmasb*Tref^3*dz(1)/(2*k(1))));
            boundary(1) = 2*beta*((1-f)*a_e + f*a_i);
      
            % Explicit Part
            dia_e(1) = 1 - (1-f)*(alpha_d(1)+(2-2*b)*beta);
            T_e = k3_e.*circshift(oldTemps,[1 0]) + k1_e.*circshift(oldTemps,[-1 0]) + dia_e.*oldTemps + boundary;

            % Implicit Part
            Amatrix_i(1,1) = 1 + f*(alpha_d(1)+(2-2*b)*beta);
            Temps(:,n) = Amatrix_i\T_e;
        
            Tsurf(n) = a_i + b*Temps(1,n);  % Uses implicit a with new T calculation- instantanous balance
            oldTemps(:) = Temps(:,n);
            Tref = Tsurf(n);
            
            frostMass = 0;
        
            if Tsurf(n) <= Tfrost
                
                frostMass = 0.25*(Tfrost-Tsurf(n))*rho(1)*c(1)*dz(1)/L_CO2;
                Tsurf(n) = Tfrost;
                Tref = Tsurf(n);
                
            end
            
        elseif frostMass > 0
            
            % In this case, a = Tfrost and b = 0
            boundary(1) = 2*beta*Tfrost;
      
            % Explicit Part
            dia_e(1) = 1 - (1-f)*(alpha_d(1)+2*beta);
            T_e = k3_e.*circshift(oldTemps,[1 0]) + k1_e.*circshift(oldTemps,[-1 0]) + dia_e.*oldTemps + boundary;

            % Implicit Part
            Amatrix_i(1,1) = 1 + f*(alpha_d(1)+2*beta);
            Temps(:,n) = Amatrix_i\T_e;
        
            Tsurf(n) = Tfrost;
            oldTemps(:) = Temps(:,n);
            Tref = Tsurf(n); 
            
            frostMass = frostMass + gamma_frost*Temps(1,n) + theta_frost(n);
            
            if frostMass < 0
                
                Tsurf(n) = Tfrost - (4*frostMass*L_CO2/(rho(1)*c(1)*dz(1)));
                Tref = Tsurf(n);
                frostMass = 0;
                
            end
            

        end

        frostMasses(n) = frostMass;
    end
    
%     fprintf('You''re %2.0f / %2.0f of the way there!\n ', yr, totalRunTime);
    
    if yr == 5
       % fprintf('Windup done!');
        windupTemp = mean(Tsurf);
        oldTemps(:) = windupTemp;
    end
    
end

% Model done running- print output
rho_CO2ice = 1600;
equivalentThicknesses = frostMasses./rho_CO2ice;
fprintf('Max CO2 frost thickness: %8.4f m\n', max(equivalentThicknesses));

fprintf('Minimum Surface Temp: %8.4f K\n', min(Tsurf));
fprintf('Maximum Surface Temp: %8.4f K\n', max(Tsurf));
fprintf('Mean Surface Temp: %8.4f K\n', mean(Tsurf));

% Find min, max and average temperatures at each depth over the last year
minT =  min(Temps(:, :), [], 2);
maxT =  max(Temps(:, :), [], 2);
averageTempDepth = mean(Temps(:, :), 2);
averageSurfTemp = mean(Temps);

T = mean(Tsurf);

%% Plot surface temperatures and CO2 frost mass throughout the year
fonts = 18;
fontn = 'Helvetica';
fontw = 'bold';

figure(1);
clf;

set(gca, 'fontname', fontn, 'fontsize', fonts, 'fontweight', fontw, 'linewidth', 1);
set(gcf, 'Color', 'w');
set(gcf,'Units','centimeters');
set(gcf,'outerposition',[2 20 25 17]); % left bottom width height

whereCrossOver360to0 = find(floor(lsWrapped)==0,1); % First index after Ls wrapped from 360 back to 0
subplot(2,1,1)
plot(lsWrapped(1:whereCrossOver360to0-1), Tsurf(1:whereCrossOver360to0-1), 'b');
hold on
plot(lsWrapped(whereCrossOver360to0:end), Tsurf(whereCrossOver360to0:end), 'b');
ylabel('Temperature (K)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');
xlim([0 360]);
% % subplot(2,1,2)
% % plot(lsWrapped(1:whereCrossOver360to0-1), frostMasses(1:whereCrossOver360to0-1), 'b');
% % hold on
% % plot(lsWrapped(whereCrossOver360to0:end), frostMasses(whereCrossOver360to0:end), 'b');
% % xlim([0 360]);
% % str = sprintf('Solar Longitude (\\circ)');
% % xlabel(str, 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal'); % In solar longitude
% % ylabel('CO_2 Frost Mass (kg/m^2)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');
% 
% %% Plot diurnal curves and min and max at each depth
% figure(2);
% clf;
% 
% set(gca, 'fontname', fontn, 'fontsize', fonts, 'fontweight', fontw, 'linewidth', 1);
% set(gcf, 'Color', 'w');
% set(gcf,'Units','centimeters');
% set(gcf,'outerposition',[2 20 25 17]); % left bottom width height
% 
% p2(1) = plot(averageTempDepth, -depthsAtMiddleOfLayers, 'k');
% set(gca,'FontSize',fonts-1);
% xlabel('Temperature (K)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');
% ylabel('Depths (m)', 'fontname', fontn, 'fontsize', fonts, 'fontweight', 'normal');
% str = sprintf('Temperature Curves with Depth');
% title(str, 'fontname', fontn, 'fontsize', fonts, 'fontweight', fontw);
% 
% hold on;
% for ii = 1:floor(nStepsInYear/10):nStepsInYear
%     plot(Temps(:,ii), -depthsAtMiddleOfLayers, 'g');
% end
% 
% % Plot min and max temps for that depth
% p2(2) = plot(minT, -depthsAtMiddleOfLayers, 'b');
% p2(3) = plot(maxT, -depthsAtMiddleOfLayers, 'r');
% 
% legend(p2, 'Annual Average', 'Minimum', 'Maximum');
% 
% ylim([-2 0]);


% %% Saves output of a flat run as scattered IR and Vis radiation to be used in a sloped case
if saveFlat == 1
   
   for n = 2:nStepsInYear
       if frostMasses(n) > 0
           flatIR(n) = emisFrost * sigmasb * Tsurf(n-1)^4;
           flatVis(n) = Afrost * sf(n);
       else
           flatIR(n) = emis * sigmasb * Tsurf(n-1)^4;
           flatVis(n) = A * sf(n);
       end
           
   end
   
   if frostMasses(1) > 0
        flatIR(1) = emisFrost * sigmasb * Tsurf(end)^4;
        flatVis(1) = Afrost * sf(1);
   else
        flatIR(1) = emis * sigmasb * Tsurf(end)^4;
        flatVis(1) = A * sf(1);
   end
    
   save('flatSaved', 'flatIR', 'flatVis');
end

% %% Save surface temperatures, if that's what you want to do, yo!
if saveTsurf == 1
    save(TsurfSaved, 'lsWrapped', 'Tsurf');
end
    
