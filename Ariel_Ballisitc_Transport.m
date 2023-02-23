%BALLISTIC TRANSPORT MODEL 
%Author: Stephanie Menten

%% INPUT PARAMETERS

g = 0.249;                              % gravity Ariel in m/s^2
r = 578900;                             % average radius Ariel in m
% m = 1.27e21;                            % mass of Charon in kg
miu = 7.3106e-26;                       % molar mass one CO2 molecule in kg
% p = 8.93e40;                            % number of molecules of CH4
sa = 4*pi*(r^2);                        % surface area of Ariel in m^2   
rho_m = 522 ;                           % density of methane ice in kg/m^3
rho_c = 1562 ;                          % density of carbon dioxide ice in kg/m^3
G = 6.67e-11;                           % Gravitational constant in m^3kg^-1s^-2
sigmasb = 5.67e-8;                      % Stefan Boltzmann constant W m^-2 K^-4
au  = 1.4959787061e11;                  % Size of one astronomical unit in meters
kb = 1.38e-23;                          % Boltzmann constant in J/K

%% One Hop
lat_i = -90 + (90-(-90)) .* rand(1,1) ;     % random Latitude between -90 and 90
lat = lat_i ;
[T] = chooseTemperature(lat, 270)  ;        %avg Charon Surface Temp
lon_i = 0 + (360-0) .* rand(1,1) ;          %random Longitude between 0 and 360 
theta = 0.5*acosd((2*(0 + (1-0) .* rand(1,1)))-1) ;  %weighted random theta value between 0 and 90 
azimuth = 0 + (360-0) .* rand(1,1) ;        %random azimuth value between 0 and 360
T0 = T ; 
int0 = 0 + (1-0) .* rand(1,1) ;             %random integral value picked from integrated s distribution
[ v ] = find_velocity(T0, int0) ;           %velocity picked from s distribution between 0 and 591 m/s (escape velocity)


% CONSTANTS

G = 6.67e-11;                               % Gravitational constant in m^3kg^-1s^-2
sigmasb = 5.67e-8;                          % Stefan Boltzmann constant W m^-2 K^-4
au  = 1.4959787061e11;                      % Size of one astronomical unit in meters
kb = 1.38e-23;                              %Boltzmann constant in J/K

% CONVERSIONS

ev = sqrt(2*g*r);               %escape velocity of Charon m/s
deg2rad = pi/180 ;
%atoau = a*au ;                 %orbital major axis in m
rad2deg = 180/pi ;
theta_r = theta*deg2rad ;       %converting to radians
azimuth_r = azimuth*deg2rad ;
lat_ir = lat_i*deg2rad ;
lon_ir = lon_i*deg2rad ;


% Mordor macula starts at 63 degrees lat, more pronounced red occurs around
% 70 or 75 degrees

a = (1/((2/r)-((v^2)/(m*G))));
%calculate CH4 orbit's semimajor axis
e = sqrt((1-((sind(theta)^2)*r*(2*a-r))/a^2));
%calculate eccentricity 

f = acos((a-(a*(e^2))-r)/(r*e));
%finding the true anomaly 
d =  ((2*pi)-(2*f))*r ;
% distance a molecule travels (through the air, not over a surface)




% Haversine Equation 
%Solving with the distance for a new lat/long coordinate

lat_fr = asin(sin(lat_ir)*cos(d/r) + cos(lat_ir)*sin(d/r)*cos(azimuth_r)) ;

lon_fr = lon_ir + atan2(sin(azimuth_r)*sin(d/r)*cos(lat_ir), cos(d/r)-sin(lat_ir)*sin(lat_fr)) ;

%atan2(a,b)    = arc tan(b/a)  

lat_f = lat_fr*rad2deg ; %changing radians to degrees

lon_f = lon_fr*rad2deg ;

%% Multiple Hops and Particles
% for loop = multiple particles
% while loop = multiple hops
latituden90 = 0;
latituden85 = 0;
latituden80 = 0;
latituden75 = 0;
latituden70 = 0;
latituden65 = 0;
latituden60 = 0;
latituden55 = 0;
latituden50 = 0;
latituden45 = 0;
latituden40 = 0;
latituden35 = 0;
latituden30 = 0;
latituden25 = 0;
latituden20 = 0;
latituden15 = 0;
latituden10 = 0;
latituden5 = 0;
latitude0 = 0;
latitude5 = 0;
latitude10 = 0;
latitude15 = 0;
latitude20 = 0;
latitude25 = 0;
latitude30 = 0;
latitude35 = 0;
latitude40 = 0;
latitude45 = 0;
latitude50 = 0;
latitude55 = 0;
latitude60 = 0;
latitude65 = 0;
latitude70 = 0;
latitude75 = 0;
latitude80 = 0;
latitude85 = 0;
latitude90 = 0;


longitude0 = 0;
longitude10 = 0;
longitude20 = 0;
longitude30 = 0;
longitude40 = 0;
longitude50 = 0;
longitude60 = 0;
longitude70 = 0;
longitude80 = 0;
longitude90 = 0;
longitude100 = 0;
longitude110 = 0;
longitude120 = 0;
longitude130 = 0;
longitude140 = 0;
longitude150 = 0;
longitude160 = 0;
longitude170 = 0;
longitude180 = 0;
longitude190 = 0;
longitude200 = 0;
longitude210 = 0;
longitude220 = 0;
longitude230 = 0;
longitude240 = 0;
longitude250 = 0;
longitude260 = 0;
longitude270 = 0;
longitude280 = 0;
longitude290 = 0;
longitude300 = 0;
longitude310 = 0;
longitude320 = 0;
longitude330 = 0;
longitude340 = 0;
longitude350 = 0;
longitude360 = 0;


escape = 0;


for z = 1:10000

    

lat=[];
lon=[];
lat(1) = -30 + (10-(-30)) .* rand(1,1) ;  
lon(1) = 120 + (210-120) .* rand(1,1);   
ev = sqrt(2*g*r); 

i = 2;
v=0 ;
while v < ev
     
[T] = chooseTemperature(lat(i-1), 330)  ;               %Charon Surface Temp
theta = 0.5*acosd((2*(0 + (1-0) .* rand(1,1)))-1) ;     %weighted random theta value between 0 and 90 
azimuth = 0 + (360-0) .* rand(1,1) ;                    %random azimuth value between 0 and 360
T0 = T ;
int0 = 0 + (1-0) .* rand(1,1) ;                         %random integral value picked from integrated s distribution
[v] = chooseVelocity(T0, int0) ;                         %velocity picked from s distribution between 0 and 591 m/s (escape velocity)

a = (1/((2/r)-((v^2)/(m*G))));                          %calculate CH4 orbit's semimajor axis
e = sqrt((1-((sind(theta)^2)*r*(2*a-r))/a^2));          %calculate eccentricity

f = acos((a-(a*(e^2))-r)/(r*e));                        %finding the true anomaly 
d =  ((2*pi)-(2*f))*r ;                                 %distance a molecule travels (through the air, not over a surface)

deg2rad = pi/180 ;
%atoau = a*au ;                  % orbital major axis in m
rad2deg = 180/pi ;
theta_r = theta*deg2rad ;        % converting to radians
azimuth_r = azimuth*deg2rad ;
lat_ir = lat(i-1)*deg2rad ;
lon_ir = lon(i-1)*deg2rad ;


% Haversine Equation 
%Solving with the distance for a new lat/long coordinate

lat_fr = asin(sin(lat_ir).*cos(d./r) + cos(lat_ir).*sin(d./r).*cos(azimuth_r)) ;

lon_fr = lon_ir + atan2(sin(azimuth_r).*sin(d./r).*cos(lat_ir), cos(d./r)-sin(lat_ir).*sin(lat_fr)) ;

%atan2(a,b)    = arc tan(b/a) 
lat(i) = lat_fr.*rad2deg ; %changing radians to degrees
lon(i) = lon_fr.*rad2deg ;

% Calculate Sublimation Rate
sinyr = 3.154e7 ;                           % Number of seconds in 1 Earth years
rho_c = 1500   ;   
p_ic = 101.3 ;

ppc = p_ic.*(exp(23.102-(3148./T))) ; %find sublimation rate for CO2
i_src = ppc.*sqrt(miu./(2.*pi.*kb.*T)) ;
i_sc = i_src./rho_c.*sinyr ;


    if v >= ev 
      % disp('Particle escapes the system') 
         break
    end
% 
    if i_sc > 1e-16
      %  disp('Particle enters Mordor') 
        break
    end
    




    i = i+1 ;
end

    if z == 100
        disp('100 runs')
    end

     if z == 500
        disp('500 runs')
     end
    
    if z == 1000
        disp('1000 runs')
    end
    
    if z == 2500
        disp('2500 runs')
    end
    
    if z == 5000
        disp('5000 runs')
    end
    
    if z == 7500
        disp('7500 runs')
    end
    
     if z == 9000
        disp('9000 runs')
     end
    
    if z == 10000
        disp('10000 runs')
    end

    % Latitude Binning
    if (lat(i) < -87.5) && (lat(i) >= -90)
       latituden90 = latituden90 +1 ;
    end
    if (lat(i) < -82.5) && (lat(i) >= -87.5)
       latituden85 = latituden85 +1 ;
    end
    if (lat(i) < -77.5) && (lat(i) >= -82.5)
       latituden80 = latituden80 +1 ;
    end
    if (lat(i) < -72.5) && (lat(i) >= -77.5)
       latituden75 = latituden75 +1 ;
    end
    if (lat(i) < -67.5) && (lat(i) >= -72.5)
       latituden70 = latituden70 +1 ;
    end
    if (lat(i) < -62.5) && (lat(i) >= -67.5)
       latituden65 = latituden65 +1 ;
    end
    if (lat(i) < -57.5) && (lat(i) >= -62.5)
       latituden60 = latituden60 +1 ;
    end
    if (lat(i) < -52.5) && (lat(i) >= -57.5)
       latituden55 = latituden55 +1 ;
    end
    if (lat(i) < -47.5) && (lat(i) >= -52.5)
       latituden50 = latituden50 +1 ;
    end
    if (lat(i) < -42.5) && (lat(i) >= -47.5)
       latituden45 = latituden45 +1 ;
    end
    if (lat(i) < -37.5) && (lat(i) >= -42.5)
       latituden40 = latituden40 +1 ;
    end
    if (lat(i) < -32.5) && (lat(i) >= -37.5)
       latituden35 = latituden35 +1 ;
    end
    if (lat(i) < -27.5) && (lat(i) >= -32.5)
       latituden30 = latituden30 +1 ;
    end
    if (lat(i) < -22.5) && (lat(i) >= -27.5)
       latituden25 = latituden25 +1 ;
    end
    if (lat(i) < -17.5) && (lat(i) >= -22.5)
       latituden20 = latituden20 +1 ;
    end
    if (lat(i) < -12.5) && (lat(i) >= -17.5)
       latituden15 = latituden15 +1 ;
    end
    if (lat(i) < -7.5) && (lat(i) >= -12.5)
       latituden10 = latituden10 +1 ;
    end
    if (lat(i) < -2.5) && (lat(i) >= -7.5)
       latituden5 = latituden5 +1 ;
    end  
    if (lat(i) < 2.5) && (lat(i) >= -2.5)
       latitude0 = latitude0 +1 ;
    end
    if (lat(i) < 7.5) && (lat(i) >= 2.5)
       latitude5 = latitude5 +1 ;
    end
    if (lat(i) < 12.5) && (lat(i) >= 7.5)
       latitude10 = latitude10 +1 ;
    end
    if (lat(i) < 17.5) && (lat(i) >= 12.5)
       latitude15 = latitude15 +1 ;
    end
    if (lat(i) < 22.5) && (lat(i) >= 17.5)
       latitude20 = latitude20 +1 ;
    end
    if (lat(i) < 27.5) && (lat(i) >= 22.5)
       latitude25 = latitude25 +1 ;
    end
    if (lat(i) < 32.5) && (lat(i) >= 27.5)
       latitude30 = latitude30 +1 ;
    end
    if (lat(i) < 37.5) && (lat(i) >= 32.5)
       latitude35 = latitude35 +1 ;
    end
    if (lat(i) < 42.5) && (lat(i) >= 37.5)
       latitude40 = latitude40 +1 ;
    end
    if (lat(i) < 47.5) && (lat(i) >= 42.5)
       latitude45 = latitude45 +1 ;
    end
    if (lat(i) < 52.5) && (lat(i) >= 47.5)
       latitude50 = latitude50 +1 ;
    end
    if (lat(i) < 57.5) && (lat(i) >= 52.5)
       latitude55 = latitude55 +1 ;
    end
    if (lat(i) < 62.5) && (lat(i) >= 57.5)
       latitude60 = latitude60 +1 ;
    end
    if (lat(i) < 67.5) && (lat(i) >= 62.5)
       latitude65 = latitude65 +1 ;
    end
    if (lat(i) < 72.5) && (lat(i) >= 67.5)
       latitude70 = latitude70 +1 ;
    end
    if (lat(i) < 77.5) && (lat(i) >= 72.5)
       latitude75 = latitude75 +1 ;
    end
    if (lat(i) < 82.5) && (lat(i) >= 77.5)
       latitude80 = latitude80 +1 ;
    end
    if (lat(i) < 87.5) && (lat(i) >= 82.5)
       latitude85 = latitude85 +1 ;
    end
    if (lat(i) <= 90) && (lat(i) >= 87.5)
       latitude90 = latitude90 +1 ;
    end

% Longitude Binning
    if (lon(i) < 5) && (lon(i) >= 0)
       longitude0 = longitude0 +1 ;
    end
    if (lon(i) < 15) && (lon(i) >= 5)
       longitude10 = longitude10 +1 ;
    end
    if (lon(i) < 25) && (lon(i) >= 15)
       longitude20 = longitude20 +1 ;
    end
    if (lon(i) < 35) && (lon(i) >= 25)
       longitude30 = longitude30 +1 ;
    end
    if (lon(i) < 45) && (lon(i) >= 35)
       longitude40 = longitude40 +1 ;
    end
    if (lon(i) < 55) && (lon(i) >= 45)
       longitude50 = longitude50 +1 ;
    end
    if (lon(i) < 65) && (lon(i) >= 55)
       longitude60 = longitude60 +1 ;
    end
    if (lon(i) < 75) && (lon(i) >= 65)
       longitude70 = longitude70 +1 ;
    end
    if (lon(i) < 85) && (lon(i) >= 75)
       longitude80 = longitude80 +1 ;
    end
    if (lon(i) < 95) && (lon(i) >= 85)
       longitude90 = longitude90 +1 ;
    end
    if (lon(i) < 105) && (lon(i) >= 95)
       longitude100 = longitude100 +1 ;
    end
    if (lon(i) < 115) && (lon(i) >= 105)
       longitude110 = longitude110 +1 ;
    end
    if (lon(i) < 125) && (lon(i) >= 115)
       longitude120 = longitude120 +1 ;
    end
    if (lon(i) < 135) && (lon(i) >= 125)
       longitude130 = longitude130 +1 ;
    end
    if (lon(i) < 145) && (lon(i) >= 135)
       longitude140 = longitude140 +1 ;
    end
    if (lon(i) < 155) && (lon(i) >= 145)
       longitude150 = longitude150 +1 ;
    end
    if (lon(i) < 165) && (lon(i) >= 155)
       longitude160 = longitude160 +1 ;
    end
    if (lon(i) < 175) && (lon(i) >= 165)
       longitude170 = longitude170 +1 ;
    end
    if (lon(i) < 185) && (lon(i) >= 175)
       longitude180 = longitude180 +1 ;
    end
    if (lon(i) < 195) && (lon(i) >= 185)
       longitude190 = longitude190 +1 ;
    end
    if (lon(i) < 205) && (lon(i) >= 195)
       longitude200 = longitude200 +1 ;
    end
    if (lon(i) < 215) && (lon(i) >= 205)
       longitude210 = longitude210 +1 ;
    end
    if (lon(i) < 225) && (lon(i) >= 215)
       longitude220 = longitude220 +1 ;
    end
    if (lon(i) < 235) && (lon(i) >= 225)
       longitude230 = longitude230 +1 ;
    end
    if (lon(i) < 245) && (lon(i) >= 235)
       longitude240 = longitude240 +1 ;
    end
    if (lon(i) < 255) && (lon(i) >= 245)
       longitude250 = longitude250 +1 ;
    end
    if (lon(i) < 265) && (lon(i) >= 255)
       longitude260 = longitude260 +1 ;
    end
    if (lon(i) < 275) && (lon(i) >= 265)
       longitude270 = longitude270 +1 ;
    end
    if (lon(i) < 285) && (lon(i) >= 275)
       longitude280 = longitude280 +1 ;
    end
    if (lon(i) < 295) && (lon(i) >= 285)
       longitude290 = longitude290 +1 ;
    end
    if (lon(i) < 305) && (lon(i) >= 295)
       longitude300 = longitude300 +1 ;
    end
    if (lon(i) < 315) && (lon(i) >= 305)
       longitude310 = longitude310 +1 ;
    end
    if (lon(i) < 325) && (lon(i) >= 315)
       longitude320 = longitude320 +1 ;
    end
    if (lon(i) < 335) && (lon(i) >= 325)
       longitude330 = longitude330 +1 ;
    end
    if (lon(i) < 345) && (lon(i) >= 335)
       longitude340 = longitude340 +1 ;
    end
    if (lon(i) < 355) && (lon(i) >= 345)
       longitude350 = longitude350 +1 ;
    end
    if (lon(i) <= 360) && (lon(i) >= 355)
       longitude360 = longitude360 +1 ;
    end


    if v >= ev
       escape = escape+1 ;
    end

end

saveMordor = 0;                     % If this is 1, the script saves mordor and escape stats to the file given below; set to 0 if you don't want to save output
 MordorSaved = 'MordorEscapeLs330.mat';    % File to store surface mordor/escape stats in


if saveMordor == 1
     save(MordorSaved, 'mordor', 'escape');
end
% plot ([lon_i lon_f],[lat_i lat_f])
% hold on
% plot (lon_i, lat_i, 'b*')
% plot (lon_f, lat_f, 'r*')
% 
% 
% 
% hold off
% 
% xlim ([0 360]);
% ylim ([-90 90]);

%% Bar Graph Latitude Plots

lats =[] ;



%% Bar Graph Longitude Plots 
lons =[] ;



%% THICKNESS CALCULATION

SA_MMkm = 242190.8 ;                         % in kilometers^2
t_Charon = 551860 ;                     % 1 Charon year in s

SA_MM = (SA_MMkm*(1000^2)) ;                       % in m^2
Methane_mol = (8.93e40)*0.6699 ;
Methane_mass = (2.38e15)*.6699 ;       % in kg
MM_vol = Methane_mass/rho_m ;           % in m^3

MM_thickness = MM_vol/SA_MM ;           % in meters

MM_duration = Methane_mol/(SA_MM*t_Charon) ;   %in mol/m^2s
MM_duration1Myr = Methane_mol/(SA_MM*(3.1536e13)) ;
MM_duration1Gyr = Methane_mol/(SA_MM*(3.1536e16)) ;

saveThickness = 1;                     % If this is 1, the script saves MM thickness, volume, and emplacement duration stats to the file given below; set to 0 if you don't want to save output
ThicknessSaved = 'MMThickness.mat';    % File to store thicknesses in


if saveThickness == 1
     save(ThicknessSaved, 'MM_thickness', 'MM_vol', 'MM_duration' , 'MM_duration1Myr', 'MM_duration1Gyr' );
end

%% SUBLIMATION PLOTS

lat = -90:90 ;


for i = 1:length(lat)
T(i) = chooseTemperature(lat(i), 270)  ;

end

%%
tp = 90.69       ;                           % triple point T in Kelvin
L_c = 8190 ;                                 %Latent heat of vaporization J/mol
R_c = 8.313 ;                               % Gas constant J/molK
p_i = 11696   ;                             %  P triple point in Pa
sinyr = 3.154e7 ;                           % Number of seconds in 1 Earth years


pp = p_i.*(exp((L_c./R_c).*((1./tp)-(1./T)))) ; 

i_sr = pp.*sqrt(miu./(2.*pi.*kb.*T)) ; 
i_s = i_sr./rho_m.*sinyr ;

saveSublimation = 0;                     % If this is 1, the script saves MM thickness, volume, and emplacement duration stats to the file given below; set to 0 if you don't want to save output
SublimationSaved = 'Ls240SublimationRate.mat';    % File to store thicknesses in


if saveSublimation == 1
     save(SublimationSaved, 'lat', 'i_s0');
end

%% Sublimation vs Temperature Plot

T = 0:200 ;                                 %Temperature in Kelvin
tp = 90.69       ;                           % triple point T in Kelvin
L_c = 8190 ;                                 %Latent heat of vaporization J/mol
R_c = 8.313 ;                               % Gas constant J/molK
p_i = 11696   ;                             %  P triple point in Pa
sinyr = 3.154e7 ;                           % Number of seconds in 1 Earth years
rho_c = 1500   ;                            % Density of CO2 ice in kg/m^3

pp = p_i.*(exp((L_c./R_c).*((1./tp)-(1./T)))) ; 

i_sr = pp.*sqrt(miu./(2.*pi.*kb.*T)) ;
i_s = i_sr./rho_m.*sinyr ;                  %convert to m/yr

p_ic = 101.3 ;


ppc = p_ic.*(exp(23.102-(3148./T))) ; %find sublimation rate for CO2
i_src = pp.*sqrt(miu./(2.*pi.*kb.*T)) ;
i_sc = i_src./rho_c.*sinyr ;


plot (T,i_s, 'LineWidth', 1)
hold on

plot (T, i_sc, 'LineWidth', 1)
hold off

%xline(35,'--')
xlabel('Temperature (K)')
ylabel('Sublimation Rate (m/yr)')
xlim ([0 200]);
ylim ([10^-10 10^10]) ;
ax = gca;
ax.YAxis.Exponent = 0;
set(gca, 'YScale', 'log');




%% Sublimation Rates by Time of Year

x = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180 185 190 195 200 205 210 215 220 225 230 235 240 245 250 255 260 265 270 275 280 285 290 295 300 305 310 315 320 325 330 335 340 345 350 355 360] ;
y = [3.9513e-14 6.7562e-5 0.0702 1.8059 13.8274 59.1228 179.7443 436.296 901.7511 1.6499e3 2740.8 4205.2 6.0317e3 8158.5 1.0473e4 12820 15016 16868 1.8202e4 18880 18825 18032 1.6567e4 14557 12204 9707.3 7278.1 5.096e3 3288.7 1920.1 986.2695 426.7098 144.4943 33.4605 3.9664 0.1097 2.2563e-5 0.00000051811 0.000000085704 0.000000022262 7.3407e-9 0.0000000027802 0.0000000011756 0.00000000053845 0.0000000002628 1.3518e-10 0.00000000007273 0.000000000040693 0.00000000002358 0.000000000014107 8.6932e-12 0.000000000005507 0.000000000003582 0.0000000000023888 1.6319e-12 0.000000000001141 0.00000000000081585 0.0000000000005961 0.0000000000004447 0.00000000000033845 2.6258e-13 0.00000000000020747 0.0000000000001668 1.3631e-13 0.00000000000011313 0.000000000000095258 8.1294e-14 0.000000000000070247 0.000000000000061403 0.000000000000054244 0.000000000000048387 0.000000000000043547 3.9513e-14] ; %North Pole
y1 = [9.2786e-7 2.7495e-8 6.3268e-9 2.3244e-9 1.0832e-9 5.8522e-10 3.5005e-10 0.00000000022562 0.00000000015397 1.0991e-10 0.000000000081349 0.000000000062026 4.8474e-11 0.000000000038678 3.1408e-11 0.00000000002589 0.000000000021617 0.000000000018251 1.5557e-11 0.000000000013371 0.000000000011575 0.000000000010083 8.8312e-12 0.0000000000065158 0.0000000000057446 0.0000000000050823 0.0000000000045096 4.0115e-12 0.0000000000035756 0.0000000000031923 2.8537e-12 0.0000000000025535 0.0000000000022863 0.0000000000020477 0.0000000000018343 0.0000000000016428 1.4708e-12 0.0012 0.3883 4.8546 21.1908 56.0074 110.6269 180.7658 258.521 334.9388 402.0969 454.4013 488.9405 505.2722 504.8433 490.2805 464.7415 431.4099 393.1529 352.3517 310.8369 269.9343 230.5483 193.2655 158.4809 126.4833 97.5457 71.9651 50.0785 32.2125 18.5905 9.1904 3.6023 0.9715 0.1333 0.0041 9.2786e-7] ; %South Pole

plot (x,y, 'LineWidth', 3)

hold on

plot (x, y1, 'LineWidth', 3)

hold off

xlabel('L_s')
ylabel('Sublimation Rate (m/yr)')
xlim ([0 360]);
%ylim ([-90 90]);
ax = gca;
ax.YAxis.Exponent = 0;
set(gca, 'YScale', 'log')

%% LSPC Crater Density Plots

x = [5 10 20 30 40 50 100] ;
y = [6.606e-5 4.129e-5 2.064e-5 1.239e-5 1.239e-5 8.258e-6 4.129e-6] ;

plot (x,y, 'LineWidth', 3)
xlabel('Diameter')
ylabel('Cumulative Crater Number (km^-^2)')
xlim ([1 100]);
%ylim ([-90 90]);
ax = gca;
ax.YAxis.Exponent = 0;
ax.XAxis.Exponent = 0;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')



%% LPSC Methane Solubility

x = [.497 .717 1.0 1.317 1.678 2.235 2.585 3.11 3.66 4.17 4.37 4.45] ;
y = [341 4.129e-5 2.064e-5 1.239e-5 1.239e-5 8.258e-6 4.129e-6] ;


plot (x,y, 'LineWidth', 1)
xlabel('Diameter')
ylabel('Cumulative Crater Number (km^-^2)')
xlim ([0 100]);
%ylim ([-90 90]);




%% VELOCITY PLOTS

T11 = 15;         %Charon Surface Temp in K

syms v

s = (4*pi)*((miu/(2*pi*kb*T11))^(3/2))*(v^2)*(exp((-1*miu*(v^2))/(2*kb*T11))); 
%finding the fraction of material that has a velocity fall above the escape
%velocity



fplot (s, [0 800],'LineWidth', 3)

hold on
T1 = 40;         %Charon Surface Temp in K

syms v1

s1 = (4*pi)*((miu/(2*pi*kb*T1))^(3/2))*(v1^2)*(exp((-1*miu*(v1^2))/(2*kb*T1))); 
%finding the fraction of material that has a velocity fall above the escape
%velocity



fplot (s1, [0 800],'LineWidth', 3)

T2 = 60;         %Charon Surface Temp in K

syms v2

s2 = (4*pi)*((miu/(2*pi*kb*T2))^(3/2))*(v2^2)*(exp((-1*miu*(v2^2))/(2*kb*T2))); 
%finding the fraction of material that has a velocity fall above the escape
%velocity


fplot (s2, [0 800],'LineWidth', 3)

 xline(591,'LineWidth', 1)

hold off

ax = gca;
ax.YAxis.Exponent = 0;

xlabel('Velocity (m/s)')
ylabel('Probability Distribution')

%T.^2 

% dot notation extremely helpful! use it 


%integrate area above escape velocity
 
%% JEANS ESCAPE FRACTION
%Equations sourced from the David Catling Textbook bc I couldn't find them
%credibly anywhere else 

%looking for the escape fraction of molecules of CH4

T=1:100;

%ef = int((4.*pi).*((miu./(2.*pi.*kb.*T)).^(3./2)).*(v.^2).*(exp((-1.*miu.*(v.^2))./(2.*kb.*T))), 591, inf) ;
syms v
ef = int((4*pi)*((miu./(2*pi*kb.*T)).^(3/2)).*(v.^2).*(exp((-1*miu.*(v.^2))./(2*kb.*T))), v, 591, inf) ;
figure

%ef1 = log10(ef) ;


plot (T, ef,'LineWidth', 3)

xlabel('Temperature (K)')
ylabel('Escape Fraction')