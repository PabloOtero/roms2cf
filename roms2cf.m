function roms2cf(grdname,romsfile,cffile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create a netcdf his or avg file CF-complaint
%       grdname: name of the grid file
%       romsfile: name of the roms-agrif his or avg file
%       cffile: name of the new file to create
%
%       Pablo Otero, June-2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;
%clear all;
%grdname='/data/Roms_simula/SimulaRaia/op_conf/Raia_grd3.nc';
%romsfile='/data/Roms_simula/SimulaRaia/op_out/avg/Raia_20120302_avg.nc';
%cffile='/data/Roms_simula/SimulaRaia/op_out/testing_CF_avg.nc';

% Did you simulate a fixed time of year (climate experiment)?
climate_experiment=0;


%%%%--------------END OF USER CONFIGURATION------------%%%


ncgrd=netcdf(grdname);
ncold=netcdf(romsfile);
nw = netcdf(cffile, 'clobber');
result = redef(nw);

theNewFillValue=1.e+37;

%
%  Create dimensions
%
disp(['Creating the CF-compliant file structure']);
nw('xi_rho') = length(ncgrd('xi_rho'));
nw('xi_u') = length(ncgrd('xi_u'));
nw('xi_v') = length(ncgrd('xi_v'));
%nw('xi_psi') = length(ncgrd('xi_psi'));
nw('eta_rho') = length(ncgrd('eta_rho'));
nw('eta_u') = length(ncgrd('eta_u'));
nw('eta_v') = length(ncgrd('eta_v'));
%nw('eta_psi') = length(ncgrd('eta_psi'));
nw('N') = length(ncold('s_rho'));
nw('s_rho') = length(ncold('s_rho'));
nw('s_w') = length(ncold('s_w'));
nw('scalar') = 1;
nw('scrum_time') = length(ncold{'scrum_time'}(:));



%
%  Create variables and attributes
%
nw{'xl'} = ncdouble('scalar');
nw{'xl'}.long_name = ncchar('domain length in the XI-direction');
nw{'xl'}.long_name = 'domain length in the XI-direction';
nw{'xl'}.units = ncchar('meter');
nw{'xl'}.units = 'meter';

nw{'el'} = ncdouble('scalar');
nw{'el'}.long_name = ncchar('domain length in the ETA-direction');
nw{'el'}.long_name = 'domain length in the ETA-direction';
nw{'el'}.units = ncchar('meter');
nw{'el'}.units = 'meter';

nw{'theta_s'} = ncdouble('scalar');
nw{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nw{'theta_s'}.long_name = 'S-coordinate surface control parameter';

nw{'theta_b'} = ncdouble('scalar');
nw{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nw{'theta_b'}.long_name = 'S-coordinate bottom control parameter';

nw{'Tcline'} = ncdouble('scalar');
nw{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nw{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nw{'Tcline'}.units = ncchar('meter');
nw{'Tcline'}.units = 'meter';

nw{'hc'} = ncdouble('scalar');
nw{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nw{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nw{'hc'}.units = ncchar('meter');
nw{'hc'}.units = 'meter';

nw{'s_rho'} = ncdouble('s_rho');
nw{'s_rho'}.long_name = ncchar('S-coordinate at RHO-points');
nw{'s_rho'}.long_name = 'S-coordinate at RHO-points';
nw{'s_rho'}.valid_min = -1.;
nw{'s_rho'}.valid_max = 0.;
nw{'s_rho'}.standard_name = ncchar('ocean_s_coordinate');
nw{'s_rho'}.standard_name = 'ocean_s_coordinate';
nw{'s_rho'}.formula_terms = ncchar('s: s_rho eta: zeta depth: h a: theta_s b: theta_b depth_c: hc');
nw{'s_rho'}.formula_terms = 's: s_rho eta: zeta depth: h a: theta_s b: theta_b depth_c: hc';
nw{'s_rho'}.field = ncchar('s_rho, scalar');
nw{'s_rho'}.field = 's_rho, scalar';

nw{'s_w'} = ncdouble('s_w');
nw{'s_w'}.long_name = ncchar('S-coordinate at W-points');
nw{'s_w'}.long_name = 'S-coordinate at W-points';
nw{'s_w'}.valid_min = -1.;
nw{'s_w'}.valid_max = 0.;
nw{'s_w'}.standard_name = ncchar('ocean_s_coordinate');
nw{'s_w'}.standard_name = 'ocean_s_coordinate';
nw{'s_w'}.formula_terms = ncchar('s: s_w eta: zeta depth: h a: theta_s b: theta_b depth_c: hc');
nw{'s_w'}.formula_terms = 's: s_w eta: zeta depth: h a: theta_s b: theta_b depth_c: hc';
nw{'s_w'}.field = ncchar('s_w, scalar');
nw{'s_w'}.field = 's_w, scalar';

nw{'Cs_r'} = ncdouble('s_rho');
nw{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nw{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nw{'Cs_r'}.valid_min = -1.;
nw{'Cs_r'}.valid_max = 0.;
nw{'Cs_r'}.field = ncchar('Cs_r, scalar');
nw{'Cs_r'}.field = 'Cs_r, scalar';

nw{'Cs_w'} = ncdouble('s_w');
nw{'Cs_w'}.long_name = ncchar('S-coordinate stretching curves at W-points');
nw{'Cs_w'}.long_name = 'S-coordinate stretching curves at W-points';
nw{'Cs_w'}.valid_min = -1.;
nw{'Cs_w'}.valid_max = 0.;
nw{'Cs_w'}.field = ncchar('Cs_w, scalar');
nw{'Cs_w'}.field = 'Cs_w, scalar';

nw{'h'} = ncdouble('eta_rho', 'xi_rho');
nw{'h'}.long_name = ncchar('bathymetry at RHO-points');
nw{'h'}.long_name = 'S-coordinate stretching curves at W-points';
nw{'h'}.units = ncchar('meter');
nw{'h'}.units = 'meter';
nw{'h'}.coordinates = ncchar('lon_rho lat_rho');
nw{'h'}.coordinates = 'lon_rho lat_rho';
nw{'h'}.field = ncchar('bath, scalar');
nw{'h'}.field = 'bath, scalar';

nw{'lon_rho'} = ncdouble('eta_rho', 'xi_rho');
nw{'lon_rho'}.long_name = ncchar('longitude of RHO-points');
nw{'lon_rho'}.long_name = 'longitude of RHO-points';
nw{'lon_rho'}.units = ncchar('degree_east');
nw{'lon_rho'}.units = 'degree_east';
nw{'lon_rho'}.field = ncchar('lon_rho, scalar');
nw{'lon_rho'}.field = 'lon_rho, scalar';

nw{'lat_rho'} = ncdouble('eta_rho', 'xi_rho');
nw{'lat_rho'}.long_name = ncchar('latitude of RHO-points');
nw{'lat_rho'}.long_name = 'latitude of RHO-points';
nw{'lat_rho'}.units = ncchar('degree_north');
nw{'lat_rho'}.units = 'degree_north';
nw{'lat_rho'}.field = ncchar('lat_rho, scalar');
nw{'lat_rho'}.field = 'lat_rho, scalar';

nw{'lon_u'} = ncdouble('eta_u', 'xi_u');
nw{'lon_u'}.long_name = ncchar('longitude of U-points');
nw{'lon_u'}.long_name = 'longitude of U-points';
nw{'lon_u'}.units = ncchar('degree_east');
nw{'lon_u'}.units = 'degree_east';
nw{'lon_u'}.field = ncchar('lon_u, scalar');
nw{'lon_u'}.field = 'lon_u, scalar';

nw{'lat_u'} = ncdouble('eta_u', 'xi_u');
nw{'lat_u'}.long_name = ncchar('latitude of U-points');
nw{'lat_u'}.long_name = 'latitude of U-points';
nw{'lat_u'}.units = ncchar('degree_north');
nw{'lat_u'}.units = 'degree_north';
nw{'lat_u'}.field = ncchar('lat_u, scalar');
nw{'lat_u'}.field = 'lat_u, scalar';

nw{'lon_v'} = ncdouble('eta_v', 'xi_v');
nw{'lon_v'}.long_name = ncchar('longitude of V-points');
nw{'lon_v'}.long_name = 'longitude of V-points';
nw{'lon_v'}.units = ncchar('degree_east');
nw{'lon_v'}.units = 'degree_east';
nw{'lon_v'}.field = ncchar('lon_v, scalar');
nw{'lon_v'}.field = 'lon_v, scalar';

nw{'lat_v'} = ncdouble('eta_v', 'xi_v');
nw{'lat_v'}.long_name = ncchar('latitude of V-points');
nw{'lat_v'}.long_name = 'latitude of V-points';
nw{'lat_v'}.units = ncchar('degree_north');
nw{'lat_v'}.units = 'degree_north';
nw{'lat_v'}.field = ncchar('lat_v, scalar');
nw{'lat_v'}.field = 'lat_v, scalar';

%nw{'lon_psi'} = ncdouble('eta_psi', 'xi_psi');
%nw{'lon_psi'}.long_name = ncchar('longitude of PSI-points');
%nw{'lon_psi'}.long_name = 'longitude of PSI-points';
%nw{'lon_psi'}.units = ncchar('degree_east');
%nw{'lon_psi'}.units = 'degree_east';
%nw{'lon_psi'}.field = ncchar('lon_psi, scalar');
%nw{'lon_psi'}.field = 'lon_psi, scalar';

%nw{'lat_psi'} = ncdouble('eta_psi', 'xi_psi');
%nw{'lat_psi'}.long_name = ncchar('latitude of PSI-points');
%nw{'lat_psi'}.long_name = 'latitude of PSI-points';
%nw{'lat_psi'}.units = ncchar('degree_north');
%nw{'lat_psi'}.units = 'degree_north';
%nw{'lat_psi'}.field = ncchar('lat_psi, scalar');
%nw{'lat_psi'}.field = 'lat_psi, scalar';

nw{'mask_rho'} = ncdouble('eta_rho', 'xi_rho');
nw{'mask_rho'}.long_name = ncchar('mask on RHO-points');
nw{'mask_rho'}.long_name = 'mask on RHO-points';
nw{'mask_rho'}.option_0 = ncchar('land');
nw{'mask_rho'}.option_0 = 'land';
nw{'mask_rho'}.option_1 = ncchar('water');
nw{'mask_rho'}.option_1 = 'water';
nw{'mask_rho'}.coordinates = ncchar('lon_rho lat_rho');
nw{'mask_rho'}.coordinates = 'lon_rho lat_rho';

nw{'mask_u'} = ncdouble('eta_u', 'xi_u');
nw{'mask_u'}.long_name = ncchar('mask on U-points');
nw{'mask_u'}.long_name = 'mask on U-points';
nw{'mask_u'}.option_0 = ncchar('land');
nw{'mask_u'}.option_0 = 'land';
nw{'mask_u'}.option_1 = ncchar('water');
nw{'mask_u'}.option_1 = 'water';
nw{'mask_u'}.coordinates = ncchar('lon_u lat_u');
nw{'mask_u'}.coordinates = 'lon_u lat_u';

nw{'mask_v'} = ncdouble('eta_v', 'xi_v');
nw{'mask_v'}.long_name = ncchar('mask on V-points');
nw{'mask_v'}.long_name = 'mask on V-points';
nw{'mask_v'}.option_0 = ncchar('land');
nw{'mask_v'}.option_0 = 'land';
nw{'mask_v'}.option_1 = ncchar('water');
nw{'mask_v'}.option_1 = 'water';
nw{'mask_v'}.coordinates = ncchar('lon_v lat_v');
nw{'mask_v'}.coordinates = 'lon_v lat_v';

%nw{'mask_psi'} = ncdouble('eta_psi', 'xi_psi');
%nw{'mask_psi'}.long_name = ncchar('mask on PSI-points');
%nw{'mask_psi'}.long_name = 'mask on PSI-points';
%nw{'mask_psi'}.option_0 = ncchar('land');
%nw{'mask_psi'}.option_0 = 'land';
%nw{'mask_psi'}.option_1 = ncchar('water');
%nw{'mask_psi'}.option_1 = 'water';
%nw{'mask_psi'}.coordinates = ncchar('lon_psi lat_psi');
%nw{'mask_psi'}.coordinates = 'lon_psi lat_psi';

nw{'scrum_time'} = ncdouble('scrum_time');
nw{'scrum_time'}.long_name = ncchar('time since initialization');
nw{'scrum_time'}.long_name = 'time since initialization';
if(climate_experiment)
  nw{'scrum_time'}.calendar = ncchar('none');
  nw{'scrum_time'}.calendar = 'none';
  nw{'scrum_time'}.units = ncchar('seconds');
  nw{'scrum_time'}.units = 'seconds';
else
  nw{'scrum_time'}.calendar = ncchar('gregorian');
  nw{'scrum_time'}.calendar = 'gregorian';
  nw{'scrum_time'}.units = ncchar('seconds since 2009-1-1 00:00:00');
  nw{'scrum_time'}.units = 'seconds since 2009-1-1 00:00:00';
end
nw{'scrum_time'}.field = ncchar('time, scalar, series');
nw{'scrum_time'}.field = 'time, scalar, series';


for count=1:length(var(ncold))
 namevar=name(var(ncold,count));
 if(strcmp(namevar,'zeta'))
  nw{'zeta'} = ncdouble('scrum_time', 'eta_rho', 'xi_rho');
  nw{'zeta'}.long_name = ncchar('free-surface elevation ');
  nw{'zeta'}.long_name = 'free-surface elevation';
  nw{'zeta'}.standard_name = ncchar('sea_surface_height_above_sea_level');
  nw{'zeta'}.standard_name = 'sea_surface_height_above_sea_level';
  nw{'zeta'}.units = ncchar('m');
  nw{'zeta'}.units = 'm';
  nw{'zeta'}.time = ncchar('scrum_time');
  nw{'zeta'}.time = 'scrum_time';
  nw{'zeta'}.coordinates = ncchar('lon_rho lat_rho scrum_time');
  nw{'zeta'}.coordinates = 'lon_rho lat_rho scrum_time';
  nw{'zeta'}.mask = ncchar('mask_rho');
  nw{'zeta'}.mask = 'mask_rho';
  nw{'zeta'}.field = ncchar('zeta, scalar, series');
  nw{'zeta'}.field = 'zeta, scalar, series';
 elseif(strcmp(namevar,'temp'))
  nw{'temp'} = ncdouble('scrum_time', 's_rho', 'eta_rho', 'xi_rho');
  nw{'temp'}.long_name = ncchar('potential temperature');
  nw{'temp'}.long_name = 'potential temperature';
  nw{'temp'}.standard_name = ncchar('sea_water_potential_temperature');
  nw{'temp'}.standard_name = 'sea_water_potential_temperature';
  nw{'temp'}.units = ncchar('Celsius');
  nw{'temp'}.units = 'Celsius';
  nw{'temp'}.time = ncchar('scrum_time');
  nw{'temp'}.time = 'scrum_time';
  nw{'temp'}.coordinates = ncchar('lon_rho lat_rho s_rho scrum_time');
  nw{'temp'}.coordinates = 'lon_rho lat_rho s_rho scrum_time';
  nw{'temp'}.mask = ncchar('mask_rho');
  nw{'temp'}.mask = 'mask_rho';
  nw{'temp'}.field = ncchar('temperature, scalar, series');
  nw{'temp'}.field = 'temperature, scalar, series';
 elseif(strcmp(namevar,'salt'))
  nw{'salt'} = ncdouble('scrum_time', 's_rho', 'eta_rho', 'xi_rho');
  nw{'salt'}.long_name = ncchar('salinity');
  nw{'salt'}.long_name = 'salinity';
  nw{'salt'}.standard_name = ncchar('sea_water_salinity');
  nw{'salt'}.standard_name = 'sea_water_salinity';
  nw{'salt'}.units = ncchar('1e-3');
  nw{'salt'}.units = '1e-3';
  nw{'salt'}.comment = ncchar('The unit of salinity is PSU, which is dimensionless.');
  nw{'salt'}.comment = 'The unit of salinity is PSU, which is dimensionless.';
  nw{'salt'}.time = ncchar('scrum_time');
  nw{'salt'}.time = 'scrum_time';
  nw{'salt'}.coordinates = ncchar('lon_rho lat_rho s_rho scrum_time');
  nw{'salt'}.coordinates = 'lon_rho lat_rho s_rho scrum_time';
  nw{'salt'}.mask = ncchar('mask_rho');
  nw{'salt'}.mask = 'mask_rho';
  nw{'salt'}.field = ncchar('salinity, scalar, series');
  nw{'salt'}.field = 'salinity, scalar, series';
 elseif(strcmp(namevar,'u'))
  nw{'u'} = ncdouble('scrum_time', 's_rho', 'eta_u', 'xi_u');
  nw{'u'}.long_name = ncchar('u-momentum component');
  nw{'u'}.long_name = 'u-momentum component';
  nw{'u'}.standard_name = ncchar('baroclinic_eastward_sea_water_velocity');
  nw{'u'}.standard_name = 'baroclinic_eastward_sea_water_velocity';
  nw{'u'}.units = ncchar('m s-1');
  nw{'u'}.units = 'm s-1';
  nw{'u'}.comment = ncchar('Positive when directed eastward (negative westward)');
  nw{'u'}.comment = 'Positive when directed eastward (negative westward)';
  nw{'u'}.time = ncchar('scrum_time');
  nw{'u'}.time = 'scrum_time';
  nw{'u'}.coordinates = ncchar('lon_u lat_u s_rho scrum_time');
  nw{'u'}.coordinates = 'lon_u lat_u s_rho scrum_time';
  nw{'u'}.mask = ncchar('mask_u');
  nw{'u'}.mask = 'mask_u';
  nw{'u'}.field = ncchar('u-velocity, scalar, series');
  nw{'u'}.field = 'u-velocity, scalar, series';
 elseif(strcmp(namevar,'v'))
  nw{'v'} = ncdouble('scrum_time', 's_rho', 'eta_v', 'xi_v');
  nw{'v'}.long_name = ncchar('v-momentum component');
  nw{'v'}.long_name = 'v-momentum component';
  nw{'v'}.standard_name = ncchar('baroclinic_northward_sea_water_velocity');
  nw{'v'}.standard_name = 'baroclinic_northward_sea_water_velocity';
  nw{'v'}.units = ncchar('m s-1');
  nw{'v'}.units = 'm s-1';
  nw{'v'}.comment = ncchar('Positive when directed northward (negative southward)');
  nw{'v'}.comment = 'Positive when directed northward (negative southward)';
  nw{'v'}.time = ncchar('scrum_time');
  nw{'v'}.time = 'scrum_time';
  nw{'v'}.coordinates = ncchar('lon_v lat_v s_rho scrum_time');
  nw{'v'}.coordinates = 'lon_v lat_v s_rho scrum_time';
  nw{'v'}.mask = ncchar('mask_v');
  nw{'v'}.mask = 'mask_v';
  nw{'v'}.field = ncchar('v-velocity, scalar, series');
  nw{'v'}.field = 'v-velocity, scalar, series';
 elseif(strcmp(namevar,'w'))
  nw{'w'} = ncdouble('scrum_time', 's_rho', 'eta_rho', 'xi_rho');
  nw{'w'}.long_name = ncchar('True vertical velocity');
  nw{'w'}.long_name = 'True vertical velocity';
  nw{'w'}.standard_name = ncchar('upward_sea_water_velocity');
  nw{'w'}.standard_name = 'upward_sea_water_velocity';
  nw{'w'}.units = ncchar('m s-1');
  nw{'w'}.units = 'm s-1';
  nw{'w'}.comment = ncchar('Positive when directed upward (negative downward)');
  nw{'w'}.comment = 'Positive when directed upward (negative downward)';
  nw{'w'}.time = ncchar('scrum_time');
  nw{'w'}.time = 'scrum_time';
  nw{'w'}.coordinates = ncchar('lon_rho lat_rho s_rho scrum_time');
  nw{'w'}.coordinates = 'lon_rho lat_rho s_rho scrum_time';
  nw{'w'}.mask = ncchar('mask_rho');
  nw{'w'}.mask = 'mask_rho';
  nw{'w'}.field = ncchar('w-velocity, scalar, series');
  nw{'w'}.field = 'w-velocity, scalar, series';
 elseif(strcmp(namevar,'omega'))
  nw{'omega'} = ncdouble('scrum_time', 's_w', 'eta_rho', 'xi_rho');
  nw{'omega'}.long_name = ncchar('Omega vertical velocity');
  nw{'omega'}.long_name = 'Omega vertical velocity';
  nw{'omega'}.standard_name = ncchar('upward_sea_water_velocity');
  nw{'omega'}.standard_name = 'upward_sea_water_velocity';
  nw{'omega'}.units = ncchar('m s-1');
  nw{'omega'}.units = 'm s-1';
  nw{'omega'}.comment = ncchar('Positive when directed upward (negative downward)');
  nw{'omega'}.comment = 'Positive when directed upward (negative downward)';
  nw{'omega'}.time = ncchar('scrum_time');
  nw{'omega'}.time = 'scrum_time';
  nw{'omega'}.coordinates = ncchar('lon_rho lat_rho s_w scrum_time');
  nw{'omega'}.coordinates = 'lon_rho lat_rho s_w scrum_time';
  nw{'omega'}.mask = ncchar('mask_rho');
  nw{'omega'}.mask = 'mask_rho';
  nw{'omega'}.field = ncchar('w-velocity, scalar, series');
  nw{'omega'}.field = 'w-velocity, scalar, series';
 elseif(strcmp(namevar,'ubar'))
  nw{'ubar'} = ncdouble('scrum_time', 'eta_u', 'xi_u');
  nw{'ubar'}.long_name = ncchar('u-barotropic velocity');
  nw{'ubar'}.long_name = 'u-barotropic velocity';
  nw{'ubar'}.standard_name = ncchar('barotropic_eastward_sea_water_velocity');
  nw{'ubar'}.standard_name = 'barotropic_eastward_sea_water_velocity';
  nw{'ubar'}.units = ncchar('m s-1');
  nw{'ubar'}.units = 'm s-1';
  nw{'ubar'}.comment = ncchar('Positive when directed eastward (negative westward)');
  nw{'ubar'}.comment = 'Positive when directed eastward (negative westward)';
  nw{'ubar'}.time = ncchar('scrum_time');
  nw{'ubar'}.time = 'scrum_time';
  nw{'ubar'}.coordinates = ncchar('lon_u lat_u scrum_time');
  nw{'ubar'}.coordinates = 'lon_u lat_u scrum_time';
  nw{'ubar'}.mask = ncchar('mask_u');
  nw{'ubar'}.mask = 'mask_u';
  nw{'ubar'}.field = ncchar('ubar, scalar, series');
  nw{'ubar'}.field = 'ubar, scalar, series';
 elseif(strcmp(namevar,'vbar'))
  nw{'vbar'} = ncdouble('scrum_time', 'eta_v', 'xi_v');
  nw{'vbar'}.long_name = ncchar('v-barotropic velocity');
  nw{'vbar'}.long_name = 'v-barotropic velocity';
  nw{'vbar'}.standard_name = ncchar('barotropic_northward_sea_water_velocity');
  nw{'vbar'}.standard_name = 'barotropic_northward_sea_water_velocity';
  nw{'vbar'}.units = ncchar('m s-1');
  nw{'vbar'}.units = 'm s-1';
  nw{'vbar'}.comment = ncchar('Positive when directed northward (negative southward)');
  nw{'vbar'}.comment = 'Positive when directed northward (negative southward)';
  nw{'vbar'}.time = ncchar('scrum_time');
  nw{'vbar'}.time = 'scrum_time';
  nw{'vbar'}.coordinates = ncchar('lon_v lat_v scrum_time');
  nw{'vbar'}.coordinates = 'lon_v lat_v scrum_time';
  nw{'vbar'}.mask = ncchar('mask_v');
  nw{'vbar'}.mask = 'mask_v';
  nw{'vbar'}.field = ncchar('vbar, scalar, series');
  nw{'vbar'}.field = 'vbar, scalar, series';
 elseif(strcmp(namevar,'rho'))
  nw{'rho'} = ncdouble('scrum_time', 's_rho', 'eta_rho', 'xi_rho');
  nw{'rho'}.long_name = ncchar('density anomaly');
  nw{'rho'}.long_name = 'density anomaly';
  nw{'rho'}.standard_name = ncchar('sea_water_sigma_theta');
  nw{'rho'}.standard_name = 'sea_water_sigma_theta';
  nw{'rho'}.units = ncchar('kg m-3');
  nw{'rho'}.units = 'kg m-3';
  nw{'rho'}.comment = ncchar('Usually referred to 1025 Kg m-3 in the ROMS code');
  nw{'rho'}.comment = 'Usually referred to 1025 Kg m-3 in the ROMS code';
  nw{'rho'}.time = ncchar('scrum_time');
  nw{'rho'}.time = 'scrum_time';
  nw{'rho'}.coordinates = ncchar('lon_rho lat_rho s_rho scrum_time');
  nw{'rho'}.coordinates = 'lon_rho lat_rho s_rho scrum_time';
  nw{'rho'}.mask = ncchar('mask_rho');
  nw{'rho'}.mask = 'mask_rho';
  nw{'rho'}.field = ncchar('rho, scalar, series');
  nw{'rho'}.field = 'rho, scalar, series';
 elseif(strcmp(namevar,'AKv'))
  nw{'AKv'} = ncdouble('scrum_time', 's_w', 'eta_rho', 'xi_rho');
  nw{'AKv'}.long_name = ncchar('vertical viscosity');
  nw{'AKv'}.long_name = 'vertical viscosity';
  nw{'AKv'}.standard_name = ncchar('ocean_vertical_momentum_diffusivity');
  nw{'AKv'}.standard_name = 'ocean_vertical_momentum_diffusivity';
  nw{'AKv'}.units = ncchar('m2 s-1');
  nw{'AKv'}.units = 'm2 s-1';
  nw{'AKv'}.comment = ncchar('Vertical component of the diffusivity of momentum due to motion which is not resolved on the grid scale of the model');
  nw{'AKv'}.comment = 'Vertical component of the diffusivity of momentum due to motion which is not resolved on the grid scale of the model';
  nw{'AKv'}.time = ncchar('scrum_time');
  nw{'AKv'}.time = 'scrum_time';
  nw{'AKv'}.coordinates = ncchar('lon_rho lat_rho s_w scrum_time');
  nw{'AKv'}.coordinates = 'lon_rho lat_rho s_w scrum_time';
  nw{'AKv'}.mask = ncchar('mask_rho');
  nw{'AKv'}.mask = 'mask_rho';
  nw{'AKv'}.field = ncchar('AKv, scalar, series');
  nw{'AKv'}.field = 'AKv, scalar, series';
 elseif(strcmp(namevar,'AKs'))
  nw{'AKs'} = ncdouble('scrum_time', 's_w', 'eta_rho', 'xi_rho');
  nw{'AKs'}.long_name = ncchar('Vertical diffusivity for salinity');
  nw{'AKs'}.long_name = 'Vertical diffusivity for salinity';
  nw{'AKs'}.standard_name = ncchar('ocean_vertical_salt_diffusivity');
  nw{'AKs'}.standard_name = 'ocean_vertical_salt_diffusivity';
  nw{'AKs'}.units = ncchar('m2 s-1');
  nw{'AKs'}.units = 'm2 s-1';
  nw{'AKs'}.comment = ncchar('Vertical component of the diffusivity of salt due to motion which is not resolved on the grid scale of the model');
  nw{'AKs'}.comment = 'Vertical component of the diffusivity of salt due to motion which is not resolved on the grid scale of the model';
  nw{'AKs'}.time = ncchar('scrum_time');
  nw{'AKs'}.time = 'scrum_time';
  nw{'AKs'}.coordinates = ncchar('lon_rho lat_rho s_w scrum_time');
  nw{'AKs'}.coordinates = 'lon_rho lat_rho s_w scrum_time';
  nw{'AKs'}.mask = ncchar('mask_rho');
  nw{'AKs'}.mask = 'mask_rho';
  nw{'AKs'}.field = ncchar('AKs, scalar, series');
  nw{'AKs'}.field = 'AKs, scalar, series';
 elseif(strcmp(namevar,'AKt'))
  nw{'AKt'} = ncdouble('scrum_time', 's_w', 'eta_rho', 'xi_rho');
  nw{'AKt'}.long_name = ncchar('Vertical diffusivity for heat');
  nw{'AKt'}.long_name = 'Vertical diffusivity for heat';
  nw{'AKt'}.standard_name = ncchar('ocean_vertical_heat_diffusivity');
  nw{'AKt'}.standard_name = 'ocean_vertical_heat_diffusivity';
  nw{'AKt'}.units = ncchar('m2 s-1');
  nw{'AKt'}.units = 'm2 s-1';
  nw{'AKt'}.comment = ncchar('Vertical component of the diffusivity of heat due to motion which is not resolved on the grid scale of the model');
  nw{'AKt'}.comment = 'Vertical component of the diffusivity of heat due to motion which is not resolved on the grid scale of the model';
  nw{'AKt'}.time = ncchar('scrum_time');
  nw{'AKt'}.time = 'scrum_time';
  nw{'AKt'}.coordinates = ncchar('lon_rho lat_rho s_w scrum_time');
  nw{'AKt'}.coordinates = 'lon_rho lat_rho s_w scrum_time';
  nw{'AKt'}.mask = ncchar('mask_rho');
  nw{'AKt'}.mask = 'mask_rho';
  nw{'AKt'}.field = ncchar('AKt, scalar, series');
  nw{'AKt'}.field = 'AKt, scalar, series';
 end
end % For loop
result = endef(nw);

%
% Create global attributes
%
nw.type = ncold.type(:);
nw.title = 'Raia configuration'; %nw.title = ncold.title(:);
nw.long_title = 'This a ROMS_AGRIF ouput file transformed to be CF-Compliant';
nw.comments = 'No comments';
nw.institution = 'IEO, A Corunha';
nw.source = 'ROMS_AGRIF';
nw.Conventions = 'CF-1.4';
nw.Conventions_help = 'http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/cf-conventions.html'; % This way I can always get back to fundamentals
nw.CreationDate = datestr(now,'yyyy-mm-dd HH:MM:SS');
nw.CreatedBy = getenv('LOGNAME');
nw.MatlabSource = version; % I like this one for tracking

%
% Write variables
%
disp(['Writing variables concerning grid structure']);
nw{'xl'}(:) = ncgrd{'xl'}(:);
nw{'el'}(:) = ncgrd{'el'}(:);
nw{'theta_s'}(:) = ncold.theta_s(:);
nw{'theta_b'}(:) = ncold.theta_b(:);
nw{'Tcline'}(:) = ncold.Tcline(:);
nw{'hc'}(:) = ncold.hc(:);
nw{'s_rho'}(:) = ncold.sc_r(:);
nw{'s_w'}(:) = ncold.sc_w(:);
nw{'Cs_r'}(:) = ncold.Cs_r(:);
nw{'Cs_w'}(:) = ncold.Cs_w(:);

nw{'h'}(:) = ncold{'h'}(:);
nw{'lon_rho'}(:) = ncgrd{'lon_rho'}(:);
nw{'lon_u'}(:) = ncgrd{'lon_u'}(:);
nw{'lon_v'}(:) = ncgrd{'lon_v'}(:);
%nw{'lon_psi'}(:) = ncgrd{'lon_psi'}(:);
nw{'lat_rho'}(:) = ncgrd{'lat_rho'}(:);
nw{'lat_u'}(:) = ncgrd{'lat_u'}(:);
nw{'lat_v'}(:) = ncgrd{'lat_v'}(:);
%nw{'lat_psi'}(:) = ncgrd{'lat_psi'}(:);
nw{'mask_rho'}(:) = ncgrd{'mask_rho'}(:);
nw{'mask_u'}(:) = ncgrd{'mask_u'}(:);
nw{'mask_v'}(:) = ncgrd{'mask_v'}(:);
%nw{'mask_psi'}(:) = ncgrd{'mask_psi'}(:);

nw{'scrum_time'}(:) = ncold{'scrum_time'}(:);

disp(['Writing variables concerning model results']);
for count=1:length(var(ncold))
 namevar=name(var(ncold,count));
 if(strcmp(namevar,'zeta'))
    kk=ncold{'zeta'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'zeta'}(:) = kk;
    fillval(nw{'zeta'}, theNewFillValue);
 elseif(strcmp(namevar,'temp'))
    kk=ncold{'temp'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'temp'}(:) = kk;
    fillval(nw{'temp'}, theNewFillValue);
 elseif(strcmp(namevar,'salt'))
    kk=ncold{'salt'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'salt'}(:) = kk;
    fillval(nw{'salt'}, theNewFillValue);
 elseif(strcmp(namevar,'u'))
    kk=ncold{'u'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'u'}(:) = kk;
    fillval(nw{'u'}, theNewFillValue);
 elseif(strcmp(namevar,'v'))
    kk=ncold{'v'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'v'}(:) = kk;
    fillval(nw{'v'}, theNewFillValue);
 elseif(strcmp(namevar,'w'))
    kk=ncold{'w'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'w'}(:) = kk;
    fillval(nw{'w'}, theNewFillValue);
 elseif(strcmp(namevar,'omega'))
    kk=ncold{'omega'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'omega'}(:) = kk;
    fillval(nw{'omega'}, theNewFillValue);
 elseif(strcmp(namevar,'ubar'))
    kk=ncold{'ubar'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'ubar'}(:) = kk;
    fillval(nw{'ubar'}, theNewFillValue);
 elseif(strcmp(namevar,'vbar'))
    kk=ncold{'vbar'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'vbar'}(:) = kk;
    fillval(nw{'vbar'}, theNewFillValue);
 elseif(strcmp(namevar,'AKv'))
    kk=ncold{'AKv'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'AKv'}(:) = kk;
    fillval(nw{'AKv'}, theNewFillValue);
 elseif(strcmp(namevar,'AKs'))
    kk=ncold{'AKs'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'AKs'}(:) = kk;
    fillval(nw{'AKs'}, theNewFillValue);
 elseif(strcmp(namevar,'AKt'))
    kk=ncold{'AKt'}(:); isee=find(kk==0); kk(isee)=theNewFillValue;
    nw{'AKt'}(:) = kk;
    fillval(nw{'AKt'}, theNewFillValue);
 end
end


close(nw);





