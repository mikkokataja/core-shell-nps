% A script for calculating the absorption and reflection coefficients of an ensemble of core-shell nanoparticles
% These models only work for particles that are much smaller than the
% incident wavelength (less than 50 nm diameter for visible light)
% 25.3.2020 Mikko Kataja

% Details of the particles

r_core = 10e-9; %Core radius
r_shell = 20e-9; %Shell radius

core_n = 'Fe3O4data.mat'; %should be in matlab file with separate variables for wl (in nm) and corresponing n and k
shell_n = 'Audata.mat'; %should be in matlab file with separate variables for wl (in nm) and corresponing n and k

wavelengths = 350:10:700; %wavelengths to calculate (in nm)
n_env = 1; %refractive index of the surrounding medium (1 = air, ~1.52 = glass)

% Number of particles and filling fraction

N_particles = 1e10;
device_area = 1e-4; %1e-4 = 1 cm^2

area_per_particle = device_area/N_particles;

%Load materialdata and calculate permittivity

load(core_n);
eps_core = (interp1(wl,n,wavelengths)+1i*interp1(wl,k,wavelengths)).^2;
load(shell_n);
eps_shell = (interp1(wl,n,wavelengths)+1i*interp1(wl,k,wavelengths)).^2;
eps_env = (n_env)^2;


% Calculate effective permittivity of the core-shell NP
% For more information see: 
% Uday K. Chettiar and Nader Engheta: Internal homogenization: Effective permittivity of a coated sphere 
% Optics Express Vol. 20, Issue 21, pp. 22976-22986 (2012)

eps_eff = eps_shell.*(((r_shell^3)*(eps_core+2*eps_shell)+2*(r_core^3)*(eps_core-eps_shell))./((r_shell^3)*(eps_core+2*eps_shell)-(r_core^3)*(eps_core-eps_shell)));

% Calculate polarizability (simple quasi-static approximation see e.g.
% Bohren & Huffman: Absorption and scattering of light by small particles
% The assumption here is that the particle is much smaller than wavelength

alpha = 4*pi*(r_shell^3)*((eps_eff-eps_env)./(eps_eff+2*eps_env));

%Calculate scattering and absorption cross sections (in square nm)

k = 2*pi./(wavelengths*1e-9);

c_abs = k.*imag(alpha);
c_sca = ((k.^4)./(6*pi)).*(abs(alpha).^2);

ext = (c_abs+c_sca)./area_per_particle; %extinction (from 1 to 0)
ref = c_sca/area_per_particle; %reflectivity (from 1 to 0)

% Plot the results
% Note that if you have too many particles (less area per particle than the
% extinction cross section) you will have transmission < 0.
figure; plot(wavelengths,1-ext)
xlabel('Wavelength (nm)')
ylabel('Transmission')

