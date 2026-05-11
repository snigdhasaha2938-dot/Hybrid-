%% ========================================================================
%  HYBRID ELECTROSTATIC-AERODYNAMIC FILTRATION BARRIER SIMULATION
%  For Urban PM2.5 and Virus-Laden Bioaerosol Mitigation
%
%  Research-Grade Multiphysics Simulation
%  Author: Computational Physics Engineer
%  Date: 2026-05-10
%
%  Citation Format:
%  "Hybrid Electrostatic–Aerodynamic Filtration Barrier for Urban PM2.5
%   and Virus-Laden Bioaerosol Mitigation" (Simulation Model v1.0)
% ========================================================================

clear all; close all; clc;
warning('off', 'MATLAB:imfilter:warnUnrecognized');

%% ========================================================================
% PHASE 1: INITIALIZE COMPUTATIONAL DOMAIN & PHYSICAL PARAMETERS
% ========================================================================

fprintf('\n');
fprintf('='*ones(1,70) + '\n');
fprintf('HYBRID ELECTROSTATIC-AERODYNAMIC FILTRATION SIMULATION v1.0\n');
fprintf('='*ones(1,70) + '\n');

% === DOMAIN GEOMETRY ===
domain.Lx = 0.5;      % Channel length [m]
domain.Ly = 0.2;      % Channel height [m] (y-direction: 0 to 0.2m)
domain.Lz = 0.15;     % Channel depth [m] (z-direction: 0 to 0.15m)

% === GRID RESOLUTION ===
domain.nx = 30;       % Points along x (flow direction)
domain.ny = 25;       % Points along y (vertical)
domain.nz = 20;       % Points along z (depth)

% Spacing
domain.dx = domain.Lx / domain.nx;
domain.dy = domain.Ly / domain.ny;
domain.dz = domain.Lz / domain.nz;

% Coordinate arrays (3D meshgrid)
[X, Y, Z] = meshgrid(linspace(0, domain.Lx, domain.nx), ...
linspace(0, domain.Ly, domain.ny), ...
linspace(0, domain.Lz, domain.nz));
X = permute(X, [2, 1, 3]);
Y = permute(Y, [2, 1, 3]);
Z = permute(Z, [2, 1, 3]);

fprintf('Domain initialized: %d x %d x %d grid\n', domain.nx, domain.ny, domain.nz);
fprintf('  Lx = %.3f m,  Ly = %.3f m,  Lz = %.3f m\n', domain.Lx, domain.Ly, domain.Lz);

% === PHYSICAL CONSTANTS ===
constants.rho_air = 1.225;        % Air density [kg/m³]
constants.mu_air = 1.81e-5;       % Dynamic viscosity [Pa·s]
constants.k_B = 1.38e-23;         % Boltzmann constant [J/K]
constants.T = 293;                % Temperature [K] (20°C)
constants.R_air = 287;            % Gas constant for air [J/(kg·K)]

% PM2.5 particle properties
particle.d_mean_pm25 = 1.0e-6;    % Mean diameter [m] (1 μm)
particle.d_std = 0.3e-6;          % Std deviation [m]
particle.rho_particle = 1500;     % Particle density [kg/m³] (dust)

% Viral aerosol properties (RSV, influenza-like)
viral.d_mean = 0.1e-6;            % Mean diameter [m] (0.1 μm)
viral.d_std = 0.05e-6;            % Std deviation [m]
viral.rho_virus = 1400;           % Viral particle density [kg/m³]
viral.charge_single = 1.6e-19;    % Single elementary charge [C]

fprintf('Physical constants set:\n');
fprintf('  ρ_air = %.3f kg/m³,  μ = %.2e Pa·s\n', constants.rho_air, constants.mu_air);
fprintf('  T = %.0f K (%.1f°C)\n', constants.T, constants.T - 273.15);

% === AIRFLOW PARAMETERS ===
airflow.U_inlet = 1.5;            % Inlet velocity [m/s]
airflow.turbulence_intensity = 0.08;  % Turbulence intensity [%]
Re = (constants.rho_air * airflow.U_inlet * domain.Ly) / constants.mu_air;
fprintf('Airflow configuration:\n');
fprintf('  U_inlet = %.2f m/s,  Re (channel) = %.1f\n', airflow.U_inlet, Re);

% === ELECTROSTATIC PARAMETERS ===
electrostatics.E_field = 5e5;     % Electric field strength [V/m]
electrostatics.V_total = electrostatics.E_field * domain.Ly;  % Total voltage
electrostatics.electrode_pos_z = 0;  % Electrode positions (simplified)

% Electrostatic number (Ne = qE / (3πμdU))
Ne_pm25 = (1 * electrostatics.E_field) / (3 * pi * constants.mu_air * particle.d_mean_pm25 * airflow.U_inlet);
fprintf('Electrostatic field:\n');
fprintf('  E = %.2e V/m,  V_total = %.1f kV\n', electrostatics.E_field, electrostatics.V_total/1000);
fprintf('  Ne (PM2.5) = %.2f\n', Ne_pm25);

% === SIMULATION TIME PARAMETERS ===
sim.t_end = 2.0;                  % Total simulation time [s]
sim.dt = 0.005;                   % Time step [s]
sim.n_steps = ceil(sim.t_end / sim.dt);
sim.save_freq = 20;               % Save state every N steps

fprintf('Simulation time: %.2f s,  dt = %.4f s,  steps = %d\n', ...
sim.t_end, sim.dt, sim.n_steps);

% === PARTICLE INITIALIZATION ===
n_particles_pm25 = 800;           % Number of PM2.5 particles
n_particles_viral = 300;          % Number of viral particles
n_particles_total = n_particles_pm25 + n_particles_viral;

fprintf('\nInitializing particle populations:\n');
fprintf('  PM2.5 particles: %d\n', n_particles_pm25);
fprintf('  Viral particles: %d\n', n_particles_viral);
fprintf('  Total: %d\n', n_particles_total);

%% ========================================================================
% PHASE 2: INITIALIZE PARTICLE STATE
% ========================================================================
% State: [x, y, z, vx, vy, vz, d, q, m, type, age, is_captured]
%   type: 1=PM2.5, 2=Viral, 3=Aggregated complex

particles = initialize_particles(n_particles_pm25, n_particles_viral, ...
domain, particle, viral, constants);

fprintf('Particle state matrix: %d x 12\n', size(particles, 1));
fprintf('  Columns: [x, y, z, vx, vy, vz, d, q, m, type, age, is_captured]\n');

%% ========================================================================
% PHASE 3: COMPUTE FLOW AND ELECTROSTATIC FIELDS
% ========================================================================

fprintf('\n--- Computing airflow field...\n');
[U_field, V_field, W_field] = airflow_solver(X, Y, Z, airflow, domain);
fprintf('Airflow field computed. Mean U = %.3f m/s\n', mean(U_field(:)));

fprintf('--- Computing electrostatic field...\n');
[E_x, E_y, E_z, Phi] = electrostatic_field(X, Y, Z, electrostatics, domain);
fprintf('Electrostatic field computed. E_z range: [%.2e, %.2e] V/m\n', ...
min(E_z(:)), max(E_z(:)));

% === DIMENSIONLESS NUMBERS ANALYSIS ===
fprintf('\n--- Dimensionless Analysis ---\n');

% Stokes number: Stk = (ρ_p * d_p²) / (18 * μ * t_scale)
t_scale = domain.Lx / airflow.U_inlet;
Stokes_pm25 = (particle.rho_particle * particle.d_mean_pm25^2) / (18 * constants.mu_air * t_scale);
Stokes_viral = (viral.rho_virus * viral.d_mean^2) / (18 * constants.mu_air * t_scale);

fprintf('Stokes number (residence time scale = %.3f s):\n', t_scale);
fprintf('  Stk_PM2.5 = %.4f (diffusive regime)\n', Stokes_pm25);
fprintf('  Stk_viral = %.4f (Brownian regime)\n', Stokes_viral);

% Peclet number: Pe = U*L / D_brownian
D_brownian_pm25 = (constants.k_B * constants.T) / (3 * pi * constants.mu_air * particle.d_mean_pm25);
D_brownian_viral = (constants.k_B * constants.T) / (3 * pi * constants.mu_air * viral.d_mean);
Pe_pm25 = (airflow.U_inlet * domain.Lx) / D_brownian_pm25;
Pe_viral = (airflow.U_inlet * domain.Lx) / D_brownian_viral;

fprintf('Peclet number:\n');
fprintf('  Pe_PM2.5 = %.2f (advection-dominated)\n', Pe_pm25);
fprintf('  Pe_viral = %.2f (Brownian-enhanced)\n', Pe_viral);

fprintf('Brownian diffusivity:\n');
fprintf('  D_PM2.5 = %.2e m²/s\n', D_brownian_pm25);
fprintf('  D_viral = %.2e m²/s\n', D_brownian_viral);

%% ========================================================================
% PHASE 4: TIME-STEPPING SIMULATION
% ========================================================================

fprintf('\n');
fprintf('='*ones(1,70) + '\n');
fprintf('BEGINNING TIME-STEPPING SIMULATION\n');
fprintf('='*ones(1,70) + '\n');

% Storage for history
history.time = [];
history.n_captured = [];
history.n_aggregated = [];
history.mean_velocity = [];
history.n_particles = [];
history.efficiency = [];

% Flags
capture_zone = (Y >= 0.15) & (Y <= 0.2);  % Top region for capture

for step = 1:sim.n_steps

t_current = step * sim.dt;  
  
% --- UPDATE PARTICLE DYNAMICS ---  
particles = particle_dynamics(particles, U_field, V_field, W_field, ...  
    E_x, E_y, E_z, X, Y, Z, domain, constants, airflow, sim);  
  
% --- BROWNIAN MOTION (Weak random walk) ---  
particles = apply_brownian_motion(particles, constants, sim.dt, domain);  
  
% --- PARTICLE AGGREGATION ---  
particles = particle_aggregation(particles, domain, constants);  
  
% --- CAPTURE DETECTION ---  
particles = detect_capture(particles, domain, capture_zone);  
  
% --- REMOVE ESCAPED PARTICLES ---  
particles = remove_escaped(particles, domain);  
  
% --- COLLECT STATISTICS ---  
n_escaped = sum(particles(:,1) == -1);  
n_active = size(particles, 1) - n_escaped;  
n_captured = sum(particles(:, 12) == 1);  
n_aggregated = sum(particles(:, 10) == 3);  
  
if mod(step, 50) == 0  
    fprintf('Step %5d / %5d | t = %.3f s | Active: %d | Captured: %d | Aggregated: %d\n', ...  
        step, sim.n_steps, t_current, n_active, n_captured, n_aggregated);  
end  
  
% Store history  
history.time = [history.time; t_current];  
history.n_captured = [history.n_captured; n_captured];  
history.n_aggregated = [history.n_aggregated; n_aggregated];  
history.n_particles = [history.n_particles; n_active];  
  
if n_active > 0  
    v_magnitude = sqrt(particles(1:n_active, 4).^2 + particles(1:n_active, 5).^2 + particles(1:n_active, 6).^2);  
    history.mean_velocity = [history.mean_velocity; mean(v_magnitude)];  
else  
    history.mean_velocity = [history.mean_velocity; 0];  
end  
  
% Compute efficiency  
if n_particles_total > 0  
    efficiency = (n_captured / n_particles_total) * 100;  
    history.efficiency = [history.efficiency; efficiency];  
end

end

fprintf('\n');
fprintf('='*ones(1,70) + '\n');
fprintf('SIMULATION COMPLETED\n');
fprintf('='*ones(1,70) + '\n');

%% ========================================================================
% PHASE 5: POST-PROCESSING & ANALYSIS
% ========================================================================

fprintf('\nPost-processing results...\n');

% Final statistics
n_escaped_final = sum(particles(:,1) == -1);
n_active_final = size(particles, 1) - n_escaped_final;
n_captured_final = sum(particles(:, 12) == 1);
efficiency_final = (n_captured_final / n_particles_total) * 100;

fprintf('\n--- FINAL RESULTS ---\n');
fprintf('Particles captured: %d / %d (%.1f%%)\n', n_captured_final, n_particles_total, efficiency_final);
fprintf('Particles escaped: %d / %d\n', n_particles_total - n_captured_final, n_particles_total);
fprintf('Active in domain: %d\n', n_active_final);

% Separate by type
pm25_captured = sum((particles(:, 12) == 1) & (particles(:, 10) == 1));
viral_captured = sum((particles(:, 12) == 1) & (particles(:, 10) == 2));
fprintf('\nBy particle type:\n');
fprintf('  PM2.5 captured: %d / %d (%.1f%%)\n', pm25_captured, n_particles_pm25, ...
(pm25_captured/n_particles_pm25)*100);
fprintf('  Viral captured: %d / %d (%.1f%%)\n', viral_captured, n_particles_viral, ...
(viral_captured/n_particles_viral)*100);

% Mean residence time
active_idx = particles(:,1) ~= -1;
if sum(active_idx) > 0
residence_times = particles(active_idx, 11);  % age column
fprintf('\nResidence time statistics:\n');
fprintf('  Mean: %.3f s,  Max: %.3f s\n', mean(residence_times), max(residence_times));
end

%% ========================================================================
% PHASE 6: VISUALIZATION SUITE
% ========================================================================

fprintf('\nGenerating publication-quality figures...\n');

% Figure 1: AIRFLOW FIELD WITH STREAMLINES
fig1 = create_airflow_visualization(X, Y, Z, U_field, V_field, W_field, domain);

% Figure 2: ELECTROSTATIC FIELD HEATMAP
fig2 = create_electrostatic_visualization(X, Y, Z, E_z, Phi, domain);

% Figure 3: PARTICLE TRAJECTORIES
fig3 = create_particle_trajectory_visualization(particles, domain);

% Figure 4: FILTRATION EFFICIENCY TIME SERIES
fig4 = create_efficiency_plot(history);

% Figure 5: PARTICLE SIZE DISTRIBUTION
fig5 = create_size_distribution_plot(particles);

% Figure 6: CAPTURE ZONE HEATMAP
fig6 = create_capture_zone_heatmap(particles, domain);

% Figure 7: VELOCITY FIELD VECTOR MAP
fig7 = create_velocity_vector_field(X, Y, Z, U_field, V_field, W_field, domain);

% Figure 8: BEFORE/AFTER PARTICLE CONCENTRATION
fig8 = create_concentration_comparison(particles, n_particles_total, domain);

fprintf('All visualizations completed.\n');

%% ========================================================================
% PHASE 7: SAVE RESULTS
% ========================================================================

fprintf('\nSaving results...\n');

% Create results structure
results.domain = domain;
results.airflow = airflow;
results.electrostatics = electrostatics;
results.particles_final = particles;
results.history = history;
results.efficiency_final = efficiency_final;
results.constants = constants;

% Save to MAT file
save('filtration_simulation_results.mat', 'results', '-v7.3');
fprintf('Results saved to: filtration_simulation_results.mat\n');

% Generate summary report
generate_summary_report(results, history, efficiency_final, ...
pm25_captured, viral_captured, n_particles_total, constants);

fprintf('\n');
fprintf('='*ones(1,70) + '\n');
fprintf('SIMULATION PIPELINE COMPLETE\n');
fprintf('='*ones(1,70) + '\n\n');

%% ========================================================================
% END OF MAIN SIMULATION SCRIPT
% ========================================================================
