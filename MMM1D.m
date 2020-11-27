% ONE-DIMENSIONAL STEADY-STATE NON-ISOTHERMAL MACRO-HOMOGENEOUS TWO-PHASE
% MASTER MEA MODEL FOR PROTON EXCHANGE MEMBRANE FUEL CELLS
%
% Authors:
% Version 2:
% Dr. Robert Herrendörfer, robert.herrendoerfer@zhaw.ch
% Prof. Dr. Juergen O. Schumacher, juergen.schumacher@zhaw.ch

% Version 1:
% Dr. Roman Vetter, roman.vetter@bsse.ethz.ch
% Prof. Dr. Juergen O. Schumacher, juergen.schumacher@zhaw.ch
%
% Copyright (c) 2017-2021 ZHAW-ICP
% Institute of Computational Physics ICP
% Zurich University of Applied Sciences ZHAW
% CH-8401 Winterthur, Switzerland
%
% This software is subject to the accompanying license agreement.

function [I,U, SOL,x, Lsum, Np, Neq, domains] = MMM1D

% CONSTANTS
M_H2O = 18e-3; % [kg/mol] molar mass of water
M_H2  = 2e-3;  % [kg/mol] molar mass of hydrogen
M_N2  = 28e-3; % [kg/mol] molar mass of nitrogen 
M_O2  = 32e-3; % [kg/mol] molar mass of nitrogen
T_0 = 273.15; % [K] zero degrees Celsius
R = 8.31446; % [J*mol/K] universal gas constant
F = 96485.333; % [C/mol] Faraday constant
P_ref = 101325; % [Pa] reference pressure
T_ref = T_0+80; % [K] reference temperature

% OPERATING CONDITIONS
U = 1.15:-0.05:0.00; % [V] list of cell voltages
P_A = 1.5e5; % [Pa] total pressure in anode gas channel
P_C = 1.5e5; % [Pa] total pressure in cathode gas channel
RH_A = 0.9; % [-] relative humidity in anode gas channel
RH_C = 0.9; % [-] relative humidity in cathode gas channel
s_C = 0.12; % [-] liquid water saturation at cathode GDL/GC interface
T_A = T_0+70; % [K] temperature of anode bipolar plate and gas channel
T_C = T_0+70; % [K] temperature of cathode bipolar plate and gas channel
alpha_O2 = 0.21; % [-] mole fraction of oxygen in dry oxidant gas

% ELECTROCHEMICAL PARAMETERS
A = @(E,T) exp(E/R*(1/T_ref-1./T)); % [-] Arrhenius correction
i_0_HOR = @(T) 0.27e4*A(16e3,T); % [A/m^2] exchange current density of HOR
i_0_ORR = @(T,P_O2) 2.47e-4*(P_O2/P_ref).^0.54.*A(67e3,T); % [A/m^2] exchange current density of ORR
beta_HOR = 0.5; % [-] HOR symmetry factor
beta_ORR = 0.5; % [-] ORR symmetry factor
DeltaH = -285.83e3; % [J/mol] enthalpy of formation of liquid water
DeltaS_HOR = 0.104; % [J/(mol*K)] reaction entropy of HOR
DeltaS_ORR = -163.3; % [J/(mol*K)] reaction entropy of ORR

% MATERIAL PARAMETERS
L = [160 10 25 10 160]*1e-6; % [m] MEA layer thicknesses
a_ACL = 1e7; % [1/m] ECSA density of ACL
a_CCL = 3e7; % [1/m] ECSA density of CCL
H_ec = 42e3; % [J/mol] molar enthalphy of evaporation/condensation
H_ad = H_ec; % [J/mol] molar enthalphy of absorption/desorption
k_GDL = 1.6; % [W/(m*K)] thermal conductivity of GDL
k_CL = 0.27; % [W/(m*K)] thermal conductivity of CL
k_PEM = 0.3; % [W/(m*K)] thermal conductivity of PEM
s_im = s_C; % [-] immobile liquid water saturation
V_m = 1.02/1.97e3; % [m^3/mol] molar volume of dry membrane (equivalent weight divided by mass density)
V_w = M_H2O/0.978e3; % [m^3/mol] molar volume of liquid water (molar mass divided by mass density)
mu_gas = 2.1e-5;  % [Pa*s] gas density 
eps_i_CL = 0.3; % [-] volume fraction of ionomer in dry CL
eps_p_GDL = 0.76; % [-] porosity of GDL
eps_p_CL = 0.4; % [-] porosity of CL
kappa_GDL = 6.15e-12; % [m^2] absolute permeability of GDL
kappa_CL = 1e-13; % [m^2] absolute permeability of CL
sigma_e_GDL = 1250; % [S/m] electrical conductivity of GDL
sigma_e_CL = 350; % [S/m] electrical conductivity of CL
tau_GDL = 1.6; % [-] pore tortuosity of GDL
tau_CL = 1.6; % [-] pore tortuosity of CL

% WATER CONSTITUTIVE RELATIONSHIPS
P_sat = @(T) exp(23.1963-3816.44./(T-46.13)); % [Pa] saturation pressure of water vapor
mu = @(T) 1e-3*exp(-3.63148+542.05./(T-144.15)); % [Pa*s] dynamic viscosity of liquid water

% MODEL PARAMETERIZATION
BV = @(i_0,a,T,beta,eta) i_0*a.*(exp(beta*2*F./(R*T).*eta)-exp(-(1-beta)*2*F./(R*T).*eta)); % [A/m^3] Butler-Volmer eq
D = @(eps_p,tau,s,P,T) eps_p/tau^2*(1-s).^3.*(T/T_ref).^1.5.*(P_ref./P); % [-] scaling factor for gas diffusivities
D_O2 = @(eps_p,tau,s,T,P) 0.28e-4*D(eps_p,tau,s,P,T); % [m^2/s] oxygen diffusion coefficient
D_H2O_A = @(eps_p,tau,s,T,P) 1.24e-4*D(eps_p,tau,s,P,T); % [m^2/s] water vapor diffusion coefficient at anode
D_H2O_C = @(eps_p,tau,s,T,P) 0.36e-4*D(eps_p,tau,s,P,T); % [m^2/s] water vapor diffusion coefficient at cathode
f = @(lambda) lambda*V_w./(V_m+lambda*V_w); % [-] volume fraction of water in ionomer
x_H2O_A = RH_A*P_sat(T_A)/P_A; % [-] mole fraction of water vapor in anode gas channel
w_H2O_A = 1/(M_H2/M_H2O*(1/x_H2O_A-1)+1); % [-] mass fraction of water vapor in anode gas channel
x_H2O_C = RH_C*P_sat(T_C)/P_C; % [-] mole fraction of water vapor in cathode gas channel
x_O2_C = alpha_O2*(1-x_H2O_C); % [-] mole fraction of oxygen in cathode gas channel
w_H2O_C = x_H2O_C*M_H2O/(M_N2-x_O2_C*(M_N2-M_O2)-x_H2O_C*(M_N2-M_H2O)); % [-] mass fraction of oxygen in anode gas channel
w_O2_C = x_O2_C*M_O2/(x_H2O_C*M_H2O)*w_H2O_C;

% AUXILIARY FUNCTIONS
iff = @(cond,a,b) cond.*a + ~cond.*b; % vectorized ternary operator

% MATERIAL CONSTITUTIVE RELATIONSHIPS
k_ad = @(lambda,lambda_eq,T) iff(lambda<lambda_eq,3.53e-5,1.42e-4).*f(lambda).*A(20e3,T); % [m/s] mass transfer coefficient for vapor sorption of Nafion
s_red = @(s) (s-s_im)/(1-s_im); % reduced liquid water saturation
gamma_ec = @(x_H2O,x_sat,s,T) 2e6*iff(x_H2O<x_sat,5e-4*s_red(s),6e-3*(1-s_red(s))).*sqrt(R*T/(2*pi*M_H2O)); % [1/s] evaporation/condensation rate
sigma_p = @(eps_i,lambda,T) eps_i^1.5*116*max(0,f(lambda)-0.06).^1.5.*A(15e3,T); % [S/m] proton conductivity of Nafion
sorption = @(RH) 0.043+17.81*RH-39.85*RH.^2+36.0*RH.^3; % [-] vapor sorption isotherm of Nafion
D_lambda = @(eps_i,lambda,T) eps_i^1.5*(3.842*lambda.^3-32.03*lambda.^2+67.74*lambda)./(lambda.^3-2.115*lambda.^2-33.013*lambda+103.37)*1e-10.*A(20e3,T); % [m^2/s] diffusion coefficient of water in Nafion
xi = @(lambda) 2.5*lambda/22; % [-] electro-osmotic drag coefficient of Nafion
dpds = @(s) 0.00011*44.02*exp(-44.02*(s-0.496))+278.3*8.103*exp(8.103*(s-0.496)); % [Pa] derivative of capillary pressure-saturation relationship of GDL
D_s = @(kappa,s,T) kappa*(1e-6+s_red(s).^3)./mu(T).*dpds(s); % [m^2/s] liquid water transport coefficient

% INITIAL MESH
Lsum = [0 cumsum(L)];
Nd = numel(L); % number of domains
x = interp1(0:Nd, Lsum, linspace(0, Nd, Nd*10+1));
x = sort([x Lsum(2:end-1)]); % duplicate interface nodes

% SOLVER PREPARATION
sol = bvpinit(x, @yinit);
options = bvpset('Vectorized', 'on', 'NMax', 1e5, 'RelTol', 1e-4, 'AbsTol', 1e-6);

% PARAMETER SWEEP
I = zeros(size(U));
SOL = cell(size(U));
Np = numel(U); % number of parameters in the sweep
Neq = size(sol.y,1)/2; % number of 2nd-order differential equations
for k = 1:Np
    sol = bvp4c(@odefun, @(ya,yb) bcfun(ya, yb, U(k)), sol, options);
    SOL{k} = sol;
    I(k) = sol.y(2,end)/1e4; % current density in [A/cm^2]
end

% POSTPROCESSING
Nref = 2; % number of refinements for smoother curve plotting
domains = [1 1 0 1 1; % domains for phi_e
           0 1 1 1 0; % domains for phi_p
           1 1 1 1 1; % domains for T
           0 1 1 1 0; % domains for lambda
           1 1 0 1 1; % domains for w_H2O
           0 0 0 1 1; % domains for w_O2
           0 0 0 1 1; % domains for s
           1 1 0 1 1];% domains for u_gas
shift = 1e-10;
for k = 1:Np
    x = [];
    for m = 1:Nd
        xa = find(SOL{k}.x==Lsum(m  ), 1, 'last' );
        xb = find(SOL{k}.x==Lsum(m+1), 1, 'first');
        N = xb-xa;
        
        % grid refinement
        x = [x interp1(linspace(0,1,N+1), SOL{k}.x(xa:xb), linspace(shift, 1-shift, N*2^Nref+1))];
        
        % fill solution on inactive domains with NaN
        SOL{k}.y(~kron(domains(:,m),ones(2,1)),xa:xb) = NaN;
    end
    [SOL{k}.y, SOL{k}.yp] = deval(SOL{k}, x);
    SOL{k}.x = x;
end

% SYSTEM OF 1ST-ORDER DIFFERENTIAL EQUATIONS
function dydx = odefun(x, y, subdomain)

    % READ POTENTIALS & FLUXES
    phi_e  = y( 1,:); j_e      = y( 2,:);
    phi_p  = y( 3,:); j_p      = y( 4,:);
    T      = y( 5,:); j_T      = y( 6,:);
    lambda = y( 7,:); j_lambda = y( 8,:);
    w_H2O  = y( 9,:); j_H2O    = y(10,:);
    w_O2   = y(11,:); j_O2     = y(12,:);
    s      = y(13,:); j_s      = y(14,:);
    P_gas  = y(15,:); rho_u_gas= y(16,:);

    % ZERO-INITIALIZE ALL DERIVATIVES
    z = zeros(size(x));
    dphi_e  = z; dj_e      = z;
    dphi_p  = z; dj_p      = z;
    dT      = z; dj_T      = z;
    dlambda = z; dj_lambda = z;
    dw_H2O  = z; dj_H2O    = z;
    dw_O2   = z; dj_O2     = z;
    ds      = z; dj_s      = z;
    dP_gas  = z; drho_u_gas= z;

    % COMPUTE DERIVATIVES
    switch subdomain
        case 1 % AGDL
            w_H2=1-w_H2O;
            Mn=1./(w_H2O/M_H2O+w_H2/M_H2);
            C = P_gas./(R*T); % total interstitial gas concentration
            rho_gas=Mn.*C; % total gas density
            % Flux definitions
            u_gas=rho_u_gas./rho_gas;
            dP_gas=-u_gas./(kappa_GDL/mu_gas);
            dphi_e = -j_e/sigma_e_GDL; % electron flux: j_e = -sigma_e*grad(phi_e)
            dT = -j_T/k_GDL; % heat flux: j_T = -k*grad(T)
            dw_H2O = -(j_H2O-rho_gas.*w_H2O.*u_gas)./(rho_gas.*D_H2O_A(eps_p_GDL,tau_GDL,s,T,P_gas)); % water vapor flux: j_H2O = -rho*D_H2O*grad(w_H2O)
            dj_T = -j_e.*dphi_e; % conservation of heat: div(j_T) = S_T
        case 2 % ACL
            w_H2=1-w_H2O;
            Mn=1./(w_H2O/M_H2O+w_H2/M_H2);
            x_H2O=w_H2O/M_H2O.*Mn;
            x_H2=1-x_H2O;
            C = P_gas./(R*T); % total interstitial gas concentration
            rho_gas=Mn.*C; % total gas density
            x_sat = P_sat(T)./P_gas; % saturation water vapor mole fraction
            lambda_eq = sorption(x_H2O./x_sat); % equilibrium water content of ionomer
            % Source terms
            S_ad = k_ad(lambda,lambda_eq,T)/(L(2)*V_m).*(lambda_eq-lambda); % absorption/desorption reaction rate
            P_H2 = x_H2.*P_gas;
            eta = phi_e-phi_p+T*DeltaS_HOR/(2*F)+R*T/(2*F).*log(P_H2/P_ref); % overpotential
            i = BV(i_0_HOR(T),a_ACL,T,beta_HOR,eta); % electrochemical reaction rate
            S_F = i/(2*F); % Faraday's law
            % Flux definitions
            u_gas=rho_u_gas./rho_gas;
            dP_gas=-u_gas./(kappa_CL/mu_gas);
            dphi_e = -j_e/sigma_e_CL; % electron flux: j_e = -sigma_e*grad(phi_e)
            dphi_p = -j_p./sigma_p(eps_i_CL,lambda,T); % proton flux: j_p = -sigma_p*grad(phi_p)
            dT = -j_T/k_CL; % heat flux: j_T = -k*grad(T)
            dlambda = (-j_lambda+xi(lambda)/F.*j_p)*V_m./D_lambda(eps_i_CL,lambda,T); % dissolved water flux: j_lambda = -D_lambda/V_m*grad(lambda)+xi/F*j_p
            dw_H2O = -(j_H2O-rho_gas.*w_H2O.*u_gas)./(rho_gas.*D_H2O_A(eps_p_CL,tau_CL,s,T,P_gas)); % water vapor flux: j_H2O = -C*D_H2O*grad(w_H2O)
            dj_e = -i; % conservation of electrons: div(j_e) = S_e
            dj_p = i; % conservation of protons: div(j_p) = S_p
            dj_T = -j_e.*dphi_e-j_p.*dphi_p+i.*eta-S_F.*T*DeltaS_HOR+H_ad*S_ad; % conservation of heat: div(j_T) = S_T
            dj_lambda = S_ad; % conservation of dissolved water: div(j_lambda) = S_lambda
            sumS=-M_H2O*S_ad-M_H2*S_F; % sum of all source term of gas species
            dj_H2O = -S_ad*M_H2O; % conservation of water vapor: div(j_H2O) = S_H2O
            drho_u_gas=sumS;
        case 3 % PEM
            dphi_p = -j_p./sigma_p(1,lambda,T); % proton flux: j_p = -sigma_p*grad(phi_p)
            dT = -j_T/k_PEM; % heat flux: j_T = -k*grad(T)
            dlambda = (-j_lambda+xi(lambda)/F.*j_p)*V_m./D_lambda(1,lambda,T); % dissolved water flux: j_lambda = -D_lambda/V_m*grad(lambda)+xi/F*j_p
            dj_T = -j_p.*dphi_p; % conservation of heat: div(j_T) = S_T
        case 4 % CCL
            w_N2=1-w_H2O-w_O2;
            Mn=1./(w_H2O/M_H2O+w_O2/M_O2+w_N2/M_N2);
            x_H2O=w_H2O/M_H2O.*Mn;
            x_O2=w_O2/M_O2.*Mn;
            C = P_gas./(R*T); % total interstitial gas concentration
            rho_gas=Mn.*C; % total gas density
            x_sat = P_sat(T)./P_gas; % saturation water vapor mole fraction
            S_ec = gamma_ec(x_H2O,x_sat,s,T).*C.*(x_H2O-x_sat); % evaporation/condensation reaction rate
            lambda_eq = sorption(x_H2O./x_sat); % equilibrium water content of ionomer
            S_ad = k_ad(lambda,lambda_eq,T)/(L(4)*V_m).*(lambda_eq-lambda); % absorption/desorption reaction rate
            P_O2 = x_O2.*P_gas; % partial pressure of oxygen
            eta = -(DeltaH-T*DeltaS_ORR)/(2*F)+R*T/(4*F).*log(P_O2/P_ref)-(phi_e-phi_p); % overpotential
            i = BV(i_0_ORR(T,P_O2),a_CCL,T,beta_ORR,eta); % electrochemical reaction rate
            S_F = i/(2*F); % Faraday's law
            % Flux definitions
            u_gas=rho_u_gas./rho_gas;
            dP_gas=-u_gas./(kappa_CL/mu_gas);
            dphi_e = -j_e/sigma_e_CL; % electron flux: j_e = -sigma_e*grad(phi_e)
            dphi_p = -j_p./sigma_p(eps_i_CL,lambda,T); % proton flux: j_p = -sigma_p*grad(phi_p)
            dT = -j_T/k_CL; % heat flux: j_T = -k*grad(T)
            dlambda = (-j_lambda+xi(lambda)/F.*j_p)*V_m./D_lambda(eps_i_CL,lambda,T); % dissolved water flux: j_lambda = -D_lambda/V_m*grad(lambda)+xi/F*j_p
            dw_H2O = -(j_H2O-rho_gas.*w_H2O.*u_gas)./(rho_gas.*D_H2O_C(eps_p_CL,tau_CL,s,T,P_gas)); % water vapor flux: j_H2O = -C*D_H2O*grad(x_H2O)
            dw_O2 = -(j_O2-rho_gas.*w_O2.*u_gas)./(rho_gas.*D_O2(eps_p_CL,tau_CL,s,T,P_gas)); % oxygen flux: j_O2 = -C*D_O2*grad(x_O2)
            ds = -j_s*V_w./D_s(kappa_CL,s,T); % liquid water flux: j_s = -D_s/V_w*grad(s)
            dj_e = i; % conservation of electrons: div(j_e) = S_e
            dj_p = -i; % conservation of protons: div(j_p) = S_p
            dj_T = -j_e.*dphi_e-j_p.*dphi_p+i.*eta-S_F.*T*DeltaS_ORR+H_ad*S_ad+H_ec*S_ec; % conservation of heat: div(j_T) = S_T
            dj_lambda = S_F+S_ad; % conservation of dissolved water: div(j_lambda) = S_lambda
            sumS=-M_H2O*(S_ec+S_ad)-M_O2*S_F/2; % sum of all source term of gas species
            dj_H2O =-M_H2O*(S_ec+S_ad); % conservation of water vapor: div(j_H2O) = S_H2O
            dj_O2 =-M_O2*S_F/2; % conservation of oxygen: div(j_O2) = S_O2
            dj_s = S_ec; % conservation of liquid water: div(j_s) = S_s
            drho_u_gas=sumS;
        case 5 % CGDL
            w_N2=1-w_H2O-w_O2;
            Mn=1./(w_H2O/M_H2O+w_O2/M_O2+w_N2/M_N2);
            x_H2O=w_H2O/M_H2O.*Mn;
            C = P_gas./(R*T); % total interstitial gas concentration
            rho_gas=Mn.*C; % total gas density
            x_sat = P_sat(T)./P_gas; % saturation water vapor mole fraction
            S_ec = gamma_ec(x_H2O,x_sat,s,T).*C.*(x_H2O-x_sat); % evaporation/condensation reaction rate
            % Flux definitions
            u_gas=rho_u_gas./rho_gas;
            dP_gas=-u_gas./(kappa_GDL/mu_gas);
            dphi_e = -j_e/sigma_e_GDL; % electron flux: j_e = -sigma_e*grad(phi_e)
            dT = -j_T/k_GDL; % heat flux: j_T = -k*grad(T)
            dw_H2O = -(j_H2O-rho_gas.*w_H2O.*u_gas)./(rho_gas.*D_H2O_C(eps_p_GDL,tau_GDL,s,T,P_gas)); % water vapor flux: j_H2O = -C*D_H2O*grad(w_H2O)
            dw_O2 = -(j_O2-rho_gas.*w_O2.*u_gas)./(rho_gas.*D_O2(eps_p_GDL,tau_GDL,s,T,P_gas)); % oxygen flux: j_O2 = -C*D_O2*grad(w_O2)
            ds = -j_s*V_w./D_s(kappa_GDL,s,T); % liquid water flux: j_s = -D_s/V_w*grad(s)
            dj_T = -j_e.*dphi_e+H_ec*S_ec; % conservation of heat: div(j_T) = S_T
            sumS=-M_H2O*S_ec;
            dj_H2O = -M_H2O*S_ec; % conservation of water vapor: div(j_H2O) = S_H2O
            dj_s = S_ec; % conservation of liquid water: div(j_s) = S_s
            drho_u_gas=sumS;
    end
    % ASSEMBLE DERIVATIVES
    dydx = [dphi_e; dj_e;
            dphi_p; dj_p;
            dT; dj_T;
            dlambda; dj_lambda;
            dw_H2O; dj_H2O;
            dw_O2; dj_O2;
            ds; dj_s;
            dP_gas; drho_u_gas];
end

% INITIAL GUESS
function y0 = yinit(x, subdomain)

    % A ROUGH ESTIMATE OF THE POTENTIALS FOR THE FIRST CELL VOLTAGE
    phi_e = iff(subdomain > 3, U(1), 0);
    phi_p = 0;
    T = (T_C+T_A)/2;
    lambda = iff(subdomain > 1 && subdomain < 5, sorption(1), 0);
    w_H2O = iff(subdomain < 3, w_H2O_A, iff(subdomain > 3, w_H2O_C, 0));
    w_O2 = iff(subdomain > 3, w_O2_C, 0);
    s = iff(subdomain > 3, s_C, 0);
    P_gas = iff(subdomain < 3, P_A, iff(subdomain > 3, P_C, 0));
    
    % ALL FLUXES ARE INITIALLY ZERO
    y0 = [phi_e; 0; phi_p; 0; T; 0; lambda; 0; w_H2O; 0; w_O2; 0; s; 0; P_gas; 0];

end

% BOUNDARY & INTERFACE CONDITIONS
function res = bcfun(ya, yb, U)

    res = ya(:); % homogeneous BC everywhere by default

    % ELECTRONS
    res(0*Neq+1) = ya(1,1); % potential zeroing
    res(0*Neq+2) = ya(2,2) - yb(2,1); % flux continuity between AGDL & ACL
    res(2*Neq+1) = ya(1,2) - yb(1,1); % potential continuity between AGDL & ACL
    res(2*Neq+2) = yb(2,2); % zero flux between ACL & PEM
    res(6*Neq+1) = ya(1,5) - yb(1,4); % potential continuity between CCL & CGDL
    res(6*Neq+2) = ya(2,4); % zero flux between PEM & CCL
    res(8*Neq+1) = yb(1,5) - U; % cell voltage
    res(8*Neq+2) = ya(2,5) - yb(2,4); % flux continuity between CCL & CGDL

    % PROTONS
    res(2*Neq+3) = ya(4,4) - yb(4,3); % flux continuity between PEM & CCL
    res(2*Neq+4) = ya(4,2); % zero flux between AGDL & ACL
    res(4*Neq+3) = ya(3,3) - yb(3,2); % potential continuity between ACL & PEM
    res(4*Neq+4) = ya(4,3) - yb(4,2); % flux continuity between ACL & PEM
    res(6*Neq+3) = ya(3,4) - yb(3,3); % potential continuity between PEM & CCL
    res(6*Neq+4) = yb(4,4); % zero flux between CCL & CGDL

    % TEMPERATURE
    for d = 2:Nd
        res(2*(d-1)*Neq+5) = ya(5,d) - yb(5,d-1); % potential continuity
        res(2*(d-1)*Neq+6) = ya(6,d) - yb(6,d-1); % flux continuity
    end
    res(0*Neq+5) = ya(5,1) - T_A; % anode boundary temperature
    res(0*Neq+6) = yb(5,5) - T_C; % cathode boundary temperature

    % DISSOLVED WATER
    res(2*Neq+7) = ya(8,4) - yb(8,3); % flux continuity between PEM & CCL
    res(2*Neq+8) = ya(8,2); % zero flux between AGDL & ACL
    res(4*Neq+7) = ya(7,3) - yb(7,2); % potential continuity between ACL & PEM
    res(4*Neq+8) = ya(8,3) - yb(8,2); % flux continuity between ACL & PEM
    res(6*Neq+7) = ya(7,4) - yb(7,3); % potential continuity between PEM & CCL
    res(6*Neq+8) = yb(8,4); % zero flux between CCL & CGDL

    % WATER VAPOR
    res(0*Neq+ 9) = ya( 9,1) - w_H2O_A; % anode GC vapor content
    res(0*Neq+10) = ya(10,2) - yb(10,1); % flux continuity between AGDL & ACL
    res(2*Neq+ 9) = ya( 9,2) - yb( 9,1); % potential continuity between AGDL & ACL
    res(2*Neq+10) = yb(10,2); % zero flux between ACL & PEM
    res(6*Neq+ 9) = ya( 9,5) - yb( 9,4); % potential continuity between CCL & CGDL
    res(6*Neq+10) = ya(10,4); % zero flux between PEM & CCL
    res(8*Neq+ 9) = yb( 9,5) - w_H2O_C; % cathode GC vapor content
    res(8*Neq+10) = ya(10,5) - yb(10,4); % flux continuity between CCL & CGDL
    
    % OXYGEN
    res(6*Neq+11) = ya(11,5) - yb(11,4); % potential continuity between CCL & CGDL
    res(6*Neq+12) = ya(12,4); % zero flux between PEM & CCL
    res(8*Neq+11) = yb(11,5) - w_O2_C; % cathode GC oxygen content
    res(8*Neq+12) = ya(12,5) - yb(12,4); % flux continuity between CCL & CGDL

    % LIQUID WATER
    res(6*Neq+13) = ya(13,5) - yb(13,4); % potential continuity between CCL & CGDL
    res(6*Neq+14) = ya(14,4); % zero flux between PEM & CCL
    res(8*Neq+13) = yb(13,5) - s_C; % cathode GC liquid water content
    res(8*Neq+14) = ya(14,5) - yb(14,4); % flux continuity between CCL & CGDL
    
    % GAS
    res(0*Neq+15) = ya(15,1) - P_A; % anode GC pressure
    res(0*Neq+16) = ya(16,2) - yb(16,1); % flux continuity between AGDL & ACL
    res(2*Neq+15) = ya(15,2) - yb(15,1); % potential continuity between AGDL & ACL
    res(2*Neq+16) = yb(16,2); % zero flux between ACL & PEM
    res(6*Neq+15) = ya(15,5) - yb(15,4); % potential continuity between CCL & CGDL
    res(6*Neq+16) = ya(16,4); % zero flux between PEM & CCL
    res(8*Neq+15) = yb(15,5) - P_C; % cathode GC vapor content
    res(8*Neq+16) = ya(16,5) - yb(16,4); % flux continuity between CCL & CGDL
    
end

end