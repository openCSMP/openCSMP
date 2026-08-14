function [result,debug] = FFSA_Aperture(params)
% *************************************************************************
% FFSA_Aperture.m
% Created by Michael Liem on 18.01.2023
% Copyright (C) 2023 IFD, ETH Zurich.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% *************************************************************************
%
% params            Struct with information about uncertainty of parameters
%   .L                Fracture length                                   [m]
%   .alpha            Fracture angle (i.e. between fracture and x-axis) [°]
%   .sigma_H          Max. principal far field stress                 [MPa]
%   .sigma_h          Min. principal far field stress                 [MPa]
%   .beta             Orientation of sigma_H                            [°]
%   .p_f              Fluid pressure                                  [MPa]
%   .JRC              Joint roughness coefficient                       [-]
%   .sigma_c          Unconfined compression strength                 [MPa]
%   .JCS              Joint wall compressive strength                   [-]
%   .K_ni             Initial normal stiffness                     [MPa/mm]
%   .vm_factor        Factor for maximum possible joint closure         [-]
%   .phi_r            Residual friction angle                           [°]
%   .E_mod            E-modulus                                       [MPa]
%   .nu               Poisson's ratio                                   [-]
%   .C_g              Proportionality between displacement and shear stress
%   .M                Damage coefficient (optional)                     [-]
%   .sigma_EFF        sigma_eff for phi_s_mob & phi_d_mob (optional)  [MPa]
%   .method           Method to calculate displacement & dilation  [string]
%
% -------------------------------------------------------------------------
%
% result            Array with resulting apertures                     [mm]
%                     [a_0 delta_n delta_s delta_d a]
%
% debug             Struct with all the variables
%
% *************************************************************************


fn= fieldnames(params);

% Find number of fractures
N_list= NaN(numel(fn),1);
for i=1:numel(fn)
    if isnumeric(params.(fn{i}))
        N_list(i)= length(params.(fn{i}));
    end
end
N_list= N_list(~isnan(N_list));
N_frac= max(N_list);

if any(N_list~=1 & N_list~=N_frac)
    error('Inconsistent number of fractures')
end

% Make vectors out of all parameters
for i=1:numel(fn)
    if isnumeric(params.(fn{i}))
        if any(isnan(params.(fn{i}))) && ~strcmp(fn{i},'M')
            error(['Input parameter ' fn{i} ' contains NaN values'])
        end
        if length(params.(fn{i})) ~= N_frac
            params.(fn{i})= params.(fn{i})*ones(N_frac,1);
        end
    end
end


result= NaN(N_frac,5);

debug.sigma_H=    NaN(N_frac,1);
debug.sigma_h=    NaN(N_frac,1);
debug.beta=       NaN(N_frac,1);
debug.p_f=        NaN(N_frac,1);
debug.sigma_n=    NaN(N_frac,1);
debug.sigma_s=    NaN(N_frac,1);
debug.sigma_eff=  NaN(N_frac,1);
debug.sigma_EFF=  NaN(N_frac,1);

debug.JRC=        NaN(N_frac,1);
debug.JCS=        NaN(N_frac,1);
debug.sigma_c=    NaN(N_frac,1);
debug.a_0=        NaN(N_frac,1);
debug.vm_factor=  NaN(N_frac,1);
debug.v_m=        NaN(N_frac,1);
debug.K_ni=       NaN(N_frac,1);
debug.delta_n=    NaN(N_frac,1);

debug.phi_r=      NaN(N_frac,1);
debug.E_mod=      NaN(N_frac,1);
debug.nu=         NaN(N_frac,1);
debug.C_g=        NaN(N_frac,1);
debug.G=          NaN(N_frac,1);
debug.K_s=        NaN(N_frac,1);
debug.delta_s=    NaN(N_frac,1);
debug.phi_s_mob=  NaN(N_frac,1);

debug.delta_peak= NaN(N_frac,1);
debug.JRC_mob=    NaN(N_frac,1);
debug.M=          NaN(N_frac,1);
debug.phi_d_mob=  NaN(N_frac,1);
debug.delta_d=    NaN(N_frac,1);
debug.a=          NaN(N_frac,1);
debug.k=          NaN(N_frac,1);


method_disp=      params.method.displacement;
method_dil=       params.method.dilation;
method_peak=      params.method.peak_displacement;
method_sigma_eff= params.method.sigma_eff;


for frac=1:N_frac
    % Extract parameters for current fracture
    L=          params.L(        frac);
    alpha=      params.alpha(    frac);
    sigma_H=    params.sigma_H(  frac);
    sigma_h=    params.sigma_h(  frac);
    beta=       params.beta(     frac);
    p_f=        params.p_f(      frac);
    JRC=        params.JRC(      frac);
    JCS=        params.JCS(      frac);
    sigma_c=    params.sigma_c(  frac);
    K_ni=       params.K_ni(     frac);
    vm_factor=  params.vm_factor(frac);
    phi_r=      params.phi_r(    frac);
    E_mod=      params.E_mod(    frac);
    nu=         params.nu(       frac);
    C_g=        params.C_g(      frac);
    
    % Angle between sigma_H and fracture normal                         [°]
    theta= 90+alpha-beta;
    
    
    % ---------------------------------------------------------------------
    % Part 1: Calculate effective stress
    sigma_n= sigma_H*cosd(theta)^2 + sigma_h*sind(theta)^2;
    sigma_s= (sigma_h-sigma_H)*cosd(theta)*sind(theta);
    sigma_eff= sigma_n - p_f;
    
    if sigma_eff<0
        error(['Negative sigma_eff = ' num2str(sigma_eff) ' MPa. Tensile opening!'])
    end
    
    % ---------------------------------------------------------------------
    % Part 2: Normal closure
    
    % Initial fracture aperture                                        [mm]
    a_0= JRC/5*(0.2*sigma_c/JCS - 0.1);
    
    % Maximum closure                                                  [mm]
    v_m= vm_factor*a_0;
        
    % Closure due to normal stress                                     [mm]
    delta_n= sigma_eff*v_m/( K_ni*v_m + sigma_eff );
    
    % ---------------------------------------------------------------------
    % Part 3: Shear displacement
    
    % Shear modulus                                                   [MPa]
    G= E_mod/(2*(1+nu));
    
    % Shear stiffness                                              [MPa/mm]
    K_s= C_g*G/(L*1e3);
    
    
    % Use a different effective stress to calculate                   [MPa]
    %    delta_peak, M, JRC_mob, phi_s_mob, phi_d_mob
    switch method_sigma_eff
        case 'end'
            % Take the value at the highest fluid pressure (default)
            sigma_EFF= sigma_eff;
        case 'onset'
            % Take the value approx. at onset of shear displacement
            sigma_EFF= sigma_n - sigma_s/tand(phi_r);
        case 'mean'
            % Take the mean of 'end' and 'onset'
            sigma_EFF= sigma_n - 0.5*(sigma_s/tand(phi_r) + p_f);
        case 'value'
            % Take the value from input file
            if ~isfield(params,'sigma_EFF')
                error("method.sigma_eff= 'value' but no input value for sigma_EFF")
            end
            sigma_EFF= params.sigma_EFF(frac);
        otherwise
            % Take the value at the highest fluid pressure (default)
            sigma_EFF= sigma_eff;
    end
    
    % Peak shear displacement                                          [mm]
    switch method_peak
        case 'Barton'
            delta_peak= L/500 * (JRC/L)^0.33 * 1e3;
        case 'Asadollahi'
            delta_peak= 0.0077*L^0.45*(sigma_EFF/JCS)^0.34*cosd( JRC*log10(JCS/sigma_EFF) )*1e3;
        otherwise
            error(['Unknown method for calculating peak shear displacment: ' method.peak_disp])
    end
    
    % Shear displacement                                               [mm]
    % Mobilized joint roughness coefficient                             [-]
    % Mobilized friciton angle                                          [°]
    [delta_s,JRC_mob,phi_s_mob] = FFSA_ShearDisplacement(...
        sigma_s,K_s,delta_peak,JRC,JCS,sigma_eff,sigma_EFF,phi_r,method_disp,0);
    
    % ---------------------------------------------------------------------
    % Part 4: Shear dilation

    % Damage coefficient                                                [-]
    if isfield(params,'M') && ~isnan(params.M(frac))
        M= params.M(frac);
    else
        M= 0.7 + JRC/(12*log10(JCS/sigma_EFF));
    end
    
    % Shear dilation                                                   [mm]
    % Mobilized dilation angle                                          [°]
    [delta_d,~,phi_d_mob] = FFSA_ShearDilation(...
        delta_s,delta_peak,JRC,JCS,sigma_EFF,phi_r,M,method_dil);
    
    % ---------------------------------------------------------------------
    % Part 5: Aperture
    a= a_0 - delta_n + delta_d;
    
    k= (a/1e3)^2 / 12;
    
    
    
    % ---------------------------------------------------------------------
    % Collecting results
    result(frac,1)= a_0;
    result(frac,2)= delta_n;
    result(frac,3)= delta_s;
    result(frac,4)= delta_d;
    result(frac,5)= a;
    
    debug.sigma_H(   frac)= sigma_H;
    debug.sigma_h(   frac)= sigma_h;
    debug.beta(      frac)= beta;
    debug.p_f(       frac)= p_f;
    debug.sigma_n(   frac)= sigma_n;
    debug.sigma_s(   frac)= sigma_s;
    debug.sigma_eff( frac)= sigma_eff;
    debug.sigma_EFF( frac)= sigma_EFF;
    
    debug.JRC(       frac)= JRC;
    debug.JCS(       frac)= JCS;
    debug.sigma_c(   frac)= sigma_c;
    debug.a_0(       frac)= a_0;
    debug.vm_factor( frac)= vm_factor;
    debug.v_m(       frac)= v_m;
    debug.K_ni(      frac)= K_ni;
    debug.delta_n(   frac)= delta_n;
    
    debug.phi_r(     frac)= phi_r;
    debug.E_mod(     frac)= E_mod;
    debug.nu(        frac)= nu;
    debug.C_g(       frac)= C_g;
    debug.G(         frac)= G;
    debug.K_s(       frac)= K_s;
    debug.delta_s(   frac)= delta_s;
    debug.phi_s_mob( frac)= phi_s_mob;
    
    debug.delta_peak(frac)= delta_peak;
    debug.JRC_mob(   frac)= JRC_mob;
    debug.M(         frac)= M;
    debug.phi_d_mob( frac)= phi_d_mob;
    debug.delta_d(   frac)= delta_d;
    debug.a(         frac)= a;
    debug.k(         frac)= k;
      
end

