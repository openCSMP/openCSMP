function [delta_s,JRC_mob,phi_s_mob] = FFSA_ShearDisplacement(...
    sigma_s,K_s,delta_peak,JRC,JCS,sigma_eff,sigma_EFF,phi_r,method,verbose)
% *************************************************************************
% FFSA_ShearDisplacement.m
% Created by Michael Liem on 22.11.2022
% Copyright (C) 2022 IFD, ETH Zurich.
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
% sigma_s           Shear stress                                      [MPa]
% K_s               Shear stiffness of fracture                    [MPa/mm]
% delta_peak        Peak shear displacement                            [mm]
% JRC               Joint roughness coefficient                         [-]
% JCS               Joint wall compressive strength                   [MPa]
% sigma_eff         Effective normal stress (for delta_tau or tau)    [MPa]
% sigma_EFF         Effective normal stress (for JRC_mob & phi_s_mob) [MPa]
% phi_r             Residual friction angle                             [°]
%
% method            Method how to calculate shear dilation         [string]
%                      'mob': friction angle is function of delta_s
%                      'old': friction angle = residual friction angle
%                             (this was used in ECMOR22 paper)
% verbose           Print some information                              [-]
%                      0: no info (default)
%                      1: command window output of convergence
%                      2: plot additionally a figure
%
% -------------------------------------------------------------------------
%
% delta_s           Shear displacement                                 [mm]
% JRC_mob           Mobilized joint roughness coefficient               [-]
% phi_s_mob         Mobilized friction angle                            [°]
%
% *************************************************************************


% Set to default value
if ~exist('verbose','var')
    verbose=0;
end


% Check if effective normal stress is positive 
% -> This code does not consider tensile opening
if sigma_eff<0
    error(['Negative sigma_eff = ' num2str(sigma_eff) ' MPa. Tensile opening!'])
end
if sigma_EFF<0
    error(['Negative sigma_EFF = ' num2str(sigma_EFF) ' MPa. Tensile opening!'])
end


% Legacy code: option to calculate delta_s as in ECMOR22 paper
if strcmp(method,'old')
    % Excess shear stress                                             [MPa]
    delta_tau= abs(sigma_s) - sigma_eff*tand(phi_r);
    
    % Shear displacement                                               [mm]
    delta_s= max(delta_tau/K_s,0);
    
    % Ratio of shear displacement over peak shear displacement          [-]
    d_by_dpeak= delta_s/delta_peak;
    
    % Mobilized JRC                                                     [-]
    JRC_mob= FFSA_JRCmob(d_by_dpeak,sigma_EFF,JRC,JCS,phi_r);
    
    % Mobilized friction angle                                          [°]
    phi_s_mob= JRC_mob * log10(JCS/sigma_EFF) + phi_r;
    
    return;
elseif ~strcmp(method,'mob')
    error(['Unknown method: ' method])
end


% Create an initial list of candidates for shear displacement          [mm]
delta_s_cand= [0:0.05:1.9 2:10 20:10:100]*delta_peak;

% Shear displacement for which relaxed shear stress is zero            [mm]
% -> This would be the solution without friction
delta_s_max= abs(sigma_s)/K_s;

% Limit candidates of delta_s so that only 1 entry is >= delta_s_max
delta_s_cand= delta_s_cand(1:min(max(sum(delta_s_cand<delta_s_max)+1,2),length(delta_s_cand)));


% Prepare figure (if desired)
if verbose>=2
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    figure('units','normalized','outerposition',[0 0 1 1])
end
if verbose >=1
    disp('FFSA_ShearDisplacement')
end

% Value for iteration
i=0;

% Max. numbers of iterations
max_iter= 10;

% Tolerance value (absolute value)                                    [MPa]
tolerance= 1e-3;

% Find delta_s such that relaxed shear stress = shear strengh
while true
    i= i+1;
    % Relaxed shear stress for each candidate                         [MPa]
    sigma_s_rel= max( abs(sigma_s) - K_s*delta_s_cand , 0);
    
    % Ratio of shear displacement over peak shear displacement          [-]
    d_by_dpeak= delta_s_cand/delta_peak;
    % Mobilized JRC                                                     [-]
    JRC_mob= FFSA_JRCmob(d_by_dpeak,sigma_EFF,JRC,JCS,phi_r);
    % Mobilized friction angle                                          [°]
    phi_s_mob= JRC_mob * log10(JCS/sigma_EFF) + phi_r;
    % Shear strength for each candidate of delta_s
    tau= sigma_eff*tand(phi_s_mob);
    
    % Absolute difference between relaxed shear stress and shear strength
    residual= abs(sigma_s_rel-tau);
    
    % Plot some figure (if desired)
    if verbose>=2
        subplot(2,ceil(max_iter/2),i)
        plot(delta_s_cand,tau)
        hold on
        plot(delta_s_cand,sigma_s_rel)
        plot(delta_s_cand,residual)
        legend({'$\tau$','$\sigma_{s}^{rel}$','$|\tau-\sigma_{s}^{rel}|$'})
        title(['iter ' num2str(i)])
        xlabel('$\delta_s$ [mm]')
        ylabel('MPa')
        set(gca,'FontSize',16)
    end
    
    % Find index where residual is minimal
    idx= find(min(residual)==residual,1,'first');
    
    % Print information about convergence (if desired)
    if verbose >=1
        disp(['   iter = ' num2str(i) ': delta_s = ' num2str(delta_s_cand(idx)) ' mm'...
            ' , residual = ' num2str(residual(idx)) ' MPa'])
    end
    
    % Break loop if residual is small enough or after 10 iterations
    if residual(idx)<=tolerance || i>max_iter
        break;
    end
    
    % Choose one index lower and higher as bounds for next iteration
    % -> Make sure that index are valid
    idx1= max(idx-1,1);
    idx2= min(idx+1,length(residual));
    
    % Create new candidate list
    delta_s_cand= [0 .25 .4 .45 .48 .5 .52 .55 .6 .75 1]*...
        (delta_s_cand(idx2)-delta_s_cand(idx1)) + delta_s_cand(idx1);
    
    % Limit candidates of delta_s so that only 1 entry is >= delta_s_max
    delta_s_cand= ...
        delta_s_cand(1:min(max(sum(delta_s_cand<delta_s_max)+1,2),length(delta_s_cand)));
end


% Shear displacement                                                   [mm]
delta_s= delta_s_cand(idx);

% Mobilized joint roughness coefficient                                 [-]
JRC_mob= JRC_mob(idx);

% Mobilized friction angle                                              [°]
phi_s_mob= phi_s_mob(idx);
