function [delta_d,JRC_mob,phi_d_mob] = ...
    FFSA_ShearDilation(delta_s,delta_peak,JRC,JCS,sigma_eff,phi_r,M,method)
% *************************************************************************
% FFSA_ShearDilation.m
% Created by Michael Liem on 21.06.2022
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
% delta_s           Shear displacement                                 [mm]
% delta_peak        Peak shear displacement                            [mm]
% JRC               Joint roughness coefficient                         [-]
% JCS               Joint wall compressive strength                   [MPa]
% sigma_eff         Effective normal stress                           [MPa]
% phi_r             Residual friction angle                             [°]
% M                 Damage coefficient                                  [-]
%
% method            Method how to calculate shear dilation         [string]
%                      'integrate'      Integrate over delta_s
%                      'integrate_pos'  Integrate only phi_d_mob>0
%                      'last'           Use last value (no integration)
%
% -------------------------------------------------------------------------
%
% delta_d           Shear dilation                                     [mm]
% JRC_mob           Mobilized joint roughness coefficient               [-]
% phi_d_mob         Mobilized dilation angle                            [°]
%
% *************************************************************************


% Check if effective normal stress is positive 
% -> This code does not consider tensile opening
if sigma_eff<0
    error(['Negative sigma_eff = ' num2str(sigma_eff) ' MPa. Tensile opening!'])
end

% Check if shear displacement is zero (or even negative?)
if delta_s<=0
    % Mobilized JRC
    JRC_mob= FFSA_JRCmob(0,sigma_eff,JRC,JCS,phi_r);

    % Mobilized dilation angle
    phi_d_mob= NaN;
    
    % Shear dilation
    delta_d= 0;
    
    return;
end


switch method
    case 'integrate'
        % Numerical integration width 
        d_delta_s= 0.1;

        % List of delta_s where JRC_mob need to be calculated
        if delta_s>= d_delta_s
            all_delta_s= d_delta_s/2:d_delta_s:delta_s;
        else
            all_delta_s= delta_s;
            d_delta_s= delta_s;
        end

        % Mobilized JRC
        all_JRC_mob= FFSA_JRCmob(all_delta_s/delta_peak,sigma_eff,JRC,JCS,phi_r);

        % Mobilized dilation angle
        all_phi_d_mob= 1/M*all_JRC_mob*log10(JCS/sigma_eff);

        % Incremental shear dilation
        d_delta_d= d_delta_s*tand(all_phi_d_mob);

        % Total shear dilation
        delta_d= sum(d_delta_d);

        % Final values
        JRC_mob= all_JRC_mob(end);
        phi_d_mob= all_phi_d_mob(end);
        
    case 'integrate_pos'
        % Numerical integration width 
        d_delta_s= 0.01;

        % List of delta_s where JRC_mob need to be calculated
        if delta_s>= d_delta_s
            all_delta_s= d_delta_s/2:d_delta_s:delta_s;
        else
            all_delta_s= delta_s;
            d_delta_s= delta_s;
        end

        % Mobilized JRC
        all_JRC_mob= FFSA_JRCmob(all_delta_s/delta_peak,sigma_eff,JRC,JCS,phi_r);

        % Mobilized dilation angle
        all_phi_d_mob= 1/M*all_JRC_mob*log10(JCS/sigma_eff);

        % Incremental shear dilation
        d_delta_d= d_delta_s*tand(all_phi_d_mob);

        % Total shear dilation
        delta_d= sum(d_delta_d(d_delta_d>=0));

        % Final values
        JRC_mob= all_JRC_mob(end);
        phi_d_mob= all_phi_d_mob(end);
        
    case 'last'
        % Mobilized JRC
        JRC_mob= FFSA_JRCmob(delta_s/delta_peak,sigma_eff,JRC,JCS,phi_r);

        % Mobilized dilation angle
        phi_d_mob= 1/M*JRC_mob*log10(JCS/sigma_eff);

        % Shear dilation
        delta_d= delta_s*tand(phi_d_mob);
    otherwise
        error(['FFSA_ShearDilation.m: Unknown method ' method])
end