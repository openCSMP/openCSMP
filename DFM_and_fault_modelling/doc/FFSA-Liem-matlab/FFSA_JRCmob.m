function [JRC_mob] = FFSA_JRCmob(d_by_dpeak,sigma_n,JRC,JCS,phi_r)
% *************************************************************************
% FFSA_JRCmob.m
% Created by Michael Liem on 28.04.2022
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
% d_by_dpeak        Ratio shear displacement to peak shear displacement [-]
%                   -> Has to be in range [0 100]
% sigma_n           Normal stress (or effective normal stress)        [MPa]
% JRC               (Peak) Joint roughness coefficient                  [-]
% JCS               Joint wall compressive strength                   [MPa]
% phi_r             Residual friction angle                             [°]
%
% -------------------------------------------------------------------------
%
% JRC_mob           Mobilized joint roughness coefficient               [-]
%
% *************************************************************************


i= JRC*log10(JCS/sigma_n); % [°]

list_d_by_dpeak=        [0         0.3  0.6   1  2     4    10    100];
list_JRCmob_by_JRCpeak= [-phi_r/i  0    0.75  1  0.85  0.7   0.5    0];

% Check if any input values for shear displacement are negative
if any(d_by_dpeak<0)
    disp('FFSA_JRCmob.m: Warning: d_by_dpeak < 0 -> Taking absolute value!')
    d_by_dpeak= abs(d_by_dpeak);
end

% If ratio >= 100 -> JRC_mob = 0 anyway
d_by_dpeak(d_by_dpeak>100)= 100;

% Calculate JRC_mob according to e.g. Barton et al. (1985)
JRCmob_by_JRCpeak= interp1(list_d_by_dpeak, list_JRCmob_by_JRCpeak, d_by_dpeak);

JRC_mob= JRCmob_by_JRCpeak*JRC;


end

