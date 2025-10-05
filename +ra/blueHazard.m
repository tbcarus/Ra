
function out = blueHazard(L_interp, spd, beam, Dm, dist)
% Индекс "blue-light hazard"
% beam - угол излучения
% Dm - диаметр излучателя
% dist - расстояние до излучателя

B = ra.blh(L_interp);

% Геометрия
A_em   = pi*(Dm/2)^2;                          % площадь свечения
Omega  = 2*pi*(1-cos(deg2rad(beam/2)));        % телесный угол
alpha  = Dm / dist;                            % угловой размер источника [rad]
isExt  = alpha > 0.1;                          % >100 mrad, то протяжённый источник

% Интенсивность и радиантность (топ-хэт):
% I_λ = Φ_λ / Ω  [W/sr/nm]
% L_λ = I_λ / A_proj ≈ I_λ / A_em   [W/m^2/sr/nm]
I_lambda = spd / Omega;
L_lambda = I_lambda / A_em;

% Освещённость на оси в точке r: E_λ = I_λ / r^2  [W/m^2/nm]
E_lambda = I_lambda / dist^2;

% Интегралы по λ
LB = trapz(L_interp, L_lambda .* B);          % [W/m^2/sr]
EB = trapz(L_interp, E_lambda .* B);          % [W/m^2]
BH_index = trapz(L_interp, spd .* B) / trapz(L_interp, spd);   % относительный индекс

out.LB_W_m2_sr   = LB;
out.EB_W_m2      = EB;
out.BH_index     = BH_index;
out.Omega_sr     = Omega;
out.A_em_m2      = A_em;
out.alpha_src_rad= alpha;
out.isExtended   = isExt;

end


