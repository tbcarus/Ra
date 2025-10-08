function xyuvOut = xyuv(lambda_nm, spd)
% Вход:
%   lambda_nm : [n x 1] вектор длин волн (нм), твоя рабочая сетка
%   spd      : [n x 1] спектр на той же сетке lambda_nm
%
% Выход:
% Выход (struct xyuvOut):
%   .X     - значение тристимула X
%   .Y     - значение тристимула Y
%   .Z     - значение тристимула Z
%   .XYZ   - вектор [X;Y;Z]
%   .xy    - координаты [x; y] (CIE 1931)
%   .uv    - координаты [u'; v'] (CIE 1976)

[x,y,z] = data.cieXyz(lambda_nm);
x(isnan(x)) = 0;
y(isnan(y)) = 0;
z(isnan(z)) = 0;

X = trapz(lambda_nm, spd.*x);
Y = trapz(lambda_nm, spd.*y);
Z = trapz(lambda_nm, spd.*z);
XYZ = [X; Y; Z];

% нормированные хроматические координаты (x,y)
den_xy = max(X + Y + Z, eps);
x = X / den_xy;
y = Y / den_xy;
xy = [x; y];

% координаты u' v' (CIE 1976)
den_uv = max(X + 15*Y + 3*Z, eps);
up = 4*X / den_uv;
vp = 9*Y / den_uv;
uv = [up; vp];

xyuvOut.X   = X;
xyuvOut.Y   = Y;
xyuvOut.Z   = Z;
xyuvOut.XYZ = XYZ;
xyuvOut.xy  = xy;
xyuvOut.uv  = uv;

end

