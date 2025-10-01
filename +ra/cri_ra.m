function out = cri_ra(lambda_nm, spd)
% Ra по CIE 13.3-1995 (упрощённо).
% Требуются спектры отражения TCS_k(λ). Вставь их в tcs_reflectance().
%
% Шаги:
% 1) Выбрать опорный источник: ниже 5000K -> Планк; выше -> D-источник
% 2) Получить XYZ для каждого TCS под тестовым SPD и опорным
% 3) Перевести в CIE 1964 U*V*W или CIE 1976 L*a*b* (в стандарте V(λ) масштабы)
% 4) ΔE_k -> R_k = 100 - 4.6*ΔE_k; Ra = mean(R_1..R_8)

% 0) u'v' и CCT для определения опорного источника
cct = ra.cct_exact(lambda_nm, spd).CCT;

% 1) Опорный SPD
if cct < 5000
    spd_ref = planck_norm(lambda_nm, cct);
else
    spd_ref = daylight_D(lambda_nm, cct);  % каркас; реализация ниже
end

% 2) Получить TCS
TCS = tcs_reflectance(lambda_nm);  % struct с полями R{k} [n x 1], k=1..8

% 3) Для каждого TCS вычислить координаты в выбранном цветовом пространстве
Rk = zeros(8,1);
for k = 1:8
    [L1,a1,b1] = lab_from_spd(lambda_nm, spd,     TCS{k});
    [L2,a2,b2] = lab_from_spd(lambda_nm, spd_ref, TCS{k});
    dE = sqrt((L1-L2)^2 + (a1-a2)^2 + (b1-b2)^2);
    Rk(k) = max(0, 100 - 4.6*dE);
end

out.Ra = mean(Rk);
out.Rk = Rk;
end

% ===== place-holders / helpers =====

function spd = planck_norm(lambda,T)
spd = (ra.cct_exact(lambda, ones(size(lambda))).CCT); %#ok<NASGU>
% просто используем локальную имплементацию:
c2 = 1.438776877e7; lam=lambda(:);
spd = lam.^(-5) ./ (exp(c2./(lam*T)) - 1);
spd = spd / trapz(lambda, spd);
end

function spd = daylight_D(lambda, CCT)
% Заглушка: позже вставим стандарт D-источников (CIE 15)
% Пока — приближение: возьмём нормированный Планк, чтобы интерфейс работал.
spd = planck_norm(lambda, CCT);
end

function TCS = tcs_reflectance(lambda)
% ЗАГЛУШКА: сюда вставь стандартные отражательные спектры TCS1..TCS14,
% интерполированные на lambda. Пока вернём пустые/константные (НЕ ДЛЯ РАБОТЫ!)
TCS = cell(8,1);
for k=1:8
    TCS{k} = ones(numel(lambda),1)*0.5; % <- замени на реальные данные!
end
end

function [L,a,b] = lab_from_spd(lambda, spd_src, R)
% SPD объекта = SPD источника .* R(λ)
spd_obj = spd_src(:) .* R(:);
% XYZ
m = ra.xyuv(lambda, spd_obj);
XYZ = m.XYZ;  % относительные
% Белая точка = опорный источник
white = ra.xyuv(lambda, spd_src).XYZ;
% XYZ -> CIE L*a*b* (через белую точку)
[L,a,b] = xyz2lab(XYZ, white);
end

function [L,a,b] = xyz2lab(XYZ, white)
% Быстрая XYZ->Lab. Предполагаем нормировку по белой точке.
f = @(t) (t>(6/29)^3).*(t.^(1/3)) + (t<=(6/29)^3).*(t/(3*(6/29)^2)+4/29);
X = XYZ(1)/max(white(1),eps);
Y = XYZ(2)/max(white(2),eps);
Z = XYZ(3)/max(white(3),eps);
fx=f(X); fy=f(Y); fz=f(Z);
L = 116*fy - 16;
a = 500*(fx - fy);
b = 200*(fy - fz);
end


