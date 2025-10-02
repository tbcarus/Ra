function cctOut = cctExact(L_interp, spd)
% CCT_EXACT  Точная коррелированная цветовая температура (Тц, CCT)
% методом минимизации расстояния до планковского локуса в u'v'.
%
% Вход:
%   L_interp : [n x 1] длины волн (нм) — твоя рабочая сетка
%   spd      : [n x 1] спектр источника (Relative SPD) на L_interp
%
% Выход (struct outCCT):
%   .CCT        - оптимальная температура (K), дающая минимальную Δu'v'
%   .duv        - минимальное расстояние до локуса в плоскости u'v' (без знака)
%   .uv_src     - u'v' твоего спектра
%   .uv_bb      - u'v' чёрного тела при CCT
%   .T_bounds   - границы поиска (K)
%   .method     - метка метода ('fminbnd on Planck locus')

% 1) u'v' для оцениваемого спектра
src = ra.xyuv(L_interp, spd);
uv_src = src.uv;

% 2) Определяем функцию: T -> расстояние Δu'v' до черного тела при T
fun = @(T) duv_to_bb(L_interp, uv_src, T);

% 3) Поиск минимума по температуре (границы можно подстроить)
T_lo = 2000;   % K
T_hi = 10000;  % K
opts = optimset('TolX',1,'Display','off');
[T_opt, duv_min] = fminbnd(fun, T_lo, T_hi, opts);

% 4) u'v' для черного тела при найденной температуре
uv_bb = ra.xyuv(L_interp, planck_spd(L_interp, T_opt)).uv;

% 5) Результат
cctOut.CCT      = T_opt; % коррелированная Тц
cctOut.duv      = duv_min; % минимальное расстояние между спектром чёрным телом в пространстве u′v′.
cctOut.uv_src   = uv_src; % координаты u′v′ источника
cctOut.uv_bb    = uv_bb; % % координаты u′v′ АЧТ
cctOut.T_bounds = [T_lo, T_hi]; % границы интервала поиска Тц
cctOut.method   = 'fminbnd on Planck locus';


function d = duv_to_bb(L_interp, uv_src, T)
    spd_bb = planck_spd(L_interp, T);          % спектр чёрного тела при T
    uv_bb  = ra.xyuv(L_interp, spd_bb).uv;        % его u'v'
    d      = hypot(uv_src(1)-uv_bb(1), uv_src(2)-uv_bb(2));  % расстояние
end

function spd = planck_spd(L_interp, T)
% Спектр чёрного тела B(λ, T) в относительных единицах.
% Формула Планка : c1/λ^5 * 1/(exp(c2/(λ*T)) - 1)
% где λ — в НАНОМЕТРАХ, T — в Кельвинах, c2 ≈ 1.44e7 нм·K.

h = 6.63e-34; % приведённая постоянная Планка Дж*с
k = 1.38e-23; % постоянная Больцмана Дж/К
c = 3e8*1e9; % скорость света нм/с
    c1 = 2*pi*h*c^2; % первая константа излучения
    c2  = h*c/k;       % вторая константа излучения
    spd = c1*L_interp.^(-5) ./ (exp(c2./(L_interp*T)) - 1);
    spd = spd / max(spd);      % нормируем по пику, не влияет на u'v'
end

end

