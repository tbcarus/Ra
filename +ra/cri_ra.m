function CRI = cri_ra(L_interp, spd)
% CRI_Ra  Считает общий индекс цветопередачи Ra и отдельные R1..R14.
% Вход:
%   L_interp : [n×1] длины волн (нм)
%   spd      : [n×1] спектр тестируемого источника (Relative SPD)
% Выход (struct CRI):
%   .Ra      — общий индекс (среднее R1..R8)
%   .Rk      — [14×1] индексы по каждому TCS (R1..R14)
%   .CCT     — оценённая Тц тестового источника (для информации)

    % 0) CCT тестового источника (для выбора опорного)
    C = ra.cctExact(L_interp, spd);
    CCT = C.CCT;

    % 1) Опорный источник
    if CCT < 5000
        spd_ref = ra.planckSpd(L_interp, CCT);
    else
        spd_ref = daylight_D(L_interp, CCT);
    end

    % === Белая поверхность (R=1) под обоими источниками ===
    ones_reflector = ones(numel(L_interp),1);
    Y_test_white = ra.xyuv(L_interp, spd .* ones_reflector).Y;
    Y_ref_white  = ra.xyuv(L_interp, spd_ref .* ones_reflector).Y;

    % === Подгоняем Y опорного под тестовый (классическая практика CRI) ===
    scale = Y_test_white / max(Y_ref_white, eps);
    spd_ref = spd_ref * scale;



    % 2) Отражательные спектры образцов R
    R = data.rObjects(L_interp);   % [n×14]



    % 3) Индексы для R1..R14
    Rk = zeros(14,1);
    for k = 1:14
        % Под данным источником
        spd_obj_test = spd.* R(:,k);
        % Под опорным источником
        spd_obj_ref  = spd_ref.* R(:,k);

        % Цветовые координаты (XYZ, xy, uv) → далее в Lab (с белой точкой = опорный источник)
        col_test = ra.xyuv(L_interp, spd_obj_test);
        col_ref  = ra.xyuv(L_interp, spd_obj_ref);

        % Белая точка = опорный источник (XYZ белой поверхности под опорным источником)
        white_ref = ra.xyuv(L_interp, spd_ref).XYZ;

        % XYZ -> Lab (упрощение; в строгом CRI используется своя шкала/методика,
        % но LAB даёт близкий результат – при желании заменим на U*V*W)
        [L1,a1,b1] = xyz2lab(col_test.XYZ, white_ref);
        [L2,a2,b2] = xyz2lab(col_ref.XYZ,  white_ref);

        dE = sqrt((L1-L2)^2 + (a1-a2)^2 + (b1-b2)^2);

        % Индекс для образца k
        Rk(k) = max(0, 100 - 4.6*dE);
    end

    % 4) Итоги
    CRI.Rk  = Rk;
    CRI.Ra  = mean(Rk(1:8));
    CRI.CCT = CCT;

    % === 5) Диагностика на случай "всё нули" ===
    if ~all(isfinite(Rk)) || all(Rk==0)
        warning('CRI: R_k contains non-finite or all zeros. Check: daylight_D, rObjects, XYZ white point, Y-normalization.');
        debug.CR_Y_test_white = Y_test_white;
        debug.CR_Y_ref_white  = Y_ref_white;
        debug.Rk              = Rk;
        CRI.debug = debug;
    end

end


function spd = daylight_D(L_interp, CCT)
% DAYLIGHT_D  CIE Daylight SPD при заданной CCT (K).
% Формула CIE: S(λ) = S0(λ) + M1*S1(λ) + M2*S2(λ),
% где M1, M2 определяются через x_D(CCT), y_D(CCT).
%
% Вход:
%   L_interp : [n×1] сетка длин волн (нм), на которой нужен SPD
%   CCT      : скаляр, требуемая коррелированная цветовая температура (K)
%
% Выход:
%   spd      : [n×1] спектральная плотность D-источника (отн.)
[S0, S1, S2] = ra.S0_S1_S2(L_interp);
    % Цветность дневного освещения x_D, y_D как функции CCT ===
    xD = xD_from_CCT(CCT);
    yD = -3.000*xD.^2 + 2.870*xD - 0.275;   % полином CIE для y_D(x_D)

    % Коэффициенты M1, M2 (зависят от xD, yD) ===
    denom = (0.0241 + 0.2562*xD - 0.7341*yD);
    M1 = (-1.3515 - 1.7703*xD + 5.9114*yD) / denom;
    M2 = ( 0.0300 -31.4424*xD +30.0717*yD) / denom;

    % === 4) Финальный SPD = S0 + M1*S1 + M2*S2 ===
    spd = S0 + M1*S1 + M2*S2;
end

function xD = xD_from_CCT(T)
% Цветовая координата x_D(CCT) CIE Daylight
% Два куска: 4000..7000 K и 7000..25000 K
    if T < 4000,  T = 4000;  end
    if T > 25000, T = 25000; end

    if T >= 4000 && T <= 7000
        xD = 0.244063 + 0.09911e3/T + 2.9678e6/T^2 - 4.6070e9/T^3;
    else % 7000..25000
        xD = 0.237040 + 0.24748e3/T + 1.9018e6/T^2 - 2.0064e9/T^3;
    end
end

function [L,a,b] = xyz2lab(XYZ, whiteXYZ)
% Быстрая XYZ->CIELAB с белой точкой whiteXYZ
    fr = @(t) (t>(6/29)^3).*t.^(1/3) + (t<=(6/29)^3).*(t/(3*(6/29)^2)+4/29);
    X = XYZ(1)/max(whiteXYZ(1),eps);
    Y = XYZ(2)/max(whiteXYZ(2),eps);
    Z = XYZ(3)/max(whiteXYZ(3),eps);
    fX=fr(X); fY=fr(Y); fZ=fr(Z);
    L = 116*fY - 16;
    a = 500*(fX - fY);
    b = 200*(fY - fZ);
end


