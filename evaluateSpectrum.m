function E = evaluateSpectrum(L_interp, spd, target)
% EVALUATESPECTRUM  Базовая оценка спектра (xy, u'v', XYZ, + сравнение с target)
%
% Вход:
%   L_interp : [n x 1] сетка длин волн (нм)
%   spd      : [n x 1] спектр смеси (например, fit)
%   target : [n x 1] эталонный спектр на той же сетке (опц.)
%
% Выход (struct E):
%   E.XYZ  : [3x1] тристимульные значения
%   E.xy   : [2x1] координаты x,y
%   E.uv   : [2x1] координаты u',v'
%   E.dxy  : скаляр, расстояние между xy(fit) и xy(target) (если задан target)
%   E.duv  : скаляр, расстояние между u'v'(fit) и u'v'(target) (если задан target)
%
% Позже сюда добавим: E.CCT, E.Ra, E.BHL и т.д.


% базовые метрики по текущему спектру xy/uv
    xyuvOut = ra.xyuv(L_interp, spd);
    E = xyuvOut;

    E.dxy = NaN;
    E.duv = NaN;
    if ~isempty(target)
        xyuvOut_t = ra.xyuv(L_interp, target);
        E.dxy = norm(xyuvOut.xy - xyuvOut_t.xy, 2);
        E.duv = norm(xyuvOut.uv - xyuvOut_t.uv, 2);
    end


end

