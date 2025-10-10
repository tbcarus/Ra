%% Led Flux and U, I (for binning conditions)
function [Params, T, F_Sun_AM15] = ledData()

% Возвращает:
%   Params [m x 4] в порядке колонок LedParam (F, I, U, P=U.*I)
%   T      table с колонками Name, F, I, U, P
%   names, name2idx — из ra.led_order (для удобства)

[names, name2idx] = ra.led_order();
m = numel(names);
% Инициализация NaN (чтобы было сразу видно, что не заполнено)
F = nan(m,1); I = nan(m,1); U = nan(m,1);

setv(ra.Channels.WARM2700,  350,  1.05, 2.79);
setv(ra.Channels.COOL6500,  430,  1.05, 2.79);
setv(ra.Channels.ROYAL_BLUE,  0.71,  0.35, 2.93);
setv(ra.Channels.BLUE,  55,  0.35, 2.85);
setv(ra.Channels.CYAN,  262,  0.7, 2.71);
setv(ra.Channels.PC_CYAN_XEG,  330,  1, 2.9);
setv(ra.Channels.PC_CYAN_XQE,  110,  0.35, 2.8);
setv(ra.Channels.LIME,  500,  1, 2.95);
setv(ra.Channels.GREEN,  184,  0.35, 2.7);
setv(ra.Channels.AMBER,  103,  0.32, 2.18);
setv(ra.Channels.PC_AMBER,  310,  1, 2.95);
setv(ra.Channels.RED,  55,  1, 3);
setv(ra.Channels.PC_RED,  83,  0.35, 2.08);
setv(ra.Channels.DEEP_RED,  0.512,  0.35, 2.09);
setv(ra.Channels.FAR_RED,  0.437,  0.35, 1.88);

P = U .* I;
Params = [F, I, U, P];

T = table(string(names(:)), F, I, U, P, ...
    'VariableNames', {'Name','F','I','U','P'});

% Валидация
missing = any(isnan(Params),2);
if any(missing)
    warnNames = strjoin(string(names(missing)), ', ');
    warning('LED config: не все параметры заполнены: %s', warnNames);
end

    function setv(chName, Fv, Iv, Uv)
        idx = name2idx(chName);
        F(idx) = Fv;
        I(idx) = Iv;
        U(idx) = Uv;
    end

F_Sun_AM15 = 9.1680e+04;

end