%% Led Flux and U, I (for binning conditions)
function [Params, T, F_Sun_AM15] = led_config()

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
setv(ra.Channels.PC_AMBER,  310,  1, 2.95);
setv(ra.Channels.RED,  83,  0.35, 2.08);
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

% F_2700 = 350;
% I_2700 = 1.05;
% U_2700 = 2.79;
% P_2700 = U_2700 * I_2700;
%
% F_6500 = 430;
% I_6500 = 1.05;
% U_6500 = 2.79;
% P_6500 = U_6500 * I_6500;
%
% F_RB = 0.712;
% I_RB = 0.35;
% U_RB = 2.93;
% P_RB = U_RB * I_RB;
%
% F_B = 55;
% I_B = 0.35;
% U_B = 2.85;
% P_B = U_B * I_B;
%
% F_Cyan = 262;
% I_Cyan = 0.7;
% U_Cyan = 2.71;
% P_Cyan = U_Cyan * I_Cyan;
%
% F_PCCyan_XEG = 330;
% I_PCCyan_XEG = 1;
% U_PCCyan_XEG = 2.9;
% P_PCCyan_XEG = U_PCCyan_XEG * I_PCCyan_XEG;
%
% F_PCCyan_XQE = 110;
% I_PCCyan_XQE = 0.35;
% U_PCCyan_XQE = 2.8;
% P_PCCyan_XQE = U_PCCyan_XQE * I_PCCyan_XQE;
%
% F_Lime = 500;
% I_Lime = 1;
% U_Lime = 2.95;
% P_Lime = U_Lime * I_Lime;
%
% F_G = 184;
% I_G = 0.35;
% U_G = 2.7;
% P_G = U_G * I_G;
%
% F_A = 310;
% I_A = 1;
% U_A = 2.95;
% P_A = U_A * I_A;
%
% F_R = 83;
% I_R = 0.35;
% U_R = 2.08;
% P_R = U_R * I_R;
%
% F_DR = 0.512;
% I_DR = 0.35;
% U_DR = 2.09;
% P_DR = U_DR * I_DR;
%
% F_FR = 0.437;
% I_FR = 0.35;
% U_FR = 1.88;
% P_FR = U_FR * I_FR;

F_Sun_AM15 = 9.1680e+04;

end