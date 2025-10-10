clearvars;
close all;

%%
% Это основной исполняемый файл

%% === ДАННЫЕ ===
% names      : {1 x m} имена каналов в порядке столбцов S
% S          : [n x m] матрица спектров LED (каждый столбец нормирован к реальному потоку)
% lambda_led : [n x 1] вектор длин волн (нм), общая сетка для LED
% lambda_t   : [nt x 1] сетка эталона
% t_raw      : [nt x 1] эталонный спектр (нормирован к реальному потоку)

[Params, T, F_Sun_AM15] = data.ledData();
[names, name2idx] = ra.led_order();
[S, lambda_led, lambda_t, t_raw, V] = data.spectra(names, name2idx, Params, F_Sun_AM15);
assert(exist('lambda_led','var')==1 && exist('S','var')==1, 'Нет S или lambda_led');
assert(exist('lambda_t','var')==1 && exist('t_raw','var')==1, 'Нет эталона');
[n, m] = size(S);

%% === Интерполяция эталона на сетку LED ===
t = interp1(lambda_t(:), t_raw(:), lambda_led(:), 'linear', 0);  % вне диапазона -> 0
if max(t) <= 0
    error('После интерполяции максимум эталона <= 0. Проверьте данные.');
end

%% === (Опц.) Весовая функция по длинам волн ===
u = ones(n,1);                 % сюда можно подставить V(?) или чувствительность сенсора
W = spdiags(u, 0, n, n);

%% === (Опц.) Ограничения на веса диодов ===
disabled_names = {};
% disabled_names = {ra.Channels.WARM2700, ra.Channels.COOL6500};
% disabled_names = {ra.Channels.COOL6500, ra.Channels.ROYAL_BLUE, ...
%     ra.Channels.BLUE, ra.Channels.CYAN, ra.Channels.PC_CYAN_XEG, ...
%     ra.Channels.PC_CYAN_XQE, ra.Channels.LIME, ra.Channels.GREEN, ...
%     ra.Channels.AMBER, ra.Channels.PC_AMBER, ra.Channels.RED, ...
%     ra.Channels.PC_RED, ra.Channels.DEEP_RED, ra.Channels.FAR_RED};
% disabled_names = {ra.Channels.LIME};

active_mask = ~ismember(names, disabled_names);
% Для контроля вывод отключённых каналов:
if any(~active_mask)
    fprintf('Отключены каналы: %s\n', strjoin(names(~active_mask), ', '));
else
    fprintf('Все каналы активны.\n');
end
BIG = 1e6;                  % “большая” граница вместо Inf
w_max = BIG * ones(m,1);
w_max(~active_mask) = 0;    % отключённым — строго ноль
%% === CVX: L2-фит без масштаба ===
cvx_begin quiet
    variable w(m) nonnegative
    minimize( norm( W*(S*w - t), 2 ) + 1e-3*norm(w,1) )  % лёгкая L1-регуляризация
    subject to
        w >= 0;
        w <= w_max;
        % ВКЛЮЧИТЕ, если хотите смесь долей (при пик-нормировке столбцов S):
%         sum(w) == 1
cvx_end

w_opt = w;
fit   = S*w_opt;
resid = fit - t;

%% === Метрики и графики ===
rmse = sqrt(mean(resid.^2));
mae  = mean(abs(resid));
maxe = max(abs(resid));
fprintf('Status: %s\n', cvx_status);
fprintf('RMSE: %.4g, MAE: %.4g, Max|e|: %.4g\n', rmse, mae, maxe);

%% === ВЫВОД ВЕСОВ С ПОДПИСЯМИ И ПОМЕТКОЙ [OFF] ===
fprintf('\nОптимальные веса по каналам:\n');
for k = 1:m
    tag = '';
    if ~active_mask(k), tag = ' [OFF]'; end
    fprintf('%-14s : %.6f%s\n', names{k}, w_opt(k), tag);
end

% Таблица (удобно копировать/сохранять)
% T = table(string(names(:)), w_opt, logical(active_mask(:)), ...
%           'VariableNames', {'LED','Weight','Active'});
% disp(T);

%% === ГРАФИКИ ===
figure;
plot(lambda_led, t, '-', 'LineWidth', 1.5); hold on;
plot(lambda_led, fit, '--', 'LineWidth', 1.5);
xlabel('\lambda, nm'); ylabel('Relative power (norm to 1)');
legend('Target (norm)', 'Mixture S*w', 'Location', 'best'); grid on;
title('Подбор смеси LED к нормированному эталону');

figure;
plot(lambda_led, resid, 'LineWidth', 1);
xlabel('\lambda, nm'); ylabel('Error'); grid on;
title('Остаток: S*w - t');

%% === ПРОВЕРКА ===
check.checkSpectrum(lambda_led, t, S, Params, names, name2idx, w_opt, V);

%% === Оценка спектра ===

% подмена результирующей кривой планком для проверки расчёта Тц (ОК) и Ra (ORIG - не ОК, ITMO - ОК).
% fit = ra.planckSpd(lambda_led, 4500);

% базовая оценка смеси
[E, eItmo] = evaluateSpectrum(lambda_led, fit, t);

% fprintf('xy=(%.4f,%.4f), u''v''=(%.4f,%.4f)\n', E.color.xy(1),E.color.xy(2), E.color.uv(1),E.color.uv(2));
% fprintf('CCT=%.0f K, Δu''v''(до локуса)=%.5f\n', E.CCT.CCT, E.CCT.duv);
% if ~isnan(E.compare.dxy)
%     fprintf('Δxy(fit vs target)=%.5f, Δu''v''(fit vs target)=%.5f\n', E.compare.dxy, E.compare.duv);
% end

fprintf('--- Смесь (fit) ---\n');
fprintf('ORIG xy = (%.4f, %.4f), u''v'' = (%.4f, %.4f)\n', E.color.xy(1), E.color.xy(2), E.color.uv(1), E.color.uv(2));
fprintf('ITMO xy = (%.4f, %.4f), u''v'' = (%.4f, %.4f)\n', eItmo.color.xy(1), eItmo.color.xy(2), eItmo.color.uv(1), eItmo.color.uv(2));
fprintf('ORIG CCT = %.0f K, Δu''v''(до локуса) = %.5f\n', E.CCT.CCT, E.CCT.duv);
fprintf('ITMO CCT = %.0f K, Δu''v''(до локуса) = %.5f\n', eItmo.CCT.CCT);

if ~isempty(E.target.color)
    fprintf('\n--- Эталон (target) ---\n');
    fprintf('xy = (%.4f, %.4f), u''v'' = (%.4f, %.4f)\n', E.target.color.xy(1), E.target.color.xy(2), E.target.color.uv(1), E.target.color.uv(2));
    fprintf('CCT = %.0f K, Δu''v''(до локуса) = %.5f\n', E.target.CCT.CCT, E.target.CCT.duv);

    fprintf('\n--- Отличия fit vs target ---\n');
    fprintf('Δxy = %.5f, Δu''v'' = %.5f\n', E.compare.dxy, E.compare.duv);
end

CRI_fit    = ra.criRa(lambda_led, fit);
CRI_target = ra.criRa(lambda_led, t);
criItmoFit = itmo.criRaItmo(lambda_led, fit);
criItmoTarget = itmo.criRaItmo(lambda_led, t);

% Планковский излучатель при CCT смеси:
CCT_fit = 6500;
SPD_bb  = utils.planckSpd(lambda_led, CCT_fit);
CRI_bb  = ra.criRa(lambda_led, SPD_bb);
criItmoBb = itmo.criRaItmo(lambda_led, SPD_bb);

fprintf('ORIG ...............Ra(fit)=%.1f, Ra(target)=%.1f, Ra(Planck@CCT_fit)=%.1f\n', ...
        CRI_fit.Ra, CRI_target.Ra, CRI_bb.Ra);
fprintf('ITMO near Plank ....Ra(fit)=%.1f, Ra(target)=%.1f, Ra(Planck@CCT_fit)=%.1f\n', ...
        criItmoFit.nearPlank, criItmoTarget.nearPlank, criItmoBb.nearPlank);
fprintf('ITMO not near Plank Ra(fit)=%.1f, Ra(target)=%.1f, Ra(Planck@CCT_fit)=%.1f\n', ...
        criItmoFit.notNearPlnk, criItmoTarget.notNearPlnk, criItmoBb.notNearPlnk);

% --- BHL ---
BLH = ra.blueHazard(lambda_led, fit, 50, 1.0, 4.0);
fprintf('Blue-hazard radiance L_B = %.3e W·m^-2·sr^-1\n', BLH.LB_W_m2_sr); % ≤ 100 W·m⁻²·sr⁻¹ при t ≤ 100 с
fprintf('Blue-hazard irradiance E_B (on-axis @ %.1fm) = %.3e W·m^-2\n', 4.0, BLH.EB_W_m2); % ≤ 1 W·m⁻² для длительных (> 10⁴ с) экспозиций
fprintf('Relative BHL index = %.5f\n', BLH.BH_index);




