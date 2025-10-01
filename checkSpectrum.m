function checkSpectrum(L, sun, S, Params, names, name2idx, w_opt, V)
% L - массив длин волн
% sun - эталонный солнечный спектр
% S - таблица спектров СИД
% Params - параметры СИД
% names - имена каналов
% name2idx - мапа имя->индекс в таблице спектров
% w_opt - расчитанные коэффициенты каналов
% 
% Проверка спектра: 
% [V] - его составляющие
% [V] - результирующая мощность потока
% [V] - световой поток
% [] - потребляемая мощность
% [] - эффективность
% Вывод графика:
% [V] - спектры СИД, умноженные на весовые коэффициенты
% [V] - суммарный спектр
% [V] - эталонный спектр

%% Вывод графиков
figure
plot(L, w_opt(name2idx(ra.Channels.WARM2700))*S(:,name2idx(ra.Channels.WARM2700)), 'Color', [220/255 202/255 47/255]);
grid on;
hold on;
% Можно ещё цвет графиков в таблицу загнать, тогда можно пробежаться по
% индексам names в цикле, а так каждый раз вручную спектры добавлять
plot(L, w_opt(name2idx(ra.Channels.COOL6500))*S(:,name2idx(ra.Channels.COOL6500)), 'Color', [44/255 165/255 222/255]);
plot(L, w_opt(name2idx(ra.Channels.CYAN))*S(:,name2idx(ra.Channels.CYAN)), 'Color', [20/255 215/255 233/255]);
plot(L, w_opt(name2idx(ra.Channels.PC_CYAN_XEG))*S(:,name2idx(ra.Channels.PC_CYAN_XEG)), 'Color', [120/255 215/255 233/255]);
plot(L, w_opt(name2idx(ra.Channels.PC_CYAN_XQE))*S(:,name2idx(ra.Channels.PC_CYAN_XQE)), 'Color', [52/255 175/255 189/255]);
plot(L, w_opt(name2idx(ra.Channels.ROYAL_BLUE))*S(:,name2idx(ra.Channels.ROYAL_BLUE)), 'Color', [0/255 11/255 153/255]);
plot(L, w_opt(name2idx(ra.Channels.BLUE))*S(:,name2idx(ra.Channels.BLUE)), 'Color', [0/255 0/255 255/255]);
plot(L, w_opt(name2idx(ra.Channels.FAR_RED))*S(:,name2idx(ra.Channels.FAR_RED)), 'Color', [90/255 0/255 0/255]);
plot(L, w_opt(name2idx(ra.Channels.DEEP_RED))*S(:,name2idx(ra.Channels.DEEP_RED)), 'Color', [155/255 0/255 0/255]);
plot(L, w_opt(name2idx(ra.Channels.RED))*S(:,name2idx(ra.Channels.RED)), 'Color', [255/255 0/255 0/255]);
plot(L, w_opt(name2idx(ra.Channels.GREEN))*S(:,name2idx(ra.Channels.GREEN)), 'Color', [0/255 255/255 0/255]);
plot(L, w_opt(name2idx(ra.Channels.LIME))*S(:,name2idx(ra.Channels.LIME)), 'Color', [85/255 220/255 85/255]);
plot(L, w_opt(name2idx(ra.Channels.PC_AMBER))*S(:,name2idx(ra.Channels.PC_AMBER)), 'Color', [255/255 198/255 0/255]);

resultSpectrum = ... % Тут ожно не перечислять, а пробежаться сразу по индексам names, а так каждый раз вручную спектры добавлять
    w_opt(name2idx(ra.Channels.WARM2700))*S(:,name2idx(ra.Channels.WARM2700)) + ...
    w_opt(name2idx(ra.Channels.COOL6500))*S(:,name2idx(ra.Channels.COOL6500)) + ...
    w_opt(name2idx(ra.Channels.CYAN))*S(:,name2idx(ra.Channels.CYAN)) + ...
    w_opt(name2idx(ra.Channels.PC_CYAN_XEG))*S(:,name2idx(ra.Channels.PC_CYAN_XEG)) + ...
    w_opt(name2idx(ra.Channels.PC_CYAN_XQE))*S(:,name2idx(ra.Channels.PC_CYAN_XQE)) + ...
    w_opt(name2idx(ra.Channels.ROYAL_BLUE))*S(:,name2idx(ra.Channels.ROYAL_BLUE)) + ...
    w_opt(name2idx(ra.Channels.BLUE))*S(:,name2idx(ra.Channels.BLUE)) + ...
    w_opt(name2idx(ra.Channels.FAR_RED))*S(:,name2idx(ra.Channels.FAR_RED)) + ...
    w_opt(name2idx(ra.Channels.DEEP_RED))*S(:,name2idx(ra.Channels.DEEP_RED)) + ...
    w_opt(name2idx(ra.Channels.RED))*S(:,name2idx(ra.Channels.RED)) + ...
    w_opt(name2idx(ra.Channels.GREEN))*S(:,name2idx(ra.Channels.GREEN)) + ...
    w_opt(name2idx(ra.Channels.LIME))*S(:,name2idx(ra.Channels.LIME)) + ...
    w_opt(name2idx(ra.Channels.PC_AMBER))*S(:,name2idx(ra.Channels.PC_AMBER));

plot(L, resultSpectrum, 'red', 'LineWidth', 2, 'LineStyle', '--');
plot(L, sun, 'black')

%% Мощность потока
Pe = trapz(L, resultSpectrum);
Pv = 683 * trapz(L, resultSpectrum.*V');
fprintf('%s: %f Вт\n', 'Мощность потока: ', Pe);
fprintf('%s: %f лм\n', 'Световой поток: ', Pv);

%% Потребляемая мощность и эффективность
% Тут ожно не перечислять, а пробежаться сразу по индексам names, а так каждый раз вручную спектры добавлять
P = w_opt(name2idx(ra.Channels.WARM2700))*ra.getParam(names, Params, ra.Channels.WARM2700, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.COOL6500))*ra.getParam(names, Params, ra.Channels.COOL6500, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.CYAN))*ra.getParam(names, Params, ra.Channels.CYAN, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.PC_CYAN_XEG))*ra.getParam(names, Params, ra.Channels.PC_CYAN_XEG, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.PC_CYAN_XQE))*ra.getParam(names, Params, ra.Channels.PC_CYAN_XQE, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.ROYAL_BLUE))*ra.getParam(names, Params, ra.Channels.ROYAL_BLUE, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.BLUE))*ra.getParam(names, Params, ra.Channels.BLUE, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.FAR_RED))*ra.getParam(names, Params, ra.Channels.FAR_RED, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.DEEP_RED))*ra.getParam(names, Params, ra.Channels.DEEP_RED, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.RED))*ra.getParam(names, Params, ra.Channels.RED, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.GREEN))*ra.getParam(names, Params, ra.Channels.GREEN, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.LIME))*ra.getParam(names, Params, ra.Channels.LIME, ra.LedParam.P) + ...
    w_opt(name2idx(ra.Channels.PC_AMBER))*ra.getParam(names, Params, ra.Channels.PC_AMBER, ra.LedParam.P);
fprintf('%s: %f Вт\n', 'Потребляемая мощность: ', P);
fprintf('%s: %f\n', 'КПД: ', Pe/P);
fprintf('%s: %f лм/Вт\n', 'Эффективность: ', Pv/P);

end

