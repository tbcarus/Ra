function val = getParam(names, Params, channelName, paramCol)
% Быстрый и безопасный доступ к параметру конкретного канала
% Пример: I = ra.get_param(names, Params, ra.Channels.CYAN, ra.LedParam.I);
idx = find(strcmp(names, char(channelName)), 1);
assert(~isempty(idx), 'Канал "%s" не найден', channelName);
val = Params(idx, paramCol);
end
