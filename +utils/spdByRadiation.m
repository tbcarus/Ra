% Приведение спектрального распределения к заданному потоку излучения
function f = SPD_by_radiation(L, I, F)
I_result = F * I / trapz(L, I);
f = I_result;
end

