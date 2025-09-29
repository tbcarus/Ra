% Приведение амплитуды спектрального распределения к потоку излучения
function f = SPD_by_radiation(L, I, F)
I_result = F * I / trapz(L, I);
f = I_result;
end

