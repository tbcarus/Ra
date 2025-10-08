% Приведение спектрального распределения к заданному световому потоку
function f = spdByFlux(L, I, F, V)
I_result = F * I / (683 * trapz(L, I.*V));
f = I_result;
end

