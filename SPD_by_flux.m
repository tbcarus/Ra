% ѕриведение амплитуды спектрального распределени€ к световому потоку
function f = SPD_by_flux(L, I, F, V)
I_result = F * I / (683 * trapz(L, I.*V));
f = I_result;
end

