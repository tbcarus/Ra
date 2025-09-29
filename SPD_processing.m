% ќбрыботка спектра:
% - исключение отрицательных значений
% - интерпол¤ци¤
% - исключение неопределЄнностей

function f = SPD_processing(L_RAW, I_RAW, L_interp)
I_RAW(I_RAW < 0) = 0;
I = interp1(L_RAW, I_RAW, L_interp, 'linear');
I(isnan(I)) = 0;
I = I/max(I);
f = I;
end