function b = gauss_markov(n, dt, sigma, tau, m)
%GAUSS_MARKOV  Марковский процесс первого порядка (медленный дрейф).
%
%   b = gauss_markov(n, dt, sigma, tau)      — n отсчётов, один канал
%   b = gauss_markov(n, dt, sigma, tau, m)   — n×m, m независимых каналов
%
%   Дискретная модель:  b(k) = phi*b(k-1) + w(k),  phi = exp(-dt/tau),
%   СКО шума w выбрано так, чтобы установившееся СКО процесса было sigma.
%   Начальное значение берётся из установившегося распределения N(0, sigma^2).
%
%   sigma — установившееся СКО; tau — время корреляции, с.

if nargin < 5, m = 1; end

phi = exp(-dt/tau);
w = sigma*sqrt(1 - phi^2)*randn(n, m);
w(1, :) = sigma*randn(1, m);          % начальное значение b(1)
b = filter(1, [1 -phi], w);           % b(k) = phi*b(k-1) + w(k) по столбцам

end
