function s = truth_at(truth, t)
%TRUTH_AT  Эталонные величины в заданные моменты времени.
%
%   s = truth_at(truth, t)
%
%   Линейная интерполяция всех полей структуры truth (см. truth_from_state)
%   на сетку времени t. Нужна потому, что каждый датчик выдаёт данные со
%   своей частотой, не совпадающей с шагом модели автомобиля.

t = t(:);
f = {'X','Y','psi','vx','vy','r','ax','ay','VX','VY','V','s'};
s.t = t;
for i = 1:numel(f)
    s.(f{i}) = interp1(truth.t, truth.(f{i}), t, 'linear');
end

end
