function [t, x] = make_test_truth(dt, T)
%MAKE_TEST_TRUTH  Тестовая траектория для отладки (без модели динамики).
%
%   [t, x] = make_test_truth(dt, T)
%
%   Кинематическая траектория без бокового скольжения (vy = 0):
%   стоянка 5 с, разгон до 15 м/с, повороты влево и вправо, змейка,
%   торможение до 8 м/с и поворот. Часть манёвров попадает в интервал
%   пропадания ГНСС 45–85 с.
%   Выход: x = [X Y psi vx vy r], vx, vy — в СВЯЗАННОЙ СК
%   (обрабатывается функцией truth_from_state).

if nargin < 1, dt = 0.005; end
if nargin < 2, T  = 150;   end

d2r = pi/180;
t = (0:dt:T).';

% продольная скорость: кусочно-линейный профиль
vx = interp1([0 5 15 105 115 T], [0 0 15 15 8 8], t, 'linear');

% угловая скорость рыскания: трапеции [начало конец скорость °/с]
turns = [ 30  40   9.0;     % поворот влево на 90°
          55  65  -9.0;     % поворот вправо на 90°
          75  85   4.5;     % змейка
          85  95  -4.5;
         120 135   6.0];    % поворот влево на 90°
r = zeros(size(t));
ramp = 1.0;                 % время нарастания угловой скорости, с
for i = 1:size(turns,1)
    up   = (t - turns(i,1) + ramp/2)/ramp;
    down = (turns(i,2) + ramp/2 - t)/ramp;
    r = r + turns(i,3)*d2r*min(max(min(up, down), 0), 1);
end

psi = cumtrapz(t, r);
X   = cumtrapz(t, vx.*cos(psi));
Y   = cumtrapz(t, vx.*sin(psi));
vy  = zeros(size(t));

x = [X Y psi vx vy r];

end
