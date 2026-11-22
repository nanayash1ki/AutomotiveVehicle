function truth = truth_from_state(t, x)
%TRUTH_FROM_STATE  Эталон для имитаторов из тестовой траектории make_test_truth.
%
%   truth = truth_from_state(t, x)
%
%   ВНИМАНИЕ: здесь vx, vy — в СВЯЗАННОЙ СК. Для модели автомобиля
%   (run_truth, где vx, vy в мировой СК) используйте truth_from_run_truth.
%
%   t — время, N×1, с;
%   x — вектор состояния модели автомобиля, N×6:
%       [X Y psi vx vy r],
%       X, Y  — координаты центра масс в навигационной СК, м;
%       psi   — курс, рад (непрерывный, без приведения к ±pi);
%       vx,vy — продольная и поперечная скорости в связанной СК, м/с;
%       r     — угловая скорость рыскания, рад/с.
%
%   Возвращает структуру с полями t, X, Y, psi, vx, vy, r, а также:
%       ax, ay — кажущиеся ускорения в связанной СК (то, что измеряет
%                акселерометр при движении по горизонтальной дороге), м/с^2;
%       VX, VY — составляющие скорости в навигационной СК, м/с;
%       V      — модуль скорости, м/с;
%       s      — пройденный путь, м.

t = t(:);
X = x(:,1);  Y = x(:,2);  psi = unwrap(x(:,3));
vx = x(:,4); vy = x(:,5); r = x(:,6);

% Ускорение в связанной (вращающейся) СК:
%   ax = dvx/dt - r*vy,  ay = dvy/dt + r*vx
truth.t   = t;
truth.X   = X;
truth.Y   = Y;
truth.psi = psi;
truth.vx  = vx;
truth.vy  = vy;
truth.r   = r;
truth.ax  = gradient(vx, t) - r.*vy;
truth.ay  = gradient(vy, t) + r.*vx;
truth.VX  = vx.*cos(psi) - vy.*sin(psi);
truth.VY  = vx.*sin(psi) + vy.*cos(psi);
truth.V   = hypot(vx, vy);
truth.s   = cumtrapz(t, abs(vx));

end
