function truth = truth_from_run_truth(tr)
%TRUTH_FROM_RUN_TRUTH  Эталон для имитаторов из выхода run_truth.
%
%   truth = truth_from_run_truth(tr)
%
%   В модели автомобиля скорости vx, vy заданы в МИРОВОЙ СК, а ax_b, ay_b —
%   ускорения центра масс в связанной СК (именно их измеряет акселерометр
%   на горизонтальной дороге). Функция приводит выход run_truth к формату,
%   который понимают имитаторы датчиков:
%       X, Y   — координаты, м;        psi — курс, рад (непрерывный);
%       vx, vy — скорости в СВЯЗАННОЙ СК, м/с (vx — продольная, её видит одометр);
%       r      — угловая скорость рыскания, рад/с;
%       ax, ay — ускорения в связанной СК, м/с^2;
%       VX, VY — скорости в мировой (навигационной) СК, м/с (их выдаёт ГНСС);
%       V      — модуль скорости, м/с;  s — пройденный путь, м.

psi = unwrap(tr.phi(:));
VX = tr.vx(:);  VY = tr.vy(:);

truth.t   = tr.t(:);
truth.X   = tr.x(:);
truth.Y   = tr.y(:);
truth.psi = psi;
truth.vx  =  VX.*cos(psi) + VY.*sin(psi);
truth.vy  = -VX.*sin(psi) + VY.*cos(psi);
truth.r   = tr.dphi(:);
truth.ax  = tr.ax_b(:);
truth.ay  = tr.ay_b(:);
truth.VX  = VX;
truth.VY  = VY;
truth.V   = hypot(VX, VY);
truth.s   = cumtrapz(truth.t, abs(truth.vx));

end
