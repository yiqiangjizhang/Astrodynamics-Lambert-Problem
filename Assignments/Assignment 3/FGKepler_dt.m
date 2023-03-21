function [rf, rf_dot, delta_E_sol, F, G] = FGKepler_dt(mu, r0, r0_dot, delta_t)

r0_norm = norm(r0); 
r0_dot_norm = norm(r0_dot);

a = mu/(2*mu/r0_norm - r0_dot_norm^2);

n = sqrt(mu/a^3);

delta_M = n*delta_t;

sigma_0 = (dot(r0,r0_dot))/sqrt(mu);

fun = @(delta_E) delta_E - (1 - r0_norm/a)*sin(delta_E) - sigma_0/sqrt(a) * (cos(delta_E)- 1) - delta_M; % function
x0 = delta_M; % initial point
delta_E_sol = fzero(fun,x0);

% Check
fun(delta_E_sol);

F = 1 - a/r0_norm * (1 - cos(delta_E_sol));
G = delta_t + sqrt(a^3/mu) * (sin(delta_E_sol) - delta_E_sol);

rf = F*r0 + G*r0_dot;
rf_norm = norm(rf);

G_dot = 1 - a/rf_norm*(1 - cos(delta_E_sol));

rf_dot = 1/G*(G_dot*rf - r0);

end