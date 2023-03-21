function [rf] = FGKepler_trA(mu, r0, r0_dot, delta_theta)

r0_norm = norm(r0);

% Angular momentum
h_vec = cross(r0, r0_dot);
h = norm(h_vec);
p = h^2/mu;

sigma0 = dot(r0,r0_dot)/sqrt(mu);

rf = zeros(3,length(delta_theta));

for i = 1:length(delta_theta)
    rf_norm = p*r0_norm/(r0_norm + (p - r0_norm)*cos(delta_theta(i)) - sqrt(p) * sigma0 * sin(delta_theta(i)));
    
    F = 1 - rf_norm/p * (1 - cos(delta_theta(i)));
    G = (rf_norm*r0_norm)/sqrt(mu*p) * sin(delta_theta(i));
    
    rf_aux = F*r0 + G*r0_dot;
    rf(1,i) = rf_aux(1);
    rf(2,i) = rf_aux(2);
    rf(3,i) = rf_aux(3);
end

end