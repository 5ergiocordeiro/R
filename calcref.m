function [Gamma,tau]=calcref(modo, omega, thetai, er1, mur1, er2, mur2)
	% Constantes
	mu0 = 4 * pi * 1e-7;		% H/m
	ep0 = 8.854 * 1e-12;		% F/m
	% Cálculos
	mu = mu0 * [mur1 mur2];
	epsilon = ep0 * [er1 er2];
	k = omega * (mu .* epsilon).^0.5;
	eta = mu * omega ./ k;
	thetat = asin(k(1)/k(2).*sin(thetai));
	if modo == "n"
		aux = eta(2) .* cos(thetai) + eta(1) .* cos(thetat);
		Gamma = (eta(2) .* cos(thetai) - eta(1) .* cos(thetat))./aux;
	else
		aux = eta(1) .* cos(thetai) + eta(2) .* cos(thetat);
		Gamma = (eta(1) .* cos(thetai) - eta(2) .* cos(thetat))./aux;
	end
	tau = (2 * eta(2) .* cos(thetai))./aux;
end