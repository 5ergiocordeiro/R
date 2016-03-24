function [meta,teta]=ploteta(omega, sigma, er, mur)
% Calcula a impedância em função da frequêmcia para diversos materiais.
% A frequência de entrada pode ser um vetor.
% As características consideradas dos materiais são a condutividade e a permissividade, que pode ser complexa.
% Retorna as impedâncias em forma polar.
	% Constantes
	mu0 = 4 * pi * 1e-7;		% H/m
	ep0 = 8.854 * 1e-12;		% F/m
	% Cálculos
	mu = mu0 * mur;
	epsilon = ep0 * er;
	n = size(sigma, 2);
	m = size(omega, 2);
	eta = zeros(n, m);
	for i = 1:n
		k = omega .* (mu(i) * (epsilon(i) + j * sigma(i) ./ omega)).^0.5;
		eta(i, :) = j * mu(i) * omega ./ k;
	end
	meta = abs(eta);
	teta = arg(eta);
end
	