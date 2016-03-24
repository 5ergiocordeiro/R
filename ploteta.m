function [meta,teta]=ploteta(omega, sigma, er, mur)
% Calcula a imped�ncia em fun��o da frequ�mcia para diversos materiais.
% A frequ�ncia de entrada pode ser um vetor.
% As caracter�sticas consideradas dos materiais s�o a condutividade e a permissividade, que pode ser complexa.
% Retorna as imped�ncias em forma polar.
	% Constantes
	mu0 = 4 * pi * 1e-7;		% H/m
	ep0 = 8.854 * 1e-12;		% F/m
	% C�lculos
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
	