%{
Calcula a parte real e a parte imaginária da constante de propagação (k) em função da frequência e das propriedades do meio.
Aceita que a permissividade e a permeabilidade sejam complexas.
Aceita um vetor de frequências na entrada.
%}
function [alpha,beta] = calcparms(omega, sigma, mur, er)
	% Constantes
	mu0 = 4 * pi * 1e-7;		% H/m
	ep0 = 8.854 * 1e-12;		% F/m
	% Cálculos
	epsilon = er * ep0;
	mu = mur * mu0;
	aux = (1 + (sigma./(omega * epsilon).^2)).^0.5;
	alpha = omega .* (mu * epsilon / 2 .* (aux - 1)).^0.5;
	beta = omega .* (mu * epsilon / 2 .* (aux + 1)).^0.5;
end
