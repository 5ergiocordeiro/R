% Par�metros do problema (em metros)
R = 1;					% raio do cilindro
lambda = [0.3 1 3];		% comprimento de onda
% Outros valores
infty = 20;				% controla quantas componentes ser�o integradas
step = 10;				% resolu��o angular
prec = 10;				% resolu��o radial
scale = pi()/step;		% convers�o de unidades
Eo = 1;					% intensidade do campo incidente

% Fun��es auxiliares
function [b,c] = Calc_bc(m, alpha)
% Calcula os coeficientes an, bn e cn, com 0 <= n <= m
	b = c = zeros(m,1);
	for n = 1:m
		a = j^(-n) * (2 * n + 1) / (n * (n + 1) );
		b(n) = - a * Jotalinha(n, alpha) / Haga2linha(n, alpha);
		c(n) = - a * Jota(n, alpha) / Haga2(n, alpha);
	end
end

function [Jlinha] = Jotalinha(n, alpha)
% Calcula a derivada da fun��o de Bessel de primeira esp�cie esf�rica modificada
	Jlinha = Jota(n - 1, alpha) - n / alpha * Jota(n, alpha);
end

function [H2linha] = Haga2linha(n, alpha)
% Calcula a derivada da fun��o de Hankel de segunda esp�cie esf�rica modificada
	H2linha = Haga2(n - 1, alpha) - n / alpha *  Haga2(n, alpha);		
end

function [H2linhalinha] = Haga2linhalinha(n, alpha)
% Calcula a derivada segunda da fun��o de Hankel de segunda esp�cie esf�rica modificada
	H2linhalinha = alpha * Haga2(n - 2, alpha) + 2 / alpha * Haga2(n - 1, alpha) + n ^2 / alpha * Haga2(n, alpha);
end

function [Hn2] = Haga2(n, alpha)
% Calcula a fun��o de Hankel de segunda esp�cie esf�rica modificada
	Hn2 = sqrt(pi/2 * alpha) * besselh(n + 1/2, 2, alpha, 0);  
end

function [Jn] = Jota(n, alpha)
% Calcula a fun��o de Bessel de primeira esp�cie esf�rica modificada
	Jn = sqrt(pi/2 * alpha) * besselj(n + 1/2, alpha, 0);  
end

function [pl] = Pl(n, theta)
% Calcula a fun��o de Legendre associada
	if n <= 0
		pl = 0;
	else
		Pn = legendre(n, cos(theta));
		pl = Pn(end);
	end
end

function [p1] = Pe1(n, theta)
% Calcula a fun��o de Legendre associada
	p1 = n / (theta^2 - 1) * (theta * Pl(n, theta) - Pl(n - 1, theta));
end

function [p1] = Pelinha1(n, theta)
% Calcula a fun��o de Legendre associada
	theta2 = theta^2;
	c1 = (n - 1) * theta2 - 1;
	c2 = 2 * (n +1) * theta;
	p1 = n / sin(theta) / (theta2 - 1)^2 * (c1 * Pl(n, theta) + c2 * Pl(n - 1, theta) + n * Pl(n - 2, theta));
end

function [E] = Calc_E(m, alpha, theta, phi, b, c)
% Calcula o campo el�trico espalhado na posi��o dada
	costheta = cos(theta);
	sentheta = sin(theta);
	cosphi = cos(phi);
	senphi = sin(phi);
	sEr = sEphi1 = sEphi2 = sEtheta1 = sEtheta2 = 0;
	for n = 1:m
		sEr = sEr + b(n) * (Haga2linhalinha(n, alpha) + Haga2(n, alpha)) * Pe1(n, theta);
		if sentheta != 0
			sEtheta1 = sEtheta1 + j * b(n) * Haga2linha(n, alpha) * Pelinha1(n, theta);
			sEtheta2 = sEtheta2 + c(n) * Haga2(n, alpha) * Pe1(n, theta);
			sEphi1 = sEphi1 + j * b(n) * Haga2linha(n, alpha) * Pe1(n, theta);
			sEphi2 = sEphi2 + c(n) * Haga2(n, alpha) * Pelinha1(n, theta);
		end
	end		
	Er = -j * cosphi * sEr;
	if sentheta != 0
		Etheta = cosphi / alpha * (sentheta * sEtheta1 - sEtheta2 / sentheta);
		Ephi = senphi / alpha * (sEphi1 / sentheta - sEphi2 * sentheta); 
	else
		Etheta = Ephi = 0;
	end
	E = [Er, Etheta, Ephi];
end

function [Er,Etheta,Ephi] = Calc_Esfera(l, prec, step, scale, m)
% Calcula as componentes do campo el�trico em cada ponto do espa�o, para a geometria dada
	Er=zeros(prec,step,step);
	Etheta=zeros(prec,step,step);
	Ephi=zeros(prec,step,step);
	for qr = 1:prec
		alpha = 2 * pi * qr * l ;
		for  qphi = 1:step
			phi = qphi * 2 * scale;
			for qtheta = 1:step
				theta = qtheta * scale;
				[b,c] = Calc_bc(m, alpha);
				E = Calc_E(m, alpha, theta, phi, b, c) ;
				Er(qr, qphi, qtheta) = E(1);
				Etheta(qr, qphi, qtheta) = E(2);
				Ephi(qr, qphi, qtheta) = E(3);				
			end
		end
	end
end

for i=1:3
	[Er,Etheta,Ephi] = Calc_Esfera(R / lambda(i), prec, step, scale, infty);
end	
