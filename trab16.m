% Constantes
epsilon = 8.85e-12; 			% F/m
mu = 1.26e-6; 					% H/m
global w = sqrt(epsilon/mu);	% ohms
% Parâmetros do problema (em metros)
global R = 1;					% raio do cilindro
lambda = [0.01 1 100];			% comprimento de onda
% Outros valores
global infty = 5;				% controla quantas componentes serão integradas
step = 1;						% resolução angular, em graus
nrhos = 30;						% número de valores diferentes de rho a calcular
Eo = Ho = 1;					% intensidade do campo incidente

function [Jlinha] = dbessel(n, arg)
% Retorna a derivada da função de Bessel de primeira espécie
	Jlinha = besselj(n - 1, arg, 0) - (n + 1) / arg * besselj(n, arg, 0);
end

function [Hlinha] = dhankel(n, arg)
% Retorna a derivada da função de Hankel de segunda espécie
	Hlinha = (besselh(n + 1, 2, arg, 0) - besselh(n - 1, 2, arg, 0)) / 2;
end

function [Hr, Hp, Ez, Ezero] = tmcalc(siz,phi,arg1,arg2,infty,Eo,Ho,l,w,rho)
% Calcula os campos para o modo TM
	sumEz = sumHr = sumHp = 0;
	for n =-infty:infty
		sumEz += j^(-n) * e^(j * n * phi) * (besselj(n,arg1,0) - besselj(n,arg2,0) * besselh(n,2,arg1,0)/besselh(n,2,arg2,0));
		sumHr += n * j^(1-n) * e^(j * n * phi) * (besselj(n,arg1,0) - besselj(n,arg2,0) * besselh(n,2,arg1,0)/besselh(n,2,arg2,0));
		sumHp += j^(-n) * e^(j * n * phi) * (dbessel(n,arg1) - besselj(n,arg2,0) * dhankel(n,arg1)/besselh(n,2,arg2,0));
	end
	Hr = j * Eo * l / (2 * pi * rho * w) * sumHr;
	Hp = - j * Eo / w * sumHp;
	Ez = Eo * sumEz;
	Ezero = 0;
end

function [Hz, Hzero, Er, Ep] = tecalc(siz,phi,arg1,arg2,infty,Eo,Ho,l,w,rho)
% Calcula os campos para o modo TE
	sumEr = sumEp = sumHz = 0;
	for n =-infty:infty
		sumEr += n * j^(1-n) * e^(j * n * phi) * (besselj(n,arg1,0) - dbessel(n,arg2) * besselh(n,2, arg1, 0) / dhankel(n,arg2));
		sumEp += j^(-n) * e^(j * n * phi) * (dbessel(n,arg1) - dbessel(n,arg2) * dhankel(n,arg1) / dhankel(n,arg2));
		sumHz += j^(-n) * e^(j * n * phi) * (besselj(n,arg1,0) - dbessel(n,arg2) * besselh(n,2,arg1,0) / dhankel(n,arg2));
	end
	Hz = Ho * sumHz;
	Hzero = 0;
	Er = - j * Ho * l * w / (2 * pi * rho) * sumEr;
	Ep = j * Ho * w * sumEp;
end

function [H1, H2, E1, E2] = calccyl(p, l, R, prec, step, Eo, Ho)
% Calcula os campos para o modo e o comprimento de onda indicados
	global infty;
	global w;
	siz = floor(360/step);
	tphi = linspace(-pi,pi,siz);
	arg2 = 2 * pi * R / l;
	H1 = H2 = E1 = E2 = zeros(prec, siz);
	% Para cada valor de rho
	for q = 1:prec
		rho = R * prec;
		arg1 = 2 * pi * rho / l;
		% Para cada valor de phi
		for k = 1:siz
			phi = tphi(k);
			if p == 1 	% modo TMz			
				[sumH1,sumH2,sumE1,sumE2] = tmcalc(siz,phi,arg1,arg2,infty,Eo,Ho,l,w,rho);
			else 		% modo TEz
				[sumH1,sumH2,sumE1,sumE2] = tecalc(siz,phi,arg1,arg2,infty,Eo,Ho,l,w,rho);
			end
			E1(q,k) = sumE1;
			E2(q,k) = sumE2;
			H1(q,k) = sumH1;
			H2(q,k) = sumH2;
		end
	end
end

function plotcyl(H1, H2, E1, E2, p, i, nrhos, step, prec)
% plota os gráficos
	global R;
	siz = floor(360/step);
	tphi = linspace(-pi,pi,siz);
	tit = [
		"Onda curta";
		"Onda media";
		"Onda longa"
		];
	modos = [
		"TMz";
		"TEz"
		];
	% Corrente na superfície do cilindro
	curdens = abs(H2(1,:));
	figure(1); clf;
		plot(tphi, curdens);
		title(sprintf("%s - Modo %s",tit(i,:),modos(p,:)));
		ylabel("Densidade de corrente superficial (A/m)");
		xlabel("Angulo (rd)");
	% Variação da intensidade do campo
	Eabs = sqrt(E1.^2 + E2.^2);
	Habs = sqrt(H1.^2 + H2.^2);
	% na superfície do cilindro
	figure(2); clf;
		plot(tphi,Eabs(1,:),tphi,Habs(1,:));
		title(sprintf("%s - Modo %s",tit(i,:),modos(p,:)));
		xlabel("Angulo (rd)");
		legend("E","H");
	% conforme a distância ao cilindro
	q = 1:prec;
	trho = R * q;
	figure(3); clf;
		plot(trho,Eabs(:,siz/2),trho,Habs(:,siz/2));
		title(sprintf("%s - Modo %s",tit(i,:),modos(p,:)));
		xlabel("Distancia (m)");
		legend("E","H");
end
		
% Para cada modo
for p = 1:2
	% Para cada comprimento de onda
	for i=1:3
		% Calcula os campos
		l = lambda(i);
		[H1, H2, E1, E2] = calccyl(p, l, R, nrhos, step, Eo, Ho);
		% Plota os resultados
		plotcyl(H1, H2, E1, E2, p, i, nrhos, step, nrhos);
	end
end