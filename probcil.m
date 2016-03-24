%function [V]=probcil(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,prec,maxiter,Nrho,Nz,wmax,wmin)
% Calcula o perfil do potencial el�trico, do campo el�trico, a distribui��o de carga e a capacit�ncia para o problema de dois cilindros finitos conc�ntricos de raios 'Rext' e 'Rint', comprimentos 'Lext' e 'Lint' (dimens�es em mil�metros).
% Assume os potenciais 'Vext' e 'Vint' nas superf�cies.
% Assume o potencial nulo no infinito ('rlim' e 'zlim').
% Usa uma grade de 'Nrho' linhas e 'Nz' colunas. Exige precis�o m�nima percentual 'prec'. Limita o n�mero de itera��es em 'maxiter', por seguran�a.
% Tenta fatores de relaxa��o na faixa 'wmin' <= 'w' <= 'wmax'.
% Exemplo de uso:
% 	probcil(20,18,20,16,50,50,200,100,1,500,200,200,1.2,0.1)
	% Calcula as coordenadas das extremidades dos cilindros na grade escolhida
	RhoNint = floor(Nrho * Rint / rlim);		% Posi��o radial do cilindro interno
	RhoNext = floor(Nrho * Rext / rlim);		% Posi��o radial do cilindro externo
	TNint= Nz * Lint / zlim;					% Comprimento do cilindro interno
	TNext= Nz * Lext / zlim; 					% Comprimento do cilindro externo
	ZNint= floor(Nz/2 + TNint * [-0.5 0.5]);	% Posi��o das extremidades do cilindro interno
	ZNext= floor(Nz/2 + TNext * [-0.5 0.5]);	% Posi��o das extremidades do cilindro externo
	% Tenta calcular o potencial 'V' (apenas em dois quadrantes, por simetria)
	[V,w,erro,iter]=solvegrid(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,prec,maxiter,Nrho,Nz,wmax,wmin);
	if w == 0
		disp("O c�lculo n�o convergiu para nenhum fator de relaxa��o na faixa indicada.");
		return
	end
	disp(sprintf("O c�lculo convergiu com fator de relaxa��o %f ap�s %d itera��es", w, iter));
	% Plota 'V' (em apenas dois quadrantes)
	%plotv(V,RhoNext,RhoNint,ZNext,ZNint,rlim,zlim,Nrho,Nz);
	% Calcula e plota o campo el�trico 'E' (em apenas dois quadrantes)
	[Erho,Ez]=calce(V,rlim,zlim);
	%plote(Erho,Ez,rlim,zlim,Nrho,Nz);
	% Calcula a carga (em nC) e a capacit�ncia (em nF)
	qint = 2 * pi * Rint * 8.854e-9 * sum(Erho(RhoNint+1,ZNint(1):ZNint(2))) * zlim/Nz;
	qext = 2 * pi * Rext * 8.854e-9 * sum(Erho(RhoNext,ZNext(1):ZNext(2))) * zlim/Nz;
	erro = 100 * abs(1 - qext/qint);
	disp(sprintf("Carga nas superf�cies: externa = %f nC, interna = %f nC. Erro = %f ", qext, qint, erro));
	% rho = 0.5 * laplacianoc(V,rdelta)/(1 + rdelta^2) ;
	% q = 2 * pi * Rint * 8.854e-12 * sum(rho(RhoNint,ZNint(1):ZNint(2))) * zlim/Nz * rlim/Nrho;
	C = abs(qint/(Vext - Vint));
	disp(sprintf("Capacit�ncia estimada: %f nF", C));
%end
