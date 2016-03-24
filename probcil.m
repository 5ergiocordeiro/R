%function [V]=probcil(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,prec,maxiter,Nrho,Nz,wmax,wmin)
% Calcula o perfil do potencial elétrico, do campo elétrico, a distribuição de carga e a capacitância para o problema de dois cilindros finitos concêntricos de raios 'Rext' e 'Rint', comprimentos 'Lext' e 'Lint' (dimensões em milímetros).
% Assume os potenciais 'Vext' e 'Vint' nas superfícies.
% Assume o potencial nulo no infinito ('rlim' e 'zlim').
% Usa uma grade de 'Nrho' linhas e 'Nz' colunas. Exige precisão mínima percentual 'prec'. Limita o número de iterações em 'maxiter', por segurança.
% Tenta fatores de relaxação na faixa 'wmin' <= 'w' <= 'wmax'.
% Exemplo de uso:
% 	probcil(20,18,20,16,50,50,200,100,1,500,200,200,1.2,0.1)
	% Calcula as coordenadas das extremidades dos cilindros na grade escolhida
	RhoNint = floor(Nrho * Rint / rlim);		% Posição radial do cilindro interno
	RhoNext = floor(Nrho * Rext / rlim);		% Posição radial do cilindro externo
	TNint= Nz * Lint / zlim;					% Comprimento do cilindro interno
	TNext= Nz * Lext / zlim; 					% Comprimento do cilindro externo
	ZNint= floor(Nz/2 + TNint * [-0.5 0.5]);	% Posição das extremidades do cilindro interno
	ZNext= floor(Nz/2 + TNext * [-0.5 0.5]);	% Posição das extremidades do cilindro externo
	% Tenta calcular o potencial 'V' (apenas em dois quadrantes, por simetria)
	[V,w,erro,iter]=solvegrid(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,prec,maxiter,Nrho,Nz,wmax,wmin);
	if w == 0
		disp("O cálculo não convergiu para nenhum fator de relaxação na faixa indicada.");
		return
	end
	disp(sprintf("O cálculo convergiu com fator de relaxação %f após %d iterações", w, iter));
	% Plota 'V' (em apenas dois quadrantes)
	%plotv(V,RhoNext,RhoNint,ZNext,ZNint,rlim,zlim,Nrho,Nz);
	% Calcula e plota o campo elétrico 'E' (em apenas dois quadrantes)
	[Erho,Ez]=calce(V,rlim,zlim);
	%plote(Erho,Ez,rlim,zlim,Nrho,Nz);
	% Calcula a carga (em nC) e a capacitância (em nF)
	qint = 2 * pi * Rint * 8.854e-9 * sum(Erho(RhoNint+1,ZNint(1):ZNint(2))) * zlim/Nz;
	qext = 2 * pi * Rext * 8.854e-9 * sum(Erho(RhoNext,ZNext(1):ZNext(2))) * zlim/Nz;
	erro = 100 * abs(1 - qext/qint);
	disp(sprintf("Carga nas superfícies: externa = %f nC, interna = %f nC. Erro = %f ", qext, qint, erro));
	% rho = 0.5 * laplacianoc(V,rdelta)/(1 + rdelta^2) ;
	% q = 2 * pi * Rint * 8.854e-12 * sum(rho(RhoNint,ZNint(1):ZNint(2))) * zlim/Nz * rlim/Nrho;
	C = abs(qint/(Vext - Vint));
	disp(sprintf("Capacitância estimada: %f nF", C));
%end
