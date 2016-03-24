function [V,w,erro,iter]=solvegrid(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,perc,maxiter,Nrho,Nz,wmax,wmin)
% Calcula o potencial 'V' para a geometria ('Rext','Rint','Lext','Lint','rlim','zlim') e as condi��es ('perc','maxiter','Nrho','Nz','wmax','wmin') dadas.
% Chamada pela fun��o principal ('probcil').
% Retorna tamb�m o maior fator de relaxa��o 'w' para o qual o c�lculo foi bem sucedido, o 'erro' incorrido e o n�mero de itera��es 'iter' necess�rio. Se o c�lculo n�o for bem sucedido para nenhum 'w' na faixa indicada, retorna 0 em 'w' e NaN para as demais sa�das.
	% Tenta calcular o potencial com valores sucessivamente menores do fator de relaxa��o at� ser bem sucedido.
	prec = perc * abs(Vext - Vint)/100;		% precis�o absoluta
	step = (wmin - wmax)/9;					% negativo
	for w=wmax:step:wmin
		[V,erro,iter] = calcv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,prec,maxiter,Nrho,Nz,w);
		if erro <= prec
			% C�lculo convergiu
			return
		end
		% C�lculo n�o convergiu. Tentar com outro fator de relaxa��o menor.
		disp(sprintf("w = %f, erro = %f, iter = %d",w,erro,iter));
	end
	% C�lculo n�o convergiu para nenhum 'w' na faixa.
	w = 0;
	V = erro = iter = NaN;
	return
end
