function [V,w,erro,iter]=solvegrid(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,perc,maxiter,Nrho,Nz,wmax,wmin)
% Calcula o potencial 'V' para a geometria ('Rext','Rint','Lext','Lint','rlim','zlim') e as condições ('perc','maxiter','Nrho','Nz','wmax','wmin') dadas.
% Chamada pela função principal ('probcil').
% Retorna também o maior fator de relaxação 'w' para o qual o cálculo foi bem sucedido, o 'erro' incorrido e o número de iterações 'iter' necessário. Se o cálculo não for bem sucedido para nenhum 'w' na faixa indicada, retorna 0 em 'w' e NaN para as demais saídas.
	% Tenta calcular o potencial com valores sucessivamente menores do fator de relaxação até ser bem sucedido.
	prec = perc * abs(Vext - Vint)/100;		% precisão absoluta
	step = (wmin - wmax)/9;					% negativo
	for w=wmax:step:wmin
		[V,erro,iter] = calcv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,prec,maxiter,Nrho,Nz,w);
		if erro <= prec
			% Cálculo convergiu
			return
		end
		% Cálculo não convergiu. Tentar com outro fator de relaxação menor.
		disp(sprintf("w = %f, erro = %f, iter = %d",w,erro,iter));
	end
	% Cálculo não convergiu para nenhum 'w' na faixa.
	w = 0;
	V = erro = iter = NaN;
	return
end
