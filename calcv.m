function [V,erro,i]=calcv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,prec,maxiter,Nrho,Nz,w)
% Tenta calcular o potencial 'V' para a geometria ('Rext','Rint','Lext','Lint','rlim','zlim','RhoNint','RhoNext','ZNint','ZNext') e as condições ('prec','maxiter','Nrho','Nz','w').
% Retorna o potencial calculado, o 'erro' incorrido e o número de iterações 'i' que foi necessário.
% Chamada pela função 'solvegrid'.
	% Inicializa 'V'
	[V,calc]=initv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,Nrho,Nz);
	% Inicializa as variáveis auxiliares
	last = 1e9;
	rdelta = ( rlim * Nz ) / ( zlim * Nrho ) ;
	% Aplica as condições de contorno de Dirichlet
	V(RhoNext,ZNext(1):ZNext(2)) = Vext; 	% Vext no cilindro externo
	V(RhoNint,ZNint(1):ZNint(2)) = Vint;	% Vint no clindro interno
	V(Nrho,:) = 0; V(:,Nz) = 0;				% 0 no "infinito"
	for i=1:maxiter
		% Aplica as condições de contorno de von Neumann
		V(Nrho-1,:) = 0; V(:,Nz-1) = 0;		% "infinito"
		V(1,:) = V(2,:) ; 					% simetria radial
		% Calcula o Laplaciano e o erro
		e = calc .* (0.5 * laplacianoc(V,rdelta)/(1 + rdelta^2));
		erro = max(max(abs(e)));
		if erro <= prec
			return
		end
		% Se o erro aumentou, para, pois o cálculo não vai convergir.
		if erro > last
			return
		end
		% Corrige os valores estimados e continua
		V = V + w * e;
		last = erro;
	end
end
