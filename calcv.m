function [V,erro,i]=calcv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,prec,maxiter,Nrho,Nz,w)
% Tenta calcular o potencial 'V' para a geometria ('Rext','Rint','Lext','Lint','rlim','zlim','RhoNint','RhoNext','ZNint','ZNext') e as condi��es ('prec','maxiter','Nrho','Nz','w').
% Retorna o potencial calculado, o 'erro' incorrido e o n�mero de itera��es 'i' que foi necess�rio.
% Chamada pela fun��o 'solvegrid'.
	% Inicializa 'V'
	[V,calc]=initv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,Nrho,Nz);
	% Inicializa as vari�veis auxiliares
	last = 1e9;
	rdelta = ( rlim * Nz ) / ( zlim * Nrho ) ;
	% Aplica as condi��es de contorno de Dirichlet
	V(RhoNext,ZNext(1):ZNext(2)) = Vext; 	% Vext no cilindro externo
	V(RhoNint,ZNint(1):ZNint(2)) = Vint;	% Vint no clindro interno
	V(Nrho,:) = 0; V(:,Nz) = 0;				% 0 no "infinito"
	for i=1:maxiter
		% Aplica as condi��es de contorno de von Neumann
		V(Nrho-1,:) = 0; V(:,Nz-1) = 0;		% "infinito"
		V(1,:) = V(2,:) ; 					% simetria radial
		% Calcula o Laplaciano e o erro
		e = calc .* (0.5 * laplacianoc(V,rdelta)/(1 + rdelta^2));
		erro = max(max(abs(e)));
		if erro <= prec
			return
		end
		% Se o erro aumentou, para, pois o c�lculo n�o vai convergir.
		if erro > last
			return
		end
		% Corrige os valores estimados e continua
		V = V + w * e;
		last = erro;
	end
end
