function e=laplacianoc(V,rdelta)
% Calcula o Laplaciano de 'V' em coordenadas cilíndricas. 'rdelta' é a razão entre Delta rho e Delta z escolhida.
% Os melhores resultados são obtidos quando rdelta = 1.
% Chamada pela função 'calcv'.
	[nrow,ncol] = size(V);
	e = zeros(nrow,ncol);
	m = 0.5 * [2:nrow-1].^(-1);
	% e(i,k) = V(i+1,k)*(1 + 1/(2*i)) + V(i-1,k)*(1-1/(2i)) + ((V(i,k+1) + V(i,k-1))*(rdelta^2) - V(i,k)*2*(1+rdelta^2);
	e(2:nrow-1,2:ncol-1) = V(3:nrow,2:ncol-1) .* (1+m) + V(1:nrow-2,2:ncol-1) .* (1-m) + (V(2:nrow-1,3:ncol) + V(2:nrow-1,1:ncol-2))*(rdelta^2) - V(2:nrow-1,2:ncol-1)*2*(1+rdelta^2);	
end
	