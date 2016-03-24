function e=laplacianok(V)
% Calcula o Laplaciano de 'V' por meio de uma grade perfeitamente quadrada ("em coordenadas cartesianas").
% Chamada pela função 'calcv'.
	[nrow,ncol] = size(V);
	e = zeros(nrow,ncol);
	% e(i,k) = V(i+1,k) + V(i-1,k) + V(i,k+1) + V(i,k-1) - 4 * V(i,k)
	e(2:nrow-1,2:ncol-1) = V(3:nrow,2:ncol-1) + V(1:nrow-2,2:ncol-1) + V(2:nrow-1,3:ncol) + V(2:nrow-1,1:ncol-2) - 4*V(2:nrow-1,2:ncol-1);
end
	