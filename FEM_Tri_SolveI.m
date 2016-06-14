function [map, XY] = FEM_Tri_SolveI(n, l)
% Inicializa os dados para um quadrado de lado 'l' e 'n' divisões por dimensão 
	% Cálculo das dimensões da grade
	nd = (n + 1) ^ 2;			% número de nós
	ne = 2 * n ^ 2;				% número de elementos
	np = 4 * n;					% número de nós fixos
	h = l / n;					% largura da grade retangular
	% Inicialização
	global de
	de = @(i,j) (j - 1) * (n + 1) + i;
	map = zeros(ne, 3);			% ... mapa dos nós pertencentes a cada elemento
	XY = zeros(nd, 3);			% ... coordenadas e tipos dos nós
	for i = 1:n+1
		for j = 1:n+1
			% ... para cada nó
			% ...... obtém os nós que limitam o quadrado adjacente
			ei = de(i, j);			% esquerda inferior
			es = de(i, j+1);		% esquerda superior
			di = de(i+1, j);		% direita inferior
			ds = de(i+1, j+1);		% direita superior
			if ((i <= n) && (j <= n))
				% ...... associa os nós aos elementos triangulares contidos no quadrado
				% ......... (sentido anti-horário)
				e = (j - 1) * 2 * n + 2 * (i - 1) + 1;
				map(e, :) = [ei, ds, es];
				map(e + 1, :) = [ei, di, ds];
			end
			% ...... verifica se o nó é fixo
			fixo = (j == 1) || (j == n+1) || (i == 1) || (i == n+1);
			id = ei * fixo;
			% ...... associa as coordenadas ao nó
			XY(ei, :) = [[i-1 j-1] * h, id];
		end
	end
end
