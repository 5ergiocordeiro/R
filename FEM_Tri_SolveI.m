function [map, XY] = FEM_Tri_SolveI(n, l)
% Inicializa os dados para um quadrado de lado 'l' e 'n' divis�es por dimens�o 
	% C�lculo das dimens�es da grade
	nd = (n + 1) ^ 2;			% n�mero de n�s
	ne = 2 * n ^ 2;				% n�mero de elementos
	np = 4 * n;					% n�mero de n�s fixos
	h = l / n;					% largura da grade retangular
	% Inicializa��o
	global de
	de = @(i,j) (j - 1) * (n + 1) + i;
	map = zeros(ne, 3);			% ... mapa dos n�s pertencentes a cada elemento
	XY = zeros(nd, 3);			% ... coordenadas e tipos dos n�s
	for i = 1:n+1
		for j = 1:n+1
			% ... para cada n�
			% ...... obt�m os n�s que limitam o quadrado adjacente
			ei = de(i, j);			% esquerda inferior
			es = de(i, j+1);		% esquerda superior
			di = de(i+1, j);		% direita inferior
			ds = de(i+1, j+1);		% direita superior
			if ((i <= n) && (j <= n))
				% ...... associa os n�s aos elementos triangulares contidos no quadrado
				% ......... (sentido anti-hor�rio)
				e = (j - 1) * 2 * n + 2 * (i - 1) + 1;
				map(e, :) = [ei, ds, es];
				map(e + 1, :) = [ei, di, ds];
			end
			% ...... verifica se o n� � fixo
			fixo = (j == 1) || (j == n+1) || (i == 1) || (i == n+1);
			id = ei * fixo;
			% ...... associa as coordenadas ao n�
			XY(ei, :) = [[i-1 j-1] * h, id];
		end
	end
end
