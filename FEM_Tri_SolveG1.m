function [Vmap] = FEM_Tri_SolveG1(n, l)
% Resolve a equação de Laplace pelo método dos elementos finitos em um quadrado de lado 'l', usando 'n' elementos triangulares.
% Os elementos são triângulos retângulos e formam uma grade retangular.
% Problema 15.29 do livro “Elementos de eletromagnetismo", de Matthew Sadiku
% Retorna um mapa do potencial em cada ponto da grade retangular.
	% Inicializa os dados para o método, de acordo com a geometria dada
	[map, XY] = FEM_Tri_SolveI(n, l);
	% Inicializa os pontos fixos
	nd = size(XY,1);
	Vp = zeros(nd, 1);			% V = 0 na base e nos lados
	Vp((nd-n):nd) = 100;		% V = 100 no topo
	% Calcula o potencial em cada nó
	V = FEM_Tri_SolveS(map, XY, Vp);
	% Exibe os resultados
	ref = @(i,j,h) FEM_Tri_SolveR1(i, j, h);
	Vmap = FEM_Tri_SolveP(V, ref, n, l/n);
end

function Vref = FEM_Tri_SolveR1(i, j, h)
% Calcula a distribuição ideal
	pih = pi * h;
	infty = 100;
	sum = 0;
	for k = 0:infty
		m = 2 * k + 1;
		sum = sum + sin(m * pih * (i-1)) * sinh(m * pih * (j-1)) / (m * sinh(m * pi));
	end
	Vref = 4 * 100 / pi * sum;
end