function [Vmap] = FEM_Tri_SolveG2(n, l)
% Resolve a equa��o de Laplace pelo m�todo dos elementos finitos em um quadrado de lado 'l', usando 'n' elementos triangulares.
% Os elementos s�o tri�ngulos ret�ngulos e formam uma grade retangular.
% Problema 15.30 do livro "", de Matthew Sadiku
% Retorna um mapa do potencial em cada ponto da grade retangular.
	% Inicializa os dados para o m�todo, de acordo com a geometria dada
	[map, XY] = FEM_Tri_SolveI(n, l);
	% Inicializa os pontos fixos
	nd = size(XY,1);
	h = l / n;
	Vp = zeros(nd, 1);							% V = 0 na base e nos lados
	Vp((nd-n):nd) = 100 * sin(pi * h * (0:n));	% V = 100 sin(pi x) no topo
	% Calcula o potencial em cada n�
	V = FEM_Tri_SolveS(map, XY, Vp);
	% Exibe os resultados
	ref = @(i,j,h) FEM_Tri_SolveR2(i, j, h);
	Vmap = FEM_Tri_SolveP(V, ref, n, h);
end

function Vref = FEM_Tri_SolveR2(i, j, h)
% Calcula a distribui��o ideal
	pih = pi * h;
	Vref = 100 * sin(pih * (i-1)) * sinh(pih * (j-1)) / sinh(pi);
end