function [map, XY, Vp] = FEM_Tri_SolveI(n, l)
	% Cálculo das dimensões da grade
	nd = (n + 1) ^ 2;			% número de nós
	ne = 2 * n ^ 2;				% número de elementos
	np = 4 * n;					% número de nós fixos
	h = l / n;					% largura da grade retangular
	% Inicialização
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
	Vp = zeros(nd, 1);			% potencial nos nós fixos
	Vp((nd-n):nd) = 1;
end

function V = FEM_Tri_SolveS(map, XY, Vp)
% Resolve o problema pelo método dos elementos finitos
	[ne, nd, np] = deal(size(map,1), size(XY,1), size(Vp,1)); 
	b = zeros(nd, 1);
	C = zeros(nd, nd);
	for i = 1:ne
		% ... para cada elemento
		% ...... obtém os nós que o limitam e suas coordenadas
		d = map(i, :);
		xy = XY(d, :);
		% ...... calcula a matriz de coeficientes locais 
		p = xy([2 3 1],2) - xy([3 1 2],2);
		q = xy([3 1 2],1) - xy([2 3 1],1);
		area = 2 * abs(p(2) * q(3) - q(2) * p(3));
		ce = (p * p' + q * q') / area;
		% ....... para cada nó
		for j = 1:3
			% ......... verifica se é nó fixo
			ir = d(j);
			p = xy(j,3);
			if (p > 0)
				% ............ nó fixo; potencial constante
				C(ir, ir) = 1;
				b(ir) = Vp(p);
			else
				% ............ nó livre
				for k = 1:3
					% ............... calcula a contribuição dos nós do elemento
					ic = d(k);
					p = xy(k,3);
					if (p > 0)
						% .................. constribuição de nó fixo
						b(ir) = b(ir) - ce(j,k) * Vp(p);
					else
						% .................. constribuição de nó livre
						C(ir, ic) = C(ir, ic) + ce(j, k);
					end
				end
			end
		end	
	end
	% Resolve o sistema
	V = C \ b;
	end
	
function V = FEM_Tri_SolveP(V, n)
	% Organiza os potenciais em uma matriz
	Vmap = zeros(n+1,n+1);
	for i = 1:n+1
		for j = 1:n+1
			ei = de(i, j);
			Vmap(i,j) = V(ei);
		end
	end
	nd = (n + 1) ^ 2;
	ne = 2 * n ^ 2;
	np = 4 * n;
	% Grava o resultado
	diary FEM_Tri_SolveS.out
	[nd, ne, np]
	Vmap'
	diary off
end