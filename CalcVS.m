function V = CalcVS(nx, ny, jd, iw, epsilonr)
% Calcula o potencial elétrico, por meio da solução de um sistema linear, na região limitada por 0 <= i <= 'nx' e 0 <= j <= 'ny', com dielátrico de permissividade 'epsilonr' na posição (0 <= i <= 'nx', 0 <= j <= 'jd') e placa condutora na posição (0 <= i <= 'iw', j = 'jd').
	% Criação dos mapas direto e inverso dos pontos
	nA = (nx - 1) * (ny - 1) - iw;
	nb = nx + ny + iw - 1;
	map = zeros(nA + nb, 3);
	imap = zeros(nx, ny);
	idx = [1, 1 + nA];
	for i = 1:nx
		for j = 1:ny
			if (j == jd) && (i <= iw)
				% placa condutora
				parms = [2 1];
			elseif (i == nx) || (j == ny)
				% entorno
				parms = [2 0];
			elseif (j == jd)
				% interface entre o substrato e o ar
				parms = [1 2];
			elseif (i == 1) || (j == 1)
				% eixos X e Y
				parms = [1 3];
			else
				% ar ou substrato
				parms = [1 4];
			end
			point = idx(parms(1));
			idx(parms(1)) = idx(parms(1)) + 1;
			map(point,:) = [i, j, parms(2)];
			imap(i,j) = point;
		end
	end
	% Montagem da matriz de coeficientes e do vetor de contribuição dos pontos fixos
	global epsilon0;
	epsilon = [1 epsilonr];
	p = epsilon / (2 * sum(epsilon));
	A = zeros(nA, nA);
	b = zeros(nA,1);
	for point = 1:nA
		parms = num2cell(map(point,:));
		[i, j, type] = deal(parms{:});
		A(point,point) = 4;
		near = [0, imap(i+1,j), 0, imap(i,j+1)];
		if (i > 1)
			near(1) = imap(i-1,j);
		else
			near(1) = near(2);
		end
		if (j > 1)
			near(3) = imap(i,j-1);
		else
			near(3) = near(4);
		end
		for l = 1:4
			pos = near(l);
			if pos <= nA
				if type == 2 && l == 3
					factor = - 4 * p(2);
				elseif type == 2 && l == 4
					factor = - 4 * p(1);
				else
					factor = - 1;
				end
				A(point,pos) = A(point,pos) + factor;
			else
				b(point) = b(point) + map(pos,3);
			end
		end
	end
	% Resolve o sistema
	x = A \ b;
	% Combina os pontos fixos e os livres
	V = zeros(nx,ny);
	for i = 1:nx
		for j = 1:ny
			point = imap(i,j);
			type = map(point, 3);
			if type == 1
				V(i,j) = 1;
			elseif point <= nA
				V(i,j) = x(point);
			end
		end
	end
	% Calcula o Lapalaciano, para testar a solução
	%{
	lap = zeros(anx-2,any-2);
	for i = 2:anx-1
		for j = 2:any-1
			lap(i,j) = V(i,j) - 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1));
		end
	end
	imagesc(lap), colorbar;
	%}
end
