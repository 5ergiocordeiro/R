function V = FEM_Tri_SolveS(map, XY, Vp)
% Resolve o problema pelo m�todo dos elementos finitos
	[ne, nd, np] = deal(size(map,1), size(XY,1), size(Vp,1)); 
	b = zeros(nd, 1);
	C = zeros(nd, nd);
	for i = 1:ne
		% ... para cada elemento
		% ...... obt�m os n�s que o limitam e suas coordenadas
		d = map(i, :);
		xy = XY(d, :);
		% ...... calcula a matriz de coeficientes locais 
		p = xy([2 3 1],2) - xy([3 1 2],2);
		q = xy([3 1 2],1) - xy([2 3 1],1);
		area = 2 * abs(p(2) * q(3) - q(2) * p(3));
		ce = (p * p' + q * q') / area;
		% ....... para cada n�
		for j = 1:3
			% ......... verifica se � n� fixo
			ir = d(j);
			p = xy(j,3);
			if (p > 0)
				% ............ n� fixo; potencial constante
				C(ir, ir) = 1;
				b(ir) = Vp(p);
			else
				% ............ n� livre
				for k = 1:3
					% ............... calcula a contribui��o dos n�s do elemento
					ic = d(k);
					p = xy(k,3);
					if (p > 0)
						% .................. constribui��o de n� fixo
						b(ir) = b(ir) - ce(j,k) * Vp(p);
					else
						% .................. constribui��o de n� livre
						C(ir, ic) = C(ir, ic) + ce(j, k);
					end
				end
			end
		end	
	end
	% Resolve o sistema
	V = C \ b;
	end
