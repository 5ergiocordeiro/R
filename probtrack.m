function L = probtrack(ne, d, W, a)
%{
Retorna a indut�ncia de um indutor plano composto por espiras quadradas conc�ntricas:
	ne: n�mero de espiras
	d: dist�ncia entre as espiras (mm)
	W: largura das trilhas (mm)
	a: tamanho do lado da placa (mm)
Limita��es:
	o tratamento dos cantos das trilhas n�o est� perfeito,
	despreza componentes do fluxo magn�tico tangenciais � placa,
	as susceptibilidades da fenolite e do cobre foram consideradas iguais a 1.	
%}
% Constantes
kmu = 1e-4;		% uH/mm
% Cria as trilhas
placa = zeros(a+1,a+1);
cond = zeros(4 * ne, 3);
for i = 1:ne
	% Calcula as caracter�sticas da espira
	lado = a - 2 * (i - 1) * (W + d);
	gapx = (i + 1) * (W + d) - 1;
	inicio = 1 + (i - 1) * (W + d);
	% Desenha a espira
	placa(inicio:inicio+lado-1,inicio:inicio+W-1) = 101;
	placa(gapx,inicio:inicio+W-1) = 0;
	if ne > i
		placa(gapx+1:gapx+W+d-1,inicio+W:inicio+W+d-1) = 99;
	end
	placa(a-inicio-W+2:a-inicio+1,inicio:inicio+lado-1) = 100;	
	placa(inicio:inicio+lado-1,a-inicio-W+2:a-inicio+1) = 103;
	placa(inicio:inicio+W-1,inicio:inicio+lado-1) = 102;
	% Armazena as caracter�sticas de cada condutor
	%(lado, dire��o, posi��o)
	ncond = (i - 1) * 4 + 1;
	cond(ncond:ncond+3,1) = lado;
	cond(ncond:ncond+3,2) = [101 100 103 102];
	cond(ncond:ncond+3,3) = [inicio+W/2, a-inicio+W/2, a-inicio+W/2, inicio+W/2];
end
% Mostra as trilhas
% pcolor(placa);
% shading flat;
% Calcula a indut�ncia
% Primeiro m�todo:
% Soma a indut�ncia correspondente a cada condutor
B1 = 0;
for i = 1:(4 * ne)
	l1 = cond(i,1);
	d1 = cond(i,2);
	p1 = cond(i,3);
	% Indut�ncias m�tuas com os condutores restantes
	for j = (i+1):(4*ne)
		l2 = cond(j,1);
		d2 = cond(j,2);
		p2 = cond(j,3);
		if (mod(d1,2) != mod(d2,2)) || (d1 == d2)
			continue
		end
		dist = abs(p2 - p1);
		l = min(l2, l1);
		dLp = 4 * l * log(2 * dist / W - 1);
		B1 = B1 + dLp;
	end
	% Indut�ncia pr�pria
	dLp = l1 * 0.5;
	B1 = B1 + dLp;
end
% Segundo m�todo:
% Calcula a densidade de fluxo magn�tico em cada ponto {i,k} da placa
% considera a circula��o de uma corrente de 1 A em cada elemento {ii,kk}
% assim, em cada condutor circular� uma corrente proporcional a W.
B2 = 0;
for i = 1:a
	for k = 1:a
		for ii = 1:a
			for kk = 1:a
				direcao = placa(ii,kk) - 100;
				if direcao < 0
					continue
				end
				% Circula corrente neste ponto
				R = [i k] - [ii kk];
				dist2 = sum(R .* R);
				if dist2 == 0
					continue
				end
				switch direcao
					case 0
						dB = - R(1) / dist2^1.5;
					case 1
						dB = R(2) / dist2^1.5;
					case 2
						dB = R(1) / dist2^1.5;
					case 3
						dB = - R(2) / dist2^1.5;
				end
				% Considera apenas o fluxo concatenado, que ser� sempre positivo nas condi��es do problema
				if (dB > 0)
					B2 = B2 + dB;
				end
			end
		end
	end
end
% Terceiro m�todo:
% Soma a indut�ncia correspondente a cada espira
B3 = 0;
for i = 1:ne
	l = cond(i,1);
	dist = cond(i,2);
	% Indut�ncia pr�pria da espira
	dB = l * ( 2  + 8 * log(2 * dist / W - 1));
	B3 = B3 + dB;
end
% Calcula a indut�ncia
L2 = kmu * B2 / W;
L1 = kmu * B1;
L3 = kmu * B3;
disp(sprintf("Ne = %d, W = %d mm, d = %d mm, L1 = %f uH, L2 = %f uH, L3 = %f uH", ne, W, d, L1, L2, L3));
L = L1;
