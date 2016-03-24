function L = probtrack(N, d, W, a)
%{
Retorna a indutância do indutor impresso numa placa quadrada com a geometria dada, onde:
	N: número de trilhas
	d: distância entre as trilhas (mm)
	W: largura das trilhas (mm)
	a: tamanho do lado da placa (mm)
Limitações:
	despreza a profundidade dos elementos,
	o tratamento dos cantos das trilhas não está perfeito,
	despreza componentes do fluxo magnético tangenciais à placa,
	a susceptibilidade da fenolite foi considerada igual a 1.
%}
k = 1.0e-7;		% H/m 
dl = 1e-3;		% Elemento diferencial de comprimento (m)
dS = 1e-6;		% Elemento diferencial de área (m2)
% Cria as trilhas
free = ones(a,a);
placa = zeros(a,a);
iniciox = inicioy = 1;
for i = 1:N
	direcao = mod(i,4);
	switch direcao
		case 1
			dist = sum(free(iniciox:a,inicioy),1);
			placa(iniciox:iniciox+dist-1,inicioy:inicioy+W-1) = 101;
			free(iniciox:iniciox+dist-1,inicioy:inicioy+W-1) = 0;
			free(iniciox:iniciox+dist-W+1,inicioy+W-1:inicioy+W+d-1) = 0;
			iniciox = iniciox+dist-1;
			inicioy = inicioy + W;
		case 2
			dist = sum(free(iniciox,inicioy:a),2);
			placa(iniciox-W+1:iniciox,inicioy:inicioy+dist-1) = 100;
			free(iniciox-W+1:iniciox,inicioy:inicioy+dist-1) = 0;
			free(iniciox-W-d+1:iniciox-W+1,inicioy:inicioy+dist-W+1) = 0;
			iniciox = iniciox - W;
			inicioy = inicioy + dist - 1;
		case 3
			dist = sum(free(1:iniciox,inicioy),1);
			placa(iniciox-dist+1:iniciox,inicioy-W+1:inicioy) = 101;
			free(iniciox-dist+1:iniciox,inicioy-W+1:inicioy) = 0;
			free(iniciox-dist+W-1:iniciox,inicioy-W-d+1:inicioy-W+1) = 0;
			iniciox = iniciox - dist +1;
			inicioy = inicioy - W;
		case 0
			dist = sum(free(iniciox,1:inicioy),2);
			placa(iniciox:iniciox+W-1,inicioy-dist+1:inicioy) = 100;
			free(iniciox:iniciox+W-1,inicioy-dist+1:inicioy) = 0;
			free(iniciox+W-1:iniciox+W+d-1,inicioy-dist+W+1:inicioy) = 0;
			iniciox = iniciox + W;
			inicioy = inicioy - dist + 1;
	end
end
% Mostra as trilhas
pcolor(placa);
shading interp;
% Calcula a densidade de fluxo magnético em cada ponto {i,k} da placa
% considera a circulação de uma corrente de 1 A em cada elemento {ii,kk}
B = zeros(a,a);
for i = 1:a
	for k = 1:a
		sB = 0;
		for ii = 1:a
			for kk = 1:a
				direcao = placa(ii,kk) - 100;
				if direcao >= 0
					% Circula corrente neste ponto
					R = [i k] - [ii kk];
					dist = sqrt(sum(R .^ 2));
					if dist == 0 
						dB = 0;
					else
						if direcao == 1 
							dB = R(2) / dist^3;
						else
							dB = R(1) / dist^3;
						end
					end
				sB = sB + dB;
				end
			end
		end
		B(i,k) = sB;
	end
end
% Mostra a densidade do campo
pcolor(B);
shading interp;
% Calcula a indutância
L = k * sum(sum(B)) * dS * dl;
disp(sprintf("N = %d, W = %d mm, d = %d mm, L = %f H",N,W,d,L));
