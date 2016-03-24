function L = probtrack2(ne, d, W, a)
%{
Retorna a indutância de um indutor plano composto por espiras quadradas concêntricas:
	ne: número de espiras
	d: distância entre as espiras (mm)
	W: largura das trilhas (mm)
	a: tamanho do lado da placa (mm)
Limitações:
	o tratamento dos cantos das trilhas não está perfeito,
	despreza componentes do fluxo magnético tangenciais à placa,
	a susceptibilidade da fenolite foi considerada igual a 1.	
%}
N = 4 * ne;		% Número de trilhas
k = 1.0e-7;		% H/m
% Cria as trilhas
placa = zeros(a,a);
iniciox = inicioy = 1;
dist = a;
for i = 1:ne
	for direcao = 0:3;
		switch direcao
			case 1
				placa(iniciox:iniciox+dist-1,inicioy:inicioy+W-1) = 101;
				placa(iniciox+W:iniciox+W+d-1,inicioy:inicioy+W-1) = 0;
				placa(iniciox+W+d:iniciox+2*W+d-1,inicioy+W-1:inicioy+2*W+d-1) = 99;				
			case 2
				placa(a-iniciox-W+2:a-iniciox+1,inicioy:inicioy+dist-1) = 100;
			case 3
				placa(iniciox:iniciox+dist-1,a-inicioy-W+1:a-inicioy+1) = 101;
			case 0
				placa(iniciox:iniciox+W-1,inicioy:inicioy+dist-1) = 100;
		end
	end
	dist = dist - 2 * ( W + d);
	iniciox = iniciox + W + d;
	inicioy = inicioy + W + d;
end
% Mostra as trilhas
pcolor(placa);
shading interp;
