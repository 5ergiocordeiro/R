% Inicializa��o dos par�metros
% Par�metros
Nrho = Nz = 4000;	% divis�es da grade
prec = 0.002;		% precis�o percentual a obter
maxiter= 40000;		% m�ximo de itera��es tolerado
rlim = zlim = 24;	% posi��o do "infinito"
% Geometria dos cilindros (em mm)
Rext = 20;
Rint = 12;
Lext = 20;
Lint = 16;
% Outros dados do problema
Vext = 0;
Vint = 100;			% V
wmax = 1;
wmin = 0.1;
% Solu��o
V = probcil(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,prec,maxiter,Nrho,Nz,wmax,wmin);
