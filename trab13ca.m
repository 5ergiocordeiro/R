% Inicialização dos parâmetros
% Parâmetros
Nrho = Nz = 4000;	% divisões da grade
prec = 0.002;		% precisão percentual a obter
maxiter= 40000;		% máximo de iterações tolerado
rlim = zlim = 24;	% posição do "infinito"
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
% Solução
V = probcil(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,prec,maxiter,Nrho,Nz,wmax,wmin);
