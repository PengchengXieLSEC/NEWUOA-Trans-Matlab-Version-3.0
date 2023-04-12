%NEWUOA-Trans-Matlab-Version-3.0 
%Copyright: Pengcheng Xie 
%Connect: xpc@lsec.cc.ac.cn

nvar = 20;
flag = 1;
FLAG_Epsilonk = 1;
X = rand(1, nvar);
Y = rand(1, nvar);
probs = textread('problems', '%s');
fid = fopen('solution.csv', 'wt');
for IP = 1:1
  prob = probs(IP);

  for I = 1:nvar
    X(I) = 1.0e0;
  end
  % print*,PROBLEM,IP
  % - NPT is the number of interpolation conditions. Its value must be in the
  %      interval [N+2,(N+1)(N+2)/2] with N=size(X).
  %      NPT=2*N+1 is usually appropriate
  NPT = 2 * nvar + 1;
  CALFUN = @(x)(testfun(prob, x));
  CALGREY = @(g, k, f, e, e1, e2, rd) ...
    (calgrey(g, nvar, k, f, e, e1, e2, rd, FLAG_Epsilonk));
  [X, NF, ~] = newuoa(CALFUN, CALGREY, nvar, NPT, X, 1.0e0, 1.0e-6, 0, 10000, flag);

  disp(['PROBLEM: ', prob]);
  %          for I = 1:nvar
  %              X(I)=X(I)/XSQ
  %          end

  %20230413 F = evalfun(prob, X);
  F = evalfun(prob, X,nvar);
  %          F=(X(1)-1)^4+100*(X(1)-1)^2
  G = 0.0e0;
  for I = 1:nvar
    G = G + (X(I)) ^ 2;
  end
  G = (G - 1) ^ 2;
  % G=X(1)^2+X(2)^2
  % F=X(1)^4+X(2)^4
  fprintf(fid, '%.0f', NF);
  fprintf(fid, '\t%.8f', F);
  for I = 1:nvar
    fprintf(fid, '\t%.8f', X(I));
  end
  fprintf(fid, '\r\n');

end
fclose(fid);
