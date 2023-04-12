%NEWUOA-Trans-Matlab-Version-3.0 
%Copyright: Pengcheng Xie 
%Connect: xpc@lsec.cc.ac.cn

function [Epsilonk] = mlaplace(beta)
  % implicit none
  % double precision::u1, u2, beta, Epsilonk

  u1 = rand(1);
  u2 = rand(1);
  if (u1 >= 0.5e0)
    Epsilonk = -beta * log(1.0e0 - u2);
  else
    Epsilonk = beta * log(u2);
  end
  %    if (ABS(Epsilonk) .GT. 0.5D0)THEN
  %            Epsilonk=0.0
  %    end
end
