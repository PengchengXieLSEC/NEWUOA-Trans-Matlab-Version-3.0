%NEWUOA-Trans-Matlab-Version-3.0 
%Copyright: Pengcheng Xie 
%Connect: xpc@lsec.cc.ac.cn

function [G, FLAG, EPS, EPS1, EPS2, RAND] = calgrey(G, N, K, FLAG, EPS, EPS1, EPS2, RAND, FLAG_Epsilonk)
  % implicit real*8 (a-h,o-z)
  % real*8 :: X(:),G,PI,PP,EPS1(:),EPS2(:),E,RAND(:),FF
  % integer*4 :: K,I,N,FLAG,INFO
  % character(len=15) :: PROBLEM
  % double precision :: Epsilonk

  if (FLAG_Epsilonk == 0)
    G = G + 0;
  elseif (FLAG_Epsilonk == 1)

    if (K > FLAG)
      FLAG = K;
      Epsilonk = mlaplace(100.0e0 / K);
      EPS(K) = Epsilonk;
    else
      Epsilonk = EPS(K);
    end

    G = G + EPS(K);
  elseif (FLAG_Epsilonk == 2)

    if (K > FLAG)
      Epsilonk = mlaplace(100.0e0 / K);
      EPS1(K) = Epsilonk * 10.0e0;

      Epsilonk = rand(1);
      Epsilonk = (Epsilonk -0.5e0) * 2 * (1.0e0 / K);
      EPS2(K) = Epsilonk;

      RAND(K) = rand(1);

      FLAG = K;
    end

    if (RAND(K) > 0.5)
      G = G + EPS1(K);
    else
      if (K > 2 * N + 1)
        Epsilonk = EPS2(K);
        G = G + EPS2(K) * (G);
      else
        G = G + EPS1(K);
      end
    end
  end
end
