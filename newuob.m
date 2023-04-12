%NEWUOA-Trans-Matlab-Version-3.0 
%Copyright: Pengcheng Xie 
%Connect: xpc@lsec.cc.ac.cn

function [NF, X, F] = newuob(Func, Gunc, N, NPT, X, RHOBEG, RHOEND, IPRINT, MAXFUN, FLAG)

  % implicit none
  % integer, intent(in) :: N, NPT, IPRINT,MAXFUN, NDIM
  % real(kind=8), intent(inout), dimension(:) :: X
  % real(kind=8), intent(in) :: RHOBEG, RHOEND
  % real(kind=8), intent(inout) :: XBASE(*), XOPT(*),XNEW(*),  &
  %    FVAL(*), GQ(*), HQ(*), PQ(*), D(*), VLAG(*), W(*)
  % real(kind=8), intent(inout) :: BMAT(NDIM,*),ZMAT(NPT,*), XPT(NPT,*)
  % integer, intent(inout) :: NF

  % interface
  %    subroutine CALFUN(i_x, o_f,i_k,i_problem)
  %       real(kind=8), dimension(:) :: i_x
  %       real(kind=8) :: o_f
  %       integer(kind=4) :: i_k
  %       character(len=15) :: i_problem
  %    end subroutine
  % end interface

  % local variables
  % integer :: NP, NH,NPTM,NFTEST, NFM, NFMM, ITEMP, JPT, IPT, &
  %    IH, IDZ, ITEST, NFSAV, KNEW, I, IP, J, JP, K, KSAVE, &
  %    KTEMP, KOPT
  % real(kind=8) :: HALF, ONE, TENTH, ZERO, RHOSQ, RECIP,RECIQ, XIPT, &
  %    XJPT, FBEG, FOPT, RHO, DELTA, DIFFA, DIFFB, &
  %    XOPTSQ, DSQ, DNORM, RATIO, TEMP, CRVMIN, TEMPQ, &
  %    ALPHA, BETA, BSUM, DETRAT, DIFF, DIFFC, DISTSQ, &
  %    DX, F, FSAVE, GISQ, GQSQ, HDIAG, SUM, SUMA, SUMB, &
  %    DSTEP, SUMZ, VQUAD
  % character(len=15) :: PROBLEM
  NDIM = NPT + N;
  ZERO = 0.0e0;
  [NFM, NFMM, JPT, IPT, IDZ, ITEST, NFSAV, KNEW, KSAVE, KOPT] = deal(0);
  [ZERO, XIPT, XJPT, FBEG, FOPT, RHO, DELTA, DIFFA, DIFFB, ...
     XOPTSQ, DSQ, DNORM, RATIO, CRVMIN, ALPHA, BETA, ...
     DIFF, DIFFC, F, FSAVE, DSTEP, VQUAD] = deal(ZERO);
  [XBASE, XNEW, XPT, FVAL, GQ, HQ, PQ, BMAT, ZMAT, D] = deal(ZERO);
  XOPT = zeros(1, N);
  W = zeros(1, 4 * N);
  VLAG = zeros(1, NPT + N);
  TEMP1 = zeros(1, NPT);
  XBEG = zeros(1, N);
  XTEMP = zeros(1, N);
  EPS = zeros(1, MAXFUN);
  EPS1 = zeros(1, MAXFUN);
  EPS2 = zeros(1, MAXFUN);
  RAND = zeros(1, MAXFUN);

  %     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
  %       to the corresponding arguments in SUBROUTINE NEWUOA.
  %     XBASE will hold a shift of origin that should reduce the contributions
  %       from rounding errors to values of the model and Lagrange functions.
  %     XOPT will be set to the displacement from XBASE of the vector of
  %       variables that provides the least calculated F so far.
  %     XNEW will be set to the displacement from XBASE of the vector of
  %       variables for the current calculation of F.
  %     XPT will contain the interpolation point coordinates relative to XBASE.
  %     FVAL will hold the values of F at the interpolation points.
  %     GQ will hold the gradient of the quadratic model at XBASE.
  %     HQ will hold the explicit second derivatives of the quadratic model.
  %     PQ will contain the parameters of the implicit second derivatives of
  %       the quadratic model.
  %     BMAT will hold the last N columns of H.
  %     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
  %       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
  %       the elements of DZ are plus or minus one, as specified by IDZ.
  %     NDIM is the first dimension of BMAT and has the value NPT+N.
  %     D is reserved for trial steps from XOPT.
  %     VLAG will contain the values of the Lagrange functions at a new point X.
  %       They are part of a product that requires VLAG to be of length NDIM.
  %     The array W will be used for working space. Its length must be at least
  %       10*NDIM = 10*(NPT+N).

  %     Set some constants.

  HALF = 0.5e0;
  ONE = 1.0e0;
  TENTH = 0.1e0;
  NP = N + 1;
  NH = (N * NP) / 2;
  NPTM = NPT - NP;
  NFTEST = max(MAXFUN, 1);

  % Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

  for J = 1:N
    XBASE(J) = X(J);
    for K = 1:NPT
      XPT(K, J) = ZERO;
    end
    for I = 1:NDIM
      BMAT(I, J) = ZERO;
    end
  end

  for IH = 1:NH
    HQ(IH) = ZERO;
  end

  for K = 1:NPT
    PQ(K) = ZERO;
    for J = 1:NPTM
      ZMAT(K, J) = ZERO;
    end
  end

  %     Begin the initialization procedure. NF becomes one more than the number
  %     of function values so far. The coordinates of the displacement of the
  %     next initial interpolation point from XBASE are set in XPT(NF,.).

  RHOSQ = RHOBEG * RHOBEG;
  RECIP = ONE / RHOSQ;
  RECIQ = sqrt(HALF) / RHOSQ;
  NF = 0;
  flag = 50;

  while (flag ~= 1111)
    switch flag
      case 50
        NFM = NF;
        NFMM = NF - N;
        NF = NF + 1;
        if (NFM <= 2 * N)
          if (NFM >= 1 && NFM <= N)
            XPT(NF, NFM) = RHOBEG;
          elseif (NFM > N)
            XPT(NF, NFMM) = -RHOBEG;
          end
        else
          ITEMP = (NFMM - 1) / N;
          JPT = NFM - ITEMP * N - N;
          IPT = JPT + ITEMP;
          if (IPT > N)
            ITEMP = JPT;
            JPT = IPT - N;
            IPT = ITEMP;
          end
          XIPT = RHOBEG;
          if (FVAL(IPT + NP) < FVAL(IPT + 1))
            XIPT = -XIPT;
          end
          XJPT = RHOBEG;
          if (FVAL(JPT + NP) < FVAL(JPT + 1))
            XJPT = -XJPT;
          end
          XPT(NF, IPT) = XIPT;
          XPT(NF, JPT) = XJPT;
        end

        %     Calculate the next value of F, label 70 being reached immediately
        %     after this calculation. The least function value so far and its index
        %     are required.

        for J = 1:N
          X(J) = XPT(NF, J) + XBASE(J);
        end
        flag = 310;

      case 70
        FVAL(NF) = F;
        if (NF == 1)
          for I = 1:N
            XBEG(I) = X(I);
          end
          FBEG = F;
          FOPT = F;
          KOPT = 1;
        elseif (F < FOPT)
          FOPT = F;
          KOPT = NF;
        end
        %     Set the nonzero initial elements of BMAT and the quadratic model in
        %     the cases when NF is at most 2*N+1.

        if (NFM <= 2 * N)
          if (NFM >= 1 && NFM <= N)
            GQ(NFM) = (F - FBEG) / RHOBEG;
            if (NPT < NF + N)
              BMAT(1, NFM) = -ONE / RHOBEG;
              BMAT(NF, NFM) = ONE / RHOBEG;
              BMAT(NPT + NFM, NFM) = -HALF * RHOSQ;
            end
          elseif (NFM > N)
            BMAT(NF - N, NFMM) = HALF / RHOBEG;
            BMAT(NF, NFMM) = -HALF / RHOBEG;
            ZMAT(1, NFMM) = -RECIQ - RECIQ;
            ZMAT(NF - N, NFMM) = RECIQ;
            ZMAT(NF, NFMM) = RECIQ;
            IH = (NFMM * (NFMM + 1)) / 2;
            TEMP = (FBEG - F) / RHOBEG;
            HQ(IH) = (GQ(NFMM) - TEMP) / RHOBEG;
            GQ(NFMM) = HALF * (GQ(NFMM) + TEMP);
          end
          %     Set the off-diagonal second derivatives of the Lagrange functions and
          %     the initial quadratic model.
        else
          IH = (IPT * (IPT - 1)) / 2 + JPT;
          if (XIPT < ZERO)
            IPT = IPT + N;
          end
          if (XJPT < ZERO)
            JPT = JPT + N;
          end
          ZMAT(1, NFMM) = RECIP;
          ZMAT(NF, NFMM) = RECIP;
          ZMAT(IPT + 1, NFMM) = -RECIP;
          ZMAT(JPT + 1, NFMM) = -RECIP;
          HQ(IH) = (FBEG - FVAL(IPT + 1) - FVAL(JPT + 1) + F) / (XIPT * XJPT);
        end

        if (NF < NPT)
          flag = 50;
        else
          %     Begin the iterative procedure, because the initial model is complete.
          RHO = RHOBEG;
          DELTA = RHO;
          IDZ = 1;
          DIFFA = ZERO;
          DIFFB = ZERO;
          ITEST = 0;
          XOPTSQ = ZERO;
          for I = 1:N
            XOPTSQ = XOPTSQ + XOPT(I) ^ 2;
            XOPT(I) = XPT(KOPT, I);
          end
          flag = 90;
        end

      case 90
        NFSAV = NF;
        flag = 100;
        %     Generate the next trust region step and test its length. Set KNEW
        %     to -1 if the purpose of the next F will be to improve the model.

      case 100
        KNEW = 0;
        [STEP_, D_, G_, HD_, HS_, CRVMIN_] = trsapp(N, NPT, XOPT, XPT, ...
          GQ, HQ, PQ, DELTA, D, W(1:N), W(NP:2 * N), W(NP + N:3 * N), ...
          W(NP + 2 * N:4 * N), CRVMIN);
        CRVMIN = CRVMIN_;
        D = STEP_;
        W(NP + 2 * N:4 * N) = HS_;
        W(NP + N:3 * N) = HD_;
        W(NP:2 * N) = G_;
        W(1:N) = D_;
        DSQ = ZERO;
        for I = 1:N
          DSQ = DSQ + D(I) ^ 2;
        end
        DNORM = min(DELTA, sqrt(DSQ));
        if (DNORM < HALF * RHO)
          KNEW = -1;
          DELTA = TENTH * DELTA;
          RATIO = -1.0e0;
          if (DELTA <= 1.5e0 * RHO)
            DELTA = RHO;
          end
          TEMP = 0.125e0 * CRVMIN * RHO * RHO;
          if (NF <= NFSAV + 2 || TEMP <= max([DIFFA, DIFFB, DIFFC]))
            flag = 460;
          else
            flag = 490;
          end
        else
          flag = 120;
        end

        %     Shift XBASE if XOPT may be too far from XBASE. First make the changes
        %     to BMAT that do not depend on ZMAT.
      case 120
        if (DSQ <= 1.0e-3 * XOPTSQ)
          TEMPQ = 0.25e0 * XOPTSQ;
          for K = 1:NPT
            SUM = ZERO;
            for I = 1:N
              SUM = SUM + XPT(K, I) * XOPT(I);
            end
            TEMP = PQ(K) * SUM;
            SUM = SUM - HALF * XOPTSQ;
            W(NPT + K) = SUM;
            for I = 1:N
              GQ(I) = GQ(I) + TEMP * XPT(K, I);
              XPT(K, I) = XPT(K, I) - HALF * XOPT(I);
              VLAG(I) = BMAT(K, I);
              W(I) = SUM * XPT(K, I) + TEMPQ * XOPT(I);
              IP = NPT + I;
              for J = 1:I
                BMAT(IP, J) = BMAT(IP, J) + VLAG(I) * W(J) + W(I) * VLAG(J);
              end
            end
          end

          %      the revisions of BMAT that depend on ZMAT are calculated.

          for K = 1:NPTM
            SUMZ = ZERO;
            for I = 1:NPT
              SUMZ = SUMZ + ZMAT(I, K);
              W(I) = W(NPT + I) * ZMAT(I, K);
            end
            for J = 1:N
              SUM = TEMPQ * SUMZ * XOPT(J);
              for I = 1:NPT
                SUM = SUM + W(I) * XPT(I, J);
              end
              VLAG(J) = SUM;
              if (K < IDZ)
                SUM = -SUM;
              end
              for I = 1:NPT
                BMAT(I, J) = BMAT(I, J) + SUM * ZMAT(I, K);
              end
            end

            for I = 1:N
              IP = I + NPT;
              TEMP = VLAG(I);
              if (K < IDZ)
                TEMP = -TEMP;
              end
              for J = 1:I
                BMAT(IP, J) = BMAT(IP, J) + TEMP * VLAG(J);
              end
            end
          end

          %    The following instructions complete the shift of XBASE, including
          %    the changes to the parameters of the quadratic model.

          IH = 0;
          for J = 1:N
            W(J) = ZERO;
            for K = 1:NPT
              W(J) = W(J) + PQ(K) * XPT(K, J);
              XPT(K, J) = XPT(K, J) - HALF * XOPT(J);
            end

            for I = 1:J
              IH = IH + 1;
              if (I < J)
                GQ(J) = GQ(J) + HQ(IH) * XOPT(I);
              end
              GQ(I) = GQ(I) + HQ(IH) * XOPT(J);
              HQ(IH) = HQ(IH) + W(I) * XOPT(J) + XOPT(I) * W(J);
              BMAT(NPT + I, J) = BMAT(NPT + J, I);
            end
          end

          for J = 1:N
            XBASE(J) = XBASE(J) + XOPT(J);
            XOPT(J) = ZERO;
          end
          XOPTSQ = ZERO;
        end

        %     Pick the model step if KNEW is positive. A different choice of D
        %     may be made later, if the choice of D by biglag causes substantial
        %     cancellation in DENOM.

        if (KNEW > 0)
          [ALPHA, D_, HCOL_, GC_, GD_, S_, W_] = ...
            biglag (N, NPT, XOPT, XPT, BMAT, ZMAT, IDZ, KNEW, DSTEP, ...
            D, VLAG(1:NPT), VLAG(NPT + 1:NPT + N), W(1:N), W(NP:2 * N), ...
            W(NP + N:3 * N));
          W(NP + N:3 * N) = W_;
          W(NP:2 * N) = S_;
          W(1:N) = GD_;
          VLAG(1:NPT) = HCOL_;
          VLAG(NPT + 1:NPT + N) = GC_;
          D = D_;
        end

        %     Calculate VLAG and BETA for the current choice of D. The first NPT
        %     components of W_check will be held in W.

        for K = 1:NPT
          SUMA = ZERO;
          SUMB = ZERO;
          SUM = ZERO;
          for J = 1:N
            SUMA = SUMA + XPT(K, J) * D(J);
            SUMB = SUMB + XPT(K, J) * XOPT(J);
            SUM = SUM + BMAT(K, J) * D(J);
          end
          W(K) = SUMA * (HALF * SUMA + SUMB);
          VLAG(K) = SUM;
        end
        BETA = ZERO;
        for K = 1:NPTM
          SUM = ZERO;
          for I = 1:NPT
            SUM = SUM + ZMAT(I, K) * W(I);
          end

          if (K < IDZ)
            BETA = BETA + SUM * SUM;
            SUM = -SUM;
          else
            BETA = BETA - SUM * SUM;
          end

          for I = 1:NPT
            VLAG(I) = VLAG(I) + SUM * ZMAT(I, K);
          end
        end
        BSUM = ZERO;
        DX = ZERO;
        for J = 1:N
          SUM = ZERO;
          for I = 1:NPT
            SUM = SUM + W(I) * BMAT(I, J);
          end
          BSUM = BSUM + SUM * D(J);
          JP = NPT + J;
          for K = 1:N
            SUM = SUM + BMAT(JP, K) * D(K);
          end
          VLAG(JP) = SUM;
          BSUM = BSUM + SUM * D(J);
          DX = DX + D(J) * XOPT(J);
        end
        BETA = DX * DX + DSQ * (XOPTSQ + DX + DX + HALF * DSQ) + BETA - BSUM;
        VLAG(KOPT) = VLAG(KOPT) + ONE;
        %     if KNEW is positive and if the cancellation in DENOM is unacceptable,
        %      bigden calculates an alternative model step, XNEW being used for
        %     working space.
        if (KNEW > 0)
          TEMP = ONE + ALPHA * BETA / VLAG(KNEW) ^ 2;
          if (abs(TEMP) <= 0.8e0)
            [W_, VLAG_, BETA_, S_] = bigden (N, NPT, ...
              XOPT, XPT, BMAT, ZMAT, IDZ, NDIM, KOPT, KNEW, D, ...
              W(1:NDIM), VLAG, BETA, XNEW);
            W(1:NDIM) = W_;
            VLAG = VLAG_;
            BETA = BETA_;
            XNEW = S_;
          end
        end
        flag = 290;

        %     Calculate the next value of the objective function.

      case 290
        for I = 1:N
          XNEW(I) = XOPT(I) + D(I);
          X(I) = XBASE(I) + XNEW(I);
        end
        NF = NF + 1;
        flag = 310;

      case 310
        if (NF > NFTEST)
          NF = NF - 1;
          if (IPRINT > 0)
            warning('return from NEWUOA because CALFUN has been' + ...
            ' called MAXFUN times.');
          end
          flag = 530;
        else
          F = 0;
          if (NF >= NPT + 1)
            TEMP_ = Func(X);
            [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
            F = F + G2;
            % fid = fopen('history_F.csv', 'at');
            % G = 0.0e0;
            % fprintf(fid, '%.16f\n', F - (G2 - G));
            % fclose(fid);

            if (NF >= NPT + 2)
              % PRINT 226, (XBEG(I), I = 1, N)
              % 226 FORMAT (/ 4X, 'XBEG-in-test =', 1PD23.15)

              for I = 1:N
                XBEG(I) = 1.0e0;
              end
              TEMP_ = Func(XBEG);
              [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 2, FLAG, EPS, EPS1, EPS2, RAND);
              [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
              FBEG = FBEG + (G2 - G1);
              for I = 1:N
                XTEMP(I) = XOPT(I) + XBASE(I);
              end
              TEMP_ = Func(XTEMP);
              [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 2, FLAG, EPS, EPS1, EPS2, RAND);
              [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
              FOPT = (G2 - G1) + FOPT;

              for J = 1:NPT
                for I = 1:N
                  XTEMP(I) = XPT(J, I) + XBASE(I);
                end

                % PRINT 276, (XTEMP(I), I = 1, N)
                % 276 FORMAT (/ 4X, 'XTEMP-test =', 1PD23.15)
                TEMP_ = Func(XTEMP);
                [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 2, FLAG, EPS, EPS1, EPS2, RAND);
                [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
                FVAL(J) = (G2 - G1) + FVAL(J);
              end
            end
          else
            TEMP_ = Func(X);
            [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NPT, FLAG, EPS, EPS1, EPS2, RAND);
            F = F + G1;

            % fid = fopen('history_F.csv', 'at');
            % G = 0.0e0;
            % fprintf(fid, '%.16f\n', F - (G1 - G));
            % fclose(fid);

          end
          % disp(['FVAL-test =', num2str(FVAL)])
          % disp(['FBEG-test =', num2str(FBEG)])
          % disp(['F-test =', num2str(F)])
          % disp(['FOPT-test =', num2str(FOPT)])
          % if (IPRINT == 3)
          %   disp(['Function number ', num2str(NF), '; F =', num2str(F)])
          %   disp('The corresponding X array is:')
          %   disp(X)
          % end
          if (NF <= NPT)
            flag = 70;
          elseif (KNEW == -1)
            flag = 530;
          else
            %     Use the quadratic model to predict the change
            %     in F due to the step D,
            %     and set DIFF to the error of this prediction.
            VQUAD = ZERO;
            IH = 0;
            for J = 1:N
              VQUAD = VQUAD + D(J) * GQ(J);
              for I = 1:J
                IH = IH + 1;
                TEMP = D(I) * XNEW(J) + D(J) * XOPT(I);
                if (I == J)
                  TEMP = HALF * TEMP;
                end
                VQUAD = VQUAD + TEMP * HQ(IH);
              end
            end
            for K = 1:NPT
              VQUAD = VQUAD + PQ(K) * W(K);
            end
            % disp('XOPT-test =')
            % disp(XOPT)
            % disp(['D-test =', num2str(D)])
            % disp(['GQ-test =', num2str(GQ)])
            % disp(['PQ-test =', num2str(PQ)])
            % disp(['HQ-test =', num2str(HQ)])
            for I = 1:N
              XTEMP(I) = XOPT(I) + XBASE(I) + D(I);
            end

            if (NF >= NPT + 1)
              TEMP_ = Func(XTEMP);
              [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF, FLAG, EPS, EPS1, EPS2, RAND);
              [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
              DIFF = F - FOPT + (G2 - G1) - VQUAD;
            end
            % disp(['DIFF-0405 =', num2str(DIFF)])
            % disp(['F-0405 =', num2str(F)])
            % disp(['FOPT-0405 =', num2str(FOPT)])
            % disp(['VQUAD-0405 =', num2str(VQUAD)])
            % disp('XOPT0405 =')
            % disp(XOPT)
            DIFFC = DIFFB;
            DIFFB = DIFFA;
            DIFFA = abs((-G2 + G1) + DIFF);
            % disp(['DIFFA-0405 =', num2str(DIFFA)])
            % disp(['DIFFB-0405 =', num2str(DIFFB)])
            % disp(['DIFFC-0405 =', num2str(DIFFC)])
            if (DNORM > RHO)
              NFSAV = NF;
            end
            % disp('XOPT++++-test=')
            % disp(num2str(XOPT))
            % disp(['VQUAD-test=', num2str(VQUAD)])
            % disp('++++++++++++++++++++++END++++++++++++++++++++\n')
            % disp('+++++++++++++++++++++BEGIN+++++++++++++++++++\n')
            %     Update FOPT and XOPT if the new F is the least value of the objective
            %     function so far. The branch when KNEW is positive occurs if D is not
            %     a trust region step.
            FSAVE = FOPT;
            if (F < FOPT)
              FOPT = F;
              XOPTSQ = ZERO;
              for I = 1:N
                XOPT(I) = XNEW(I);
                XOPTSQ = XOPTSQ + XOPT(I) ^ 2;
              end
            end
            KSAVE = KNEW;
            if (KNEW > 0)
              flag = 410;
              %     Pick the next value of DELTA after a trust region step.
            elseif (VQUAD >= ZERO)
              if (IPRINT > 0)
                warning('return from NEWUOA because a trust' + ...
                ' region step has failed to reduce Q.')
              end
              flag = 530;
            else
              RATIO = (F - FSAVE) / VQUAD;
              % disp(['F-0405peng=', num2str(F), ', FSAVE-0405peng=', num2str(FSAVE)])
              % disp(['VQUAD-0403peng=', num2str(VQUAD)])
              % disp(['RATIO-0403peng=', num2str(RATIO)])
              if (RATIO <= TENTH)
                DELTA = HALF * DNORM;
              elseif (RATIO <= 0.7e0)
                DELTA = max(HALF * DELTA, DNORM);
              else
                DELTA = max(HALF * DELTA, DNORM + DNORM);
              end
              if (DELTA <= 1.5e0 * RHO)
                DELTA = RHO;
              end
              %     Set KNEW to the index of the next interpolation point to be deleted.
              RHOSQ = max(TENTH * DELTA, RHO) ^ 2;
              KTEMP = 0;
              DETRAT = ZERO;
              if (F >= FSAVE)
                KTEMP = KOPT;
                DETRAT = ONE;
              end
              for K = 1:NPT
                HDIAG = ZERO;
                for J = 1:NPTM
                  TEMP = ONE;
                  if (J < IDZ)
                    TEMP = -ONE;
                  end
                  HDIAG = HDIAG + TEMP * ZMAT(K, J) ^ 2;
                end
                TEMP = abs(BETA * HDIAG + VLAG(K) ^ 2);
                DISTSQ = ZERO;
                for J = 1:N
                  DISTSQ = DISTSQ + (XPT(K, J) - XOPT(J)) ^ 2;
                end
                if (DISTSQ > RHOSQ)
                  TEMP = TEMP * (DISTSQ / RHOSQ) ^ 3;
                end
                if (TEMP > DETRAT && K ~= KTEMP)
                  DETRAT = TEMP;
                  KNEW = K;
                end
              end
              if (KNEW == 0)
                flag = 460;
              else
                flag = 410;
              end
            end
          end
        end

        %     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
        %     can be moved. Begin the updating of the quadratic model, starting
        %     with the explicit second derivative term.
      case 410
        [BMAT_, ZMAT_, IDZ_, VLAG_, W_] = ...
          update (N, NPT, BMAT, ZMAT, IDZ, VLAG, BETA, KNEW, W);
        W = W_;
        VLAG = VLAG_;
        BMAT = BMAT_;
        ZMAT = ZMAT_;
        IDZ = IDZ_;
        FVAL(KNEW) = F;
        IH = 0;
        for I = 1:N
          TEMP = PQ(KNEW) * XPT(KNEW, I);
          for J = 1:I
            IH = IH + 1;
            HQ(IH) = HQ(IH) + TEMP * XPT(KNEW, J);
          end
        end
        PQ(KNEW) = ZERO;
        %     Update the other second derivative parameters, and  the gradient
        %     vector of the model. Also include the new interpolation point.
        for J = 1:NPTM
          for JJ = 1:NPT
            if (JJ == KNEW)
              TEMP1(JJ) = DIFF * ZMAT(JJ, J);
            else
              for K = 1:N
                XTEMP(K) = XBASE(K) + XPT(JJ, K);
              end

              if (NF == NPT + 1)
                TEMP_ = Func(XTEMP);
                [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF, FLAG, EPS, EPS1, EPS2, RAND);
                [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
                TEMP1(JJ) = (G2 - G1) * ZMAT(JJ, J);
              else
                TEMP1(JJ) = 0;
              end
            end

            if (J < IDZ)
              TEMP1(JJ) = -TEMP1(J);
            end
            for K = 1:NPT
              PQ(K) = PQ(K) + TEMP1(JJ) * ZMAT(K, J);
            end
          end
        end
        GQSQ = ZERO;
        for I = 1:N
          XPT(KNEW, I) = XNEW(I);
          for J = 1:NPT
            if (J == KNEW)
              GQ(I) = GQ(I) + DIFF * BMAT(J, I);
              % DO K = 1, N
              % XTEMP(K) = XBASE(K) + XPT(J, K)
              % END DO
              % Gunc(XTEMP, NF - 1, G2)
              % Gunc(XTEMP, NF - 2, G1)
            else
              for K = 1:N
                XTEMP(K) = XBASE(K) + XPT(J, K);
              end

              if (NF >= NPT + 1)
                TEMP_ = Func(XTEMP);
                [G2, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF, FLAG, EPS, EPS1, EPS2, RAND);
                [G1, FLAG, EPS, EPS1, EPS2, RAND] = Gunc(TEMP_, NF - 1, FLAG, EPS, EPS1, EPS2, RAND);
                GQ(I) = GQ(I) + (G2 - G1) * BMAT(J, I);
              end

            end
          end
          GQSQ = GQSQ + GQ(I) ^ 2;
          XPT(KNEW, I) = XNEW(I);
        end
        %     if a trust region step makes a small change to the objective function,
        %      calculate the gradient of the least Frobenius norm interpolant at
        %     XBASE, and store it in W, using VLAG for a vector of right hand sides.
        if (KSAVE == 0 && DELTA == RHO)
          if (abs(RATIO) > 1.0e-2)
            ITEST = 0;
          else
            for K = 1:NPT
              VLAG(K) = FVAL(K) - FVAL(KOPT);
            end
            GISQ = ZERO;
            for I = 1:N
              SUM = ZERO;
              for K = 1:NPT
                SUM = SUM + BMAT(K, I) * VLAG(K);
              end
              GISQ = GISQ + SUM * SUM;
              W(I) = SUM;
            end

            %     Test whether to replace the new quadratic model by the least Frobenius
            %     norm interpolant, making the replacement if the test is satisfied.

            ITEST = ITEST + 1;
            if (GQSQ < 1.0e2 * GISQ)
              ITEST = 0;
            end
            if (ITEST >= 3)
              for I = 1:N
                GQ(I) = W(I);
              end
              for IH = 1:NH
                HQ(IH) = ZERO;
              end
              for J = 1:NPTM
                W(J) = ZERO;
                for K = 1:NPT
                  W(J) = W(J) + VLAG(K) * ZMAT(K, J);
                end
                if (J < IDZ)
                  W(J) = -W(J);
                end
              end
              for K = 1:NPT
                PQ(K) = ZERO;
                for J = 1:NPTM
                  PQ(K) = PQ(K) + ZMAT(K, J) * W(J);
                end
              end
              ITEST = 0;
            end
          end
        end
        if (F < FSAVE)
          KOPT = KNEW;
        end
        %     if a trust region step has provided a sufficient decrease in F,
        %     branch for another trust region calculation. The case KSAVE>0 occurs
        %     when the new function value was calculated by a model step.
        if (F <= FSAVE + TENTH * VQUAD || KSAVE > 0)
          flag = 100;
          %     Alternatively, find out if the interpolation points are close enough
          %     to the best point so far.
        else
          KNEW = 0;
          flag = 460;
        end

      case 460
        DISTSQ = 4.0e0 * DELTA * DELTA;
        for K = 1:NPT
          SUM = ZERO;
          for J = 1:N
            SUM = SUM + (XPT(K, J) - XOPT(J)) ^ 2;
          end
          if (SUM > DISTSQ)
            KNEW = K;
            DISTSQ = SUM;
          end
        end
        %     if KNEW is positive,  set DSTEP, and branch back for the next
        %     iteration, which will generate a "model step".
        if (KNEW > 0)
          DSTEP = max(min(TENTH * sqrt(DISTSQ), HALF * DELTA), RHO);
          DSQ = DSTEP * DSTEP;
          flag = 120;
        elseif (RATIO > ZERO || max(DELTA, DNORM) > RHO)
          flag = 100;
        else
          flag = 490;
        end
        %     The calculations with the current value of RHO are complete. Pick the
        %     next values of RHO and DELTA.

      case 490
        if (RHO > RHOEND)
          % disp(['kankan RHO:', num2str(RHO), ', kankan Delta:', num2str(DELTA)])
          DELTA = HALF * RHO;
          RATIO = RHO / RHOEND;
          if (RATIO <= 16.0e0)
            RHO = RHOEND;
          elseif (RATIO <= 250.0e0)
            RHO = sqrt(RATIO) * RHOEND;
          else
            RHO = TENTH * RHO;
          end
          DELTA = max(DELTA, RHO);
          % if (IPRINT >= 2)
          %   if (IPRINT >= 3)
          %     disp("(5X)")
          %   end
          %   disp(['New RHO =', num2str(RHO), ', Current number of function evaluations =', num2str(NF)])
          %   disp(['Least value of F = ', num2str(FOPT), ', The corresponding X array is:']); disp(XBASE + XOPT)
          % end
          flag = 90;
        else
          %     return from the calculation, after another Newton-Raphson step, if
          %     it is too short to have been tried before.
          % disp(['NF = ', num2str(NF)])
          if (KNEW == -1)
            flag = 290;
          else
            flag = 530;
          end
        end

      otherwise

        if (FOPT <= F)
          for I = 1:N
            X(I) = XBASE(I) + XOPT(I);
          end
          F = FOPT;
        end
        % if (IPRINT >= 1)
        % disp('At the return from NEWUOA ')
        % disp(['Total times of function evaluations =', num2str(NF)])
        % disp(['Least value of F =', num2str(F)])
        % disp('The corresponding X array is:'); disp(X)
        % end
        flag = 1111;
    end
  end

end
