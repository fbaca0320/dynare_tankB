function [ys, params, check] = TANK_A_steadystate(ys, exo, M_, options_)
% Steady state para tu TANK_A.mod (dos hogares/sectores, Rotemberg en log)
% - No cambia NINGUNA ecuación del .mod.
% - Resuelve Kb,Kg,L_B,L_G con Eulers + FOCs laborales + daño climático.
% - Fija pi=1, Qb=Qg=1 en SS.
% - Ajusta eps_p en el SS para que la NKPC en log quede EXACTA (mc_ss=(eps_p-1)/eps_p).
% - Construye C_G, C_B cumpliendo recursos (split por cuota salarial).
%
% Requisitos: Dynare 6.x, MATLAB Optim Toolbox (fsolve).

  if nargin < 4, options_ = struct(); end %#ok<NASGU>
  check  = 0;
  ys     = zeros(M_.endo_nbr,1);
  params = M_.params;

  % ===== Helpers =====
  EN   = cellstr(M_.endo_names);
  PN   = cellstr(M_.param_names);
  geti = @(nm) find(strcmp(nm, EN), 1);
  getp = @(nm) get_param_by_name(nm);
  setp = @(nm,val) assign_param(nm,val);

  % ===== Lee parámetros relevantes =====
  betta  = getp('betta');   ksi    = getp('ksi');    omega  = getp('omega'); %#ok<NASGU>
  alfab  = getp('alfab');   alfag  = getp('alfag');
  rhoy   = getp('rhoy');    alfayb = getp('alfayb'); alfayg = getp('alfayg');
  deltab = getp('deltab');  deltag = getp('deltag');
  d0     = getp('d0');      d1     = getp('d1');     d2     = getp('d2');
  deltax = getp('deltax');  erow   = getp('erow');
  Abar   = getp('Abar');

  % (phi_p, eps_p se manejarán al final; phi de inversión no afecta el SS)
  R_ss   = 1/betta;
  pi_ss  = 1;   P_ss = 1;   Qb_ss = 1;  Qg_ss = 1;

  % ===== Multistart para el sistema 4x4 (Kb, Kg, L_B, L_G) =====
  base_guess = [12.7; 12.0; 0.60; 0.40];
  GUESS_SET = [ ...
     1.00  1.00  1.00  1.00 ;
     1.10  0.90  1.05  0.95 ;
     0.90  1.10  0.95  1.05 ;
     1.20  1.20  0.90  0.90 ;
     0.80  0.80  1.10  1.10 ];
  opts = optimoptions('fsolve','Display','off','TolFun',1e-14,'TolX',1e-14,...
                      'MaxIter',2000,'MaxFunEvals',8000);

  solved = false; z = base_guess;
  for k=1:size(GUESS_SET,1)
    z0 = base_guess .* GUESS_SET(k,:).';
    [z_try, ~, exitflag] = fsolve(@(z) sys4(z, Abar, alfab, alfag, ...
        rhoy, alfayb, alfayg, deltab, deltag, d0, d1, d2, deltax, erow, ...
        betta), z0, opts);
    if exitflag > 0 && all(isfinite(z_try))
      z = z_try; solved = true; break;
    end
  end
  if ~solved
    check = 1; return;  % deja que Dynare reporte fallo si no converge
  end

  % ===== Variables clave resueltas =====
  Kb = max(1e-8, z(1));
  Kg = max(1e-8, z(2));
  L_B= min(0.9999, max(1e-8, z(3)));
  L_G= min(0.9999, max(1e-8, z(4)));

  % Construye todo el nivel (Yb,Yg,Y, precios, salarios, clima)
  [Yb,Yg,Y,pb,pg,P,wb,wg,X,dX] = build_levels(Kb,Kg,L_B,L_G, ...
                                  Abar,alfab,alfag,rhoy,alfayb,alfayg, ...
                                  d0,d1,d2,deltax,erow);

  % Inversión y retornos
  Ib  = deltab*Kb;  Ig  = deltag*Kg;  I = Ib+Ig;
  Rkb = R_ss;       Rkg = R_ss;       emis = Yb;

  % ===== Ajuste de eps_p para que la NKPC en log sea EXACTA en SS =====
  % log(pi)=beta*log(pi(+1)) + kappa*(mc - (eps-1)/eps)  con pi=1 => 0 = kappa*(mc - (eps-1)/eps)
  % => (eps-1)/eps = mc  =>  eps = 1/(1-mc)
  mc = (wb*L_B + wg*L_G) / max(1e-12, (P*Y));   % ULC / PY
  mc = min(0.98, max(0.02, mc));                % guarda márgenes
  eps_new = 1 / (1 - mc);                       % >1 asegurado con el clamp
  % Escribe el parámetro eps_p en 'params' de Dynare
  idx_eps = find(strcmp('eps_p', PN), 1);
  if ~isempty(idx_eps)
    params(idx_eps) = eps_new;
  end

  % ===== Consumo total y reparto (recursos; Rotemberg=0 con pi=1) =====
  Ctot = Y - I;
  wagebill = wb*L_B + wg*L_G;
  thetaC   = min(0.95, max(0.05, wagebill/max(1e-12, (P*Y))));
  C_G = thetaC * Ctot;
  C_B = (1 - thetaC) * Ctot;

  % SDF en SS (por tus ecuaciones de M_i)
  M_G = betta;  M_B = betta;

  % ===== Asignación a ys (TODAS las endógenas) =====
  ys(geti('C_G')) = C_G;    ys(geti('C_B')) = C_B;
  ys(geti('L_B')) = L_B;    ys(geti('L_G')) = L_G;
  ys(geti('M_G')) = M_G;    ys(geti('M_B')) = M_B;

  ys(geti('L_b')) = L_B;    ys(geti('L_g')) = L_G;
  ys(geti('wb'))  = wb;     ys(geti('wg'))  = wg;

  ys(geti('R')) = R_ss;     ys(geti('pi')) = pi_ss;   ys(geti('P')) = P;

  ys(geti('Qb')) = Qb_ss;   ys(geti('Qg')) = Qg_ss;
  ys(geti('Rkb'))= Rkb;     ys(geti('Rkg'))= Rkg;

  ys(geti('Kb')) = Kb;      ys(geti('Kg')) = Kg;
  ys(geti('Ib')) = Ib;      ys(geti('Ig')) = Ig;      ys(geti('I')) = I;

  ys(geti('Yb')) = Yb;      ys(geti('Yg')) = Yg;      ys(geti('Y')) = Y;
  ys(geti('pb')) = pb;      ys(geti('pg')) = pg;

  ys(geti('emis')) = emis;  ys(geti('X'))  = X;       ys(geti('dX')) = dX;
  ys(geti('A'))    = Abar;

  sg = (Qg_ss*Kg)/(Qb_ss*Kb + Qg_ss*Kg);
  ys(geti('spreadb')) = 0;
  ys(geti('spreadg')) = 0;
  ys(geti('sg'))      = sg;

  % Señal de fallo si aparece algo no finito
  if any(~isfinite(ys)) || any(~isfinite(params))
    check = 1;
  end

  % ===== Nested: sistema 4x4 =====
  function F = sys4(z, Abar, alfab, alfag, rhoy, alfayb, alfayg, ...
                    deltab, deltag, d0, d1, d2, deltax, erow, betta)
    Kb_ = max(1e-10, z(1));  Kg_ = max(1e-10, z(2));
    Lb_ = min(0.9999, max(1e-10, z(3)));
    Lg_ = min(0.9999, max(1e-10, z(4)));

    [Yb_,Yg_,~,pb_,pg_,~,wb_,wg_,~,~] = build_levels(Kb_,Kg_,Lb_,Lg_, ...
        Abar,alfab,alfag,rhoy,alfayb,alfayg, d0,d1,d2,deltax,erow);

    Rss = 1/betta;
    % Euler por sector en SS: α_i p_i (Y_i/K_i) + (1-δ_i) = R_ss
    e1 = alfab*pb_*(Yb_/Kb_) + (1 - deltab) - Rss;
    e2 = alfag*pg_*(Yg_/Kg_) + (1 - deltag) - Rss;
    % FOC laboral hogar=firma: w = ω L^ξ  vs. w = (1-α)p Y/L
    e3 = wb_ - (getp('omega')*(Lb_^getp('ksi')));  % w_b - ω L_B^ξ
    e4 = wg_ - (getp('omega')*(Lg_^getp('ksi')));  % w_g - ω L_G^ξ
    F  = [e1; e2; e3; e4];
  end

  % ===== Construye niveles (K,L)->(Y,p,w,P,X) =====
  function [Yb,Yg,Y,pb,pg,P,wb,wg,X,dX] = build_levels(Kb,Kg,Lb,Lg, ...
      Abar,alfab,alfag,rhoy,alfayb,alfayg,d0,d1,d2,deltax,erow)

    A = Abar;

    % Daño endógeno con X = (Yb + erow)/(1 - deltax); dX = d0 + d1 X + d2 X^2
    den = (1 - deltax);
    a2 = d2/(den^2);
    a1 = d1/den + (2*d2*erow)/(den^2);
    a0 = d0 + (d1*erow)/den + (d2*(erow^2))/(den^2);

    T_b = A * (Kb^alfab) * (Lb^(1 - alfab));
    qa = a2*T_b; qb = (1 + a1*T_b); qc = -T_b*(1 - a0);
    disc = qb^2 - 4*qa*qc; if disc < 0, disc = 0; end
    if abs(qa) < 1e-14
      Yb = (T_b*(1 - a0)) / max(1e-12,(1 + a1*T_b));
    else
      Yb1 = (-qb + sqrt(disc))/(2*qa);
      Yb2 = (-qb - sqrt(disc))/(2*qa);
      Yb  = max(Yb1, Yb2);
    end
    Yb = max(1e-12, Yb);

    X  = (Yb + erow) / (1 - deltax);
    dX = d0 + d1*X + d2*(X^2);

    T_g = A * (Kg^alfag) * (Lg^(1 - alfag));
    Yg  = max(1e-12, (1 - dX) * T_g);

    term = ( (alfayb^(1/rhoy))*Yb^((rhoy-1)/rhoy) + ...
             (alfayg^(1/rhoy))*Yg^((rhoy-1)/rhoy) );
    Y    = max(1e-12, term^(rhoy/(rhoy-1)));

    pb = (alfayb * Y / Yb)^(1/rhoy);
    pg = (alfayg * Y / Yg)^(1/rhoy);
    P  = ( alfayb*pb^(1 - rhoy) + alfayg*pg^(1 - rhoy) )^( 1/(1 - rhoy) );

    wb = (1 - alfab) * pb * (Yb / Lb);
    wg = (1 - alfag) * pg * (Yg / Lg);
  end

  % ===== set param helper (escribe en 'params' por nombre) =====
  function assign_param(nm, val)
    j = find(strcmp(nm, PN), 1);
    if ~isempty(j), params(j) = val; end
  end
end
