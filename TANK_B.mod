@#define use_dll = 0
options_.use_dll = 0;

warning off;

// ======================
// VARIABLES ENDÓGENAS
// ======================
var 
    // Hogares separados
    C_G C_B            // consumos de hogares G y B
    L_B L_G            // horas de hogares G y B
    M_G M_B            // SDF (GHH) de hogares G y B

    // Mercado laboral por sector (usadas por firmas)
    L_b L_g              // horas demandadas por firmas brown/green
    wb wg              // salarios por sector

    // Política monetaria y precios agregados
    R pi P

    // Precios de inversión y retornos del capital por sector
    Qb Qg Rkb Rkg

    // Producción, precios sectoriales y agregados
    Yb Yg Y pb pg

    // Capitales e inversión
    Kb Kg Ib Ig I

    // Medio ambiente
    emis X dX A

    // Métricas financieras/auxiliares
    spreadb spreadg sg
;

// ======================
// SHOCKS EXÓGENOS
// ======================
varexo eA eR;

// ======================
// PARÁMETROS
// ======================
parameters 
    betta eta ksi omega rhoL 
    alfab alfag rhoy alfayb alfayg deltab deltag phi phi_pi
    d0 d1 d2 damage_scale deltax erow 
    rhoA sigmaA Abar
    rho_i sigmaR
     phi_p eps_p tau
;

// ===== Valores de parámetros =====
betta  = 0.9975;        // (si quieres algo más estándar, usa 0.9975)
eta    = 2;
ksi    = 1;
omega  = 1.01;
rhoL   = 1;

tau = 0.05;

alfab  = 0.35;
alfag  = 0.33;

rhoy   = 2;
alfayg = 0.668;
alfayb = 1 - alfayg;

deltab = 0.025;
deltag = 0.025;

phi    = 10;

phi_pi = 2;        // respuesta a inflación (Taylor >= 1)
rho_i  = 0.75;       // suavizamiento
sigmaR = 0.002;

damage_scale = 0.9860;
d0     = -0.026;
d1     = 3.61e-5;
d2     = 1.44e-8;
// reescalado por damage_scale
d2 = d2 / (damage_scale^2);
d1 = d1 / damage_scale;

deltax = 0.9965;
erow   = 3.3705;

Abar   = 1;
rhoA   = 0.95;
sigmaA = 0.007;

phi_p = 50;     // costo cuadrático de ajustar precios (↑ phi_p = precios más pegajosos)
eps_p = 6;      // elasticidad de sustitución entre variedades (markup deseado = eps_p/(eps_p-1))


model;

// === Hogares separados por sector (G y B) ===

// 1) Enlaces horas (alias para que firmas usen Lb, Lg)
L_b = L_B;
L_g = L_G;

// 2) Factores de descuento (GHH) por hogar
M_G = betta * ((C_G - omega*(L_G^(1+ksi))/(1+ksi))^-eta)
             / ((C_G(-1) - omega*(L_G(-1)^(1+ksi))/(1+ksi))^-eta);

M_B = betta * ((C_B - omega*(L_B^(1+ksi))/(1+ksi))^-eta)
             / ((C_B(-1) - omega*(L_B(-1)^(1+ksi))/(1+ksi))^-eta);

// Restricciones presupuestarias
C_B + Ib + tau * Yb = wb * L_B + Rkb * Qb(-1) * Kb(-1);
C_G + Ig       = wg * L_G + Rkg * Qg(-1) * Kg(-1) + tau * Yb;

// 3) Índice de precio CES (dual) y 4) Inflación
P  = ( alfayb*pb^(1 - rhoy) + alfayg*pg^(1 - rhoy) )^( 1/(1 - rhoy) );
pi = P / P(-1);

// 5) Regla monetaria (π-only, smoothing, ancla 1/β)
log(R) = (1 - rho_i)*log(1/betta) + rho_i*log(R(-1)) + phi_pi*log(pi) + sigmaR*eR;

// --- (NKPC Rotemberg agregada) ---
// π_t = β E_t[π_{t+1}] + ((ε-1)/φ_p) * (mc_t - 1),
// con mc_t aproximado por el coste laboral unitario real (ULC) sobre PY.
// NKPC Rotemberg agregada (referencia correcta al mc_ss = (eps_p-1)/eps_p )
// Reemplazo sugerido (misma κ):
log(pi) = betta*log(pi(+1)) + ((eps_p - 1)/phi_p) * ( (wb*L_b + wg*L_g)/(P*Y) - (eps_p-1)/eps_p );


// 6–7) FOCs laborales separadas por hogar/sector
// Nota: omega*L_i^(ksi-rhoL)*L_i^rhoL = omega*L_i^ksi  (puedes dejar la forma expandida si prefieres)
omega * L_B^(ksi-rhoL) * L_B^rhoL = wb;
omega * L_G^(ksi-rhoL) * L_G^rhoL = wg;


// 10–11) Contaminación y daño
X  = deltax*X(-1) + emis + erow;
dX = d0 + d1*X + d2*X^2;

// 12) Emisiones (sin abatimiento): proporcionales a Yb
emis = tau * Yb;

// 13–15) Sector brown (firmas)
Yb  = (1-tau)*(1 - dX)*A*(Kb(-1)^alfab)*(L_b^(1 - alfab));
wb  = (1 - alfab)*pb*(Yb / L_b);
Rkb = (alfab*pb*(Yb/Kb(-1)) + (1 - deltab)*Qb) / Qb(-1);

// 16–18) Sector green (firmas)
Yg  = (1 - dX)*A*(Kg(-1)^alfag)*(L_g^(1 - alfag));
wg  = (1 - alfag)*pg*Yg / L_g;
Rkg = (pg*alfag*Yg/Kg(-1) + (1 - deltag)*Qg) / Qg(-1);

// 19–21) Agregación y precios (CES)
Y  = ((alfayb^(1/rhoy))*Yb^((rhoy-1)/rhoy) + (alfayg^(1/rhoy))*Yg^((rhoy-1)/rhoy))^(rhoy/(rhoy-1));
pb = (alfayb * Y / Yb)^(1 / rhoy);
pg = (alfayg * Y / Yg)^(1 / rhoy);

// 22–24) Acumulación de capital e inversión
Kb = (1 - deltab) * Kb(-1) + Ib;
Kg = (1 - deltag) * Kg(-1) + Ig;
I  = Ib + Ig;

// 25–26) Precios de inversión con coste de ajuste (descuento por hogar relevante)
Qb = 1 + (phi/2)*(Ib/Ib(-1) - 1)^2
       + phi*(Ib/Ib(-1) - 1)*(Ib/Ib(-1))
       - M_B(+1)*phi*(Ib(+1)/Ib - 1)*(Ib(+1)/Ib)^2;

Qg = 1 + (phi/2)*(Ig/Ig(-1) - 1)^2
       + phi*(Ig/Ig(-1) - 1)*(Ig/Ig(-1))
       - M_G(+1)*phi*(Ig(+1)/Ig - 1)*(Ig(+1)/Ig)^2;

// 27) Recurso agregado
// 27) Recurso agregado (añadimos coste Rotemberg (phi_p/2)*(pi-1)^2 * Y)
Y = C_G + C_B + I
  + (phi/2)*((Ib/Ib(-1) - 1)^2)*Ib
  + (phi/2)*((Ig/Ig(-1) - 1)^2)*Ig
  + (phi_p/2)*(pi - 1)^2 * Y;


// 28) TFP
log(A) = (1 - rhoA)*log(Abar) + rhoA*log(A(-1)) + sigmaA*eA;

// 29–31) Spreads y participación green (definiciones)
spreadb = Rkb(+1) - R;
spreadg = Rkg(+1) - R;
sg = Qg*Kg / (Qb*Kb + Qg*Kg);


end;

// ======================
// VALORES INICIALES
// ======================
initval;
  A  = Abar;  P = 1;
  // No fijes pi=1 si mantienes NKPC en niveles; si la pasas a logs, puedes dejar pi=1
  pi = 1;
  R  = 1/betta;

  // Horas (oferta y alias)
  L_B = 0.5;   L_G = 0.5;
  L_b = L_B;   L_g = L_G;

  // Sectores (arranque)
  Yb = 1; Yg = 1;  pb = 1; pg = 1;

  // Capitales de Euler
  Kb = 12.73;
  Kg = 12.00;
  Ib = deltab*Kb;   Ig = deltag*Kg;   I = Ib + Ig;

  // Salarios por FOC de firma
  wb = (1-alfab)*pb*(Yb/L_b);
  wg = (1-alfag)*pg*(Yg/L_g);

  // Agregado (puedes omitir esta línea y que lo cierre el CES)
  // Y = ((alfayb^(1/rhoy))*Yb^((rhoy-1)/rhoy) + (alfayg^(1/rhoy))*Yg^((rhoy-1)/rhoy))^(rhoy/(rhoy-1));

  // Clima
  emis = Yb;
  X    = (emis + erow) / (1 - deltax);
  dX   = d0 + d1*X + d2*X^2;

  // Consumos coherentes con recursos
  // (ajusta para que C_G + C_B = Y - I exactamente)
  C_G = 0.300;
  C_B = 0.300;

  // Qs
  Qb = 1; Qg = 1;

  spreadb = 0; spreadg = 0;
  sg = Qg*Kg / (Qb*Kb + Qg*Kg);
    

   M_G = betta; M_B = betta;

end;



shocks;
  var eA; stderr +1;
  var eR; stderr +1;
end;


% Opciones del steady (si no usas SS externo, igual ayudan)
options_.solve_algo   = 4;
options_.steady.maxit = 500;
options_.steady.tol   = 1e-12;

resid;
steady;     % usa TANK_A_steadystate.m automáticamente
check;

graph_format (pdf, png);  % guarda ambos formatos
nodisplay;                % no abre ventanas, solo guarda archivos

% IRFs (ajusta horizons a gusto)
stoch_simul(order=1, irf=40);
% o si quieres además series simuladas:
% stoch_simul(order=1, periods=300, irf=40);