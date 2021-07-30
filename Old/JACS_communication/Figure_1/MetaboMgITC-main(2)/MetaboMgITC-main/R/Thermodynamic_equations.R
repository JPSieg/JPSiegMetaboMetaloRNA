#'Function that describes heat produced by each injection in an ITC experiment
#'
#'Function that describes heat produced by each injection in an ITC experiment when there is 1 to 1 binding of ligand to a macromolecule and a competator.
#'
#'@param V Cell volume in Liters
#'@param X Ligand concentration in a cell
#'@param M Macromolecule in a cell
#'@param H Enthalpy of a ligand binding to a macromolecule
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@return Heat produced by an injection
#' @export
Wiseman.isotherm.competative = function(V, X, M, C, H = 2534, K = 2880, Hc = 2534, Kc = 2880, c = 0, h = 0){
  x = X/M
  Kapp = K*(1 + c*Kc*C)/(1 + Kc*C)
  Happ = H - ((Hc*Kc*C)/(1+Kc*C)) + ((Hc + h)*(x*Kc*C)/(1 + c*Kc*C))
  r = 1/(Kapp*M)
  y = 0.5 + ((1 - x - r)/(2*sqrt(((1 + x + r)^2)-(4*x))))
  dQ.dX = Happ*V*y
}

#'Function that describes heat produced by each injection in an ITC experiment
#'
#'Function that describes heat produced by each injection in an ITC experiment when there is 1 to 1 binding of ligand to a macromolecule.
#'
#'@param V Cell volume in Liters
#'@param X Ligand concentration in a cell
#'@param M Macromolecule in a cell
#'@param H Enthalpy of a ligand binding to a macromolecule
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@return Heat produced by an injection
#' @export
Wiseman.isotherm = function(V, X, M, H = 2534, K = 2880){
  x = X/M
  r = 1/(K*M)
  y = 0.5 + ((1 - x - r)/(2*sqrt(((1 + x + r)^2)-(4*x))))
  dQ.dX = H*V*y
}

#'Function that calculates the Gibbs free energy of a chemical reaction
#'
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@param Temperature Temperature in celcius
#'@param R Gas constant. Default = 0.0019872 kcal/mol/K
#'@return Gibbs free energy of a chemical reaction
#' @export
Gibbs.free.energy = function(K, Temperature, R = 0.0019872){
  output = -R*(Temperature + 273.15)*log(K)
}

#'Function that calculates the Entropy of a chemical reaction
#'
#'@param dG Gibbs free energy of a chemical reaction
#'@param dH Enthalpy of a chemical reaction
#'@param Temperature Temperature in celcius
#'@return Entropy of a chemical reaction
#' @export
Entropy = function(dG, dH, Temperature){
  output = (dH - dG)/Temperature
}

#'Function that calculates the saturation of a macromolecule in an ITC experiment
#'
#'Function that describes heat produced by each injection in an ITC experiment when there is 1 to 1 binding of ligand to a macromolecule.
#'
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@param M Macromolecule in a cell
#'@param X Ligand concentration in a cell
#'@param V Cell volume in Liters
#'@param V0 Original cell volume
#'@return Percent saturation of the macromolecule by bound ligand
#' @export
saturation = function(K, M, X, V, V0){
  M = (V*M)/(1.4 + 0.282)
  X = (0.282*X)/(1.4 + 0.282)
  a = K
  b = K*X - K*M + 1
  c = -M
  M.free = (-b + sqrt(((b^2)-(4*a*c))))/(2*a)
  output = 100 - (100*M.free/M)
  print(output)
}

#'Function that calculates the log10 of a metal ion affinity constant for Ionic strength = 0 using the Debye-Huckel equation
#'
#'@param log10K Molar affinity constant of a metabolite binding to a metal ion
#'@param I ionic strength
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param MX.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.0.calculator = function(log10K, I, M.charge, X.charge, XM.charge, A = 0.524){
  y.XM = 10^(0.1*(XM.charge^2)*I - (A*(XM.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.M = 10^(0.1*(M.charge^2)*I - (A*(M.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  log10K.0 = log10K - log10(y.XM/(y.M*y.X))
}

#'Function that calculates the log10 of a metal ion affinity constant for any Ionic strength  using the Debye-Huckel equation
#'
#'@param log10K.ref Molar affinity constant of a metabolite binding to a metal ion that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param XM.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.IS.calculator = function(log10K.ref, I.ref, I, M.charge, X.charge, XM.charge, A = 0.524){
  y.XM = 10^(0.1*(XM.charge^2)*I.ref - (A*(XM.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.M = 10^(0.1*(M.charge^2)*I.ref - (A*(M.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.X = 10^(0.1*(X.charge^2)*I.ref - (A*(X.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  log10K.0 = log10K.ref - log10(y.XM/(y.M*y.X))
  y.XM = 10^(0.1*(XM.charge^2)*I - (A*(XM.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.M = 10^(0.1*(M.charge^2)*I - (A*(M.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  log10K = log10K.0 + log10(y.XM/(y.M*y.X))
}

#'Function that calculates the pKa of any Ionic strength using the Debye-Huckel equation
#'
#'@param pKa.ref pKa for the reference that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param X.charge charge of the base
#'@param HX.charge Charge of acid
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
pKa.IS.calculator = function(pKa.ref, I.ref, I, X.charge, HX.charge, A = 0.524){
  y.HX = 10^(0.1*(HX.charge^2)*I.ref - (A*(HX.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.H = 10^(0.1*I.ref - (A*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.X = 10^(0.1*(X.charge^2)*I.ref - (A*(X.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  pKa.0 = pKa.ref - log10(y.HX/(y.H*y.X))
  print("pKa.0")
  print(pKa.0)
  y.HX = 10^(0.1*(HX.charge^2)*I - (A*(HX.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.H = 10^(0.1*I - (A*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  pKa = pKa.0 + log10(y.HX/(y.H*y.X))
  print("pKa")
  print(pKa)
}

#'Function that calculates the log10 of a metal ion affinity constant for any Ionic strength  using the Debye-Huckel equation
#'
#'@param log10K.ref Molar affinity constant of a metabolite binding to a metal ion that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param XM.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.IS.calculator.2 = function(log10K.ref, I.ref, I, M.charge, X.charge, XM.charge, A = 0.51){
  z = ((X.charge)^2 + (M.charge)^2 - (XM.charge)^2)
  log10K.0 = log10K.ref + (z*A*sqrt(I.ref))/(1 + 1.5*sqrt(I.ref))
  log10K = log10K.0 - + (z*A*sqrt(I))/(1 + 1.5*sqrt(I))
}

#'Function that calculates the mM apparent disassociation constant (Kd)
#'
#'Corrects for ionic strength using the using Specific Interaction Theory (SIT)(1) and corrects for
#'pH using pKa's to calculate populations of metal ion binding competent protonation states(2).
#'
#'[1] 	Scatchard, George. Concentrated Solutions of Strong Electrolytes. Chem. Rev. 1936, 19 (3), 309–327. https://doi.org/10.1021/cr60064a008.
#'[2] 	Mattocks, J. A.; Tirsch, J. L.; Cotruvo, J. A. Chapter Two - Determination of Affinities of Lanthanide-Binding Proteins Using Chelator-Buffered Titrations. In Methods in Enzymology; Cotruvo, J. A., Ed.; Rare-Earth Element Biochemistry: Characterization and Applications of Lanthanide-Binding Biomolecules; Academic Press, 2021; Vol. 651, pp 23–61. https://doi.org/10.1016/bs.mie.2021.01.044.
#'
#'@param metabolite The metabolite you want to model an apparant Kd for
#'@param pH The pH you want to model an apparant Kd for
#'@param I The ionic strength you want to model a Kd for
#'@param M.charge The charge of the metal ion
#'@param constants.path Path to the file containing critical stability constants for metal ion and proton binding
#'@param A Constant for SIT at 25 degC
#'@return An apparant Kd in mM (M/1000)
#' @export
Kd.app.calc = function(metabolite = "Glutathione",
                       pH = 7.5,
                       I = 0.165,
                       M.charge = 2,
                       A = 0.524,
                       constants.path = "Binding_constant_concentration_data/210525_Metaboites_binding_Mg_thermodynamics.csv"){
  df = read.csv(constants.path)
  df = subset(df, df$Metabolite == metabolite)
  if (length(df$Metabolite) == 0){
    print(paste(metabolite, "is not availible in our database"))
  }else{
    ####Correct pKas for ionic strength####
    df.is = data.frame("IS" = c(0.05, 0.10, 0.15, 0.2, 0.5, 1.0, 2.0, 3.0),
                       "Correction" = c(0.09, 0.11, 0.12, 0.13, 0.15, 0.14, 0.11, 0.07))
    correction = df.is$Correction[findInterval(I, df.is$IS)]
    df$Log_K[-c(which(df$Equilibrium == "logK_ML/M.L"), which(df$Equilibrium == "logK_MHL/M.HL"))] = df$Log_K[-c(which(df$Equilibrium == "logK_ML/M.L"), which(df$Equilibrium == "logK_MHL/M.HL"))] + correction
    ####Determine if monoprotic####

    df.pKa = df[-c(which(df$Equilibrium == "logK_ML/M.L"), which(df$Equilibrium == "logK_MHL/M.HL")),]
    if (length(which(is.na(df.pKa$Log_K))) >= 2){
      monoprotic = TRUE
    }else{
      monoprotic = FALSE
    }

    ####Calculate a.L & a.HL####

    if (monoprotic){
      pKa = df.pKa$Log_K[which(df.pKa$Equilibrium == "logK_HL/H.L")]
      a.L = 1/(1 + (10^(pKa - pH)))
      a.HL = 1 - a.L
    }else{
      pKa.HL = df.pKa$Log_K[which(df.pKa$Equilibrium == "logK_HL/H.L")]
      pKa.H2L = df.pKa$Log_K[which(df.pKa$Equilibrium == "logK_H2L/H.HL")]
      a.L = 1/(1 + 10^(pKa.HL - pH) + 10^(pKa.H2L + pKa.HL - (2*pH)))
      a.HL = 1/(1 + 10^(pH - pKa.HL) + 10^(pKa.H2L - pH))
    }

    ####Specific Interaction Theory (SIT) for IS####

    L.charge = df$X.charge[1]
    HL.charge = L.charge + 1

    I.ref.L = as.numeric(as.character(df$Ionic.strength[which(df$Equilibrium == "logK_ML/M.L")]))
    I.ref.HL = as.numeric(as.character(df$Ionic.strength[which(df$Equilibrium == "logK_MHL/M.HL")]))

    log10K.ref.ML.over.M.L = as.numeric(as.character(df$Log_K[which(df$Equilibrium == "logK_ML/M.L")]))
    log10K.ref.MHL.over.M.HL = as.numeric(as.character(df$Log_K[which(df$Equilibrium == "logK_MHL/M.HL")]))

    z.L = ((L.charge)^2 + (M.charge)^2 - (L.charge + M.charge)^2)
    z.HL = ((HL.charge)^2 + (M.charge)^2 - (L.charge + M.charge)^2)

    log10K.ML.over.M.L.0 = log10K.ref.ML.over.M.L + (z.L*A*sqrt(I.ref.L))/(1 + 1.5*sqrt(I.ref.L))
    log10K.MHL.over.M.HL.0 = log10K.ref.MHL.over.M.HL + (z.L*A*sqrt(I.ref.HL))/(1 + 1.5*sqrt(I.ref.HL))

    log10K.ML.over.M.L = log10K.ML.over.M.L.0 - (z.L*A*sqrt(I))/(1 + 1.5*sqrt(I))
    log10K.MHL.over.M.HL = log10K.MHL.over.M.HL.0 - (z.L*A*sqrt(I))/(1 + 1.5*sqrt(I))


    ####Calculate Kd.app in mM####

    Kd.app.ML.over.M.L = 1000/(a.L*(10^log10K.ML.over.M.L))
    Kd.app.MHL.over.M.HL = 1000/(a.HL*(10^log10K.MHL.over.M.HL))

    ####Print results####

    df.result = data.frame("Equillibrium" = c("ML/M.L", "MHL/M.HL"),
                           "mole fraction" = c(a.L, a.HL),
                           "Kd app mM" = c(Kd.app.ML.over.M.L, Kd.app.MHL.over.M.HL))

    print(df.result)

    ####Choose which Kd.app to use based on a####

    if (length(which(is.na(df.result$Kd.app.mM))) == 2){
      print(paste("Kd.app was not calculated because", metabolite, "binding to Mg2+ is not in the data set or is not predicted to significant affinity to Mg2+"))
      K.app = NA
    }else{
      if (length(which(is.na(df.result$Kd.app.mM))) == 1){
        K.app = df.result$Kd.app.mM[which(is.na(df.result$Kd.app.mM) != TRUE)]
      }else{
        if (a.L >= a.HL){
          K.app = Kd.app.ML.over.M.L
        }else{
          K.app = Kd.app.MHL.over.M.HL
        }
      }
    }
    ####Return####

    output = K.app

  }
}
