# MetaboMgITC

## Metabolite binding to divalent Mg cations at biologically relevant pH and ionic strength

### Authors: Jacob P. Sieg, Ryota Yamagami, and Philip C. Bevilacqua

## 1. What is MetaboMgITC?

#### RNA regulates myriad cellular events such as transcription, translation, and splicing. To perform these essential functions, RNA folds into complex tertiary structures where a negatively charged ribose-phosphate backbone interacts with metal ions. Magnesium, the most abundant divalent metal ion in cells, neutralizes the backbone thereby playing essential roles in RNA folding and function. In the cell, some Mg2+ ions are chelated by metabolites. Recently, we have demonstrated that this metabolite chelated magnesium (MCM) pool can effect RNA folding and function, depending on the strength of the interaction between Mg2+ and metabolites. The binding of Mg2+ to metabolites is highly dependent on pH and ionic strength, which changes the population of Mg2+ binding competent protonation states for a metabolite ligand and changes the activity of ions in aqueose solution respectivly, further complicating the actual Mg2+ status in cells. In an effort towards gaining insight into the role of biologically chelated magnesium ions in RNA chemistry and biology, we have:

#### (1) Currated a database of 292 absolute metabolite concentrations from E. coli, Yeast, and mouse iMBK cells.
#### (2) Currated a database of pKa's, Mg2+ binding constants, ionic strengths, and temperatures for 120 metabolite ligands. This data, 582 constants in total, is used to calculate apparant metabolite/Mg2+ dissacociation constants at specific pH's and ionic strengths.
#### (3) R functions to provide a user with a relatively facile way to acess and interperate the data. Most importantly, we provide "Kd.app.calc", which enables the user to calculate apparent metabolite/Mg2+ dissacociation constants at aany pH and ionic strength.
#### (4) Isothermal titration calorimetry (ITC) data analysis tools to experimentally determine apparent metabolite/Mg2+ dissacociation constants at a pH and ionic strength.

#### An example implementation of MetaboMgITC is dementrated below, where MCM levels are projected for E. coli, Yeast, and mouse iMBK cells at pH 7.5 and ionic strength = 0.15 M.

![Final](https://user-images.githubusercontent.com/63312483/120546207-96225600-c3bd-11eb-941d-49a7ace89fa9.png)

# 2. Using MetaboMgITC on your own console

## Dependencies

#### R version version 4.1.0
#### devtools (R package)
#### pracma (R package)
#### dplyr (R package)
#### ggplot2 (R package)

## Recomended packages and programs

#### R studio
#### tidyverse (R package)

## Video Tutorials

#### 1.) Reconstituting MetaboMgITC on your own console

https://youtu.be/OvVh7j-uIcE

#### 2.) Using Kd.app.calc

https://youtu.be/ViZcHrtGpZs 

#### 3.) ITC data analysis

https://youtu.be/hlL0VxYa2vw

# 3. Methods

## 3.1 Apparant dissacociation constant approximation

#### The apparent disassociation constant (K<sub>D</sub>') for a metal ion  binding to a metabolite at a pH and ionic strength is:

K<sub>D</sub>' = [M][L]/[ML]  Equation 1

#### Where [M] is the concentration of metal ions. [L] and [ML] are the sum of the concentration of all protonation states for the metabolite and the metabolite magnesium complex respectively. The K<sub>D</sub>' for each ligand at a given ionic strength and pH was calculated using absolute metal ion binding constants for a given ligand protonation state and protonation constants (pKa’s) from the literature using equation 2.

K'<sub>D</sub> = 1/(a<sub>i</sub>10<sup>log<sub>10</sub>K'</sup>)  Equation 2

#### Let’s first consider a<sub>i</sub>. The fraction of fully deprotonated (L) and singly protonated (HL) ligand species (i.e. the two species most likely to bind metal ion) were determined from protonation constants (pKa’s) from the literature using equations 3 and 4. For most metabolites, the metal ion binding-competent state is completely deprotonated. The major exceptions are L-Glutamic acid, L-Aspartic acid, and L-Arginine, which have affinity for divalent magnesium cations in the singly protonated state. The mole fraction of ligand in protonation states L and HL was calculated using equations 3 and 4 respectively.<sup>2</sup>

a<sub>L</sub> = 1/(1 + 10<sup>pKa<sub>HL</sub>-pH</sup> + 10<sup>pKa<sub>HL</sub> + pKa<sub>H<sub>2</sub>L</sub> - 2pH</sup>) Equation 3

a<sub>HL</sub> =  10<sup>pKa<sub>HL</sub>-pH</sup>/(1 + 10<sup>pKa<sub>HL</sub>-pH</sup> + 10<sup>pKa<sub>HL</sub> + pKa<sub>H<sub>2</sub>L</sub> - 2pH</sup>) Equation 4

#### Now let’s consider K'.  This is the ionic strength-corrected binding constant for a ligand in a single protonation state<sup>2</sup>.  It was corrected using the simplified Debye-Hückel equation<sup>3</sup>:

log<sub>10</sub>K' = log<sub>10</sub>K<sup>o</sup> - zAI<sup>0.5</sup>/(1 + 1.5I<sup>0.5</sup>) Equation 5

#### where K<sup>o</sup> is the binding constant at an ionic strength of zero and A is an empirical constant of 0.524 for aqueous solutions at 25 °C. The constant z is calculated using equation 6, where z<sub>L</sub> is the charge of the ligand in the relevant protonation state, z<sub>M</sub> is the charge of the free metal ion, and z<sub>ML</sub> is the charge of the metal ion-ligand complex.

z = z<sub>L</sub><sup>2</sup> + z<sub>M</sub><sup>2</sup> - z<sub>ML</sub><sup>2</sup> Equation 6

#### We compiled a database of protonation constants and divalent magnesium binding constants for cellular metabolites (194 constants).

#### Protonation constants were pulled from Arthur E. Martell’s and Robert M. Smith’s body of work known as Critical Stability Constants<sup>4–9</sup>, unless otherwise stated (Supplementary Table 3).<sup>10–13</sup> Protonation constants were corrected for ionic strength using the directions in the introduction of Critical Stability Constants Volume 1. Critical Stability Constants Volume 1, Volume 2, Volume 3, Volume 5, and Volume 6 were accessed from the Internet Archive (https://archive.org) and NIST46 Critically Selected Stability Constants of Metal Complexes. Version 8.0 was accessed from the National Institute of Standards and Technology (https://www.nist.gov/srd/nist46) and ran on a computer with a Windows XP operating system. Critical Stability Constants does not specify an ionic strength correction so we used equation 5. The source for all constants are specified by "Binding_constant_concentration_data/210525_Metaboites_binding_Mg_thermodynamics.csv".

#### Data and calculations were coded into a bare-bones R based program called “Kd.app.calc” and used to calculate the K<sub>D</sub>'’s at any ionic strength and pH. Usage information and help can be found in Video Tutorial 2. Source code is housed in "R/Thermodynamic_equations.R".

## 3.2 Calculating MCM levels in cells

#### MCM concentrations can be calculated using the free Mg2+ concentration in cells, about 0.5 to 3 mM, and the K<sub>D</sub>' by rearranging equation 1. In the example analysis shown here, we use 2 mM free Mg2+ for E. coli cells and 0.5 mM free Mg2+ for Yeast and iMBK cells.

[ML] = [L]<sub>Total</sub>*[M]/(K<sub>D</sub> + [M]) Equation 7

#### This assumes one-to-one binding between all metabolites and Mg ions, which is reasonable given that there is little evidence of multivalent metabolite interections with Mg in the liturature, where multivalent interactions are insoluble and occur at molar concentrations, well above the mM concentrations found in vivo.

## 3.3 ITC data collection and analysis

#### Isothermal Titration Calorimetry. Buffers were prepared by dissolving high purity (>99.99%) salts in purified (18 mΩ) water in volumetric flasks to the final concentrations described in Supplementary Table 1. A standard 1 M (±0.01) magnesium chloride solution from Millipore Sigma was diluted to ensure accurate magnesium chloride concentrations. Samples were degassed at 25°C using a ThermoVac (MicroCal, LLC) before loading into a VP-ITC MicroCalorimeter (MicroCal, LLC) according to the manufactures recommendations with pure (18 mΩ) water in the reference cell. Titration was performed at 25°C with a 10 µcal/sec reference power and a stirring speed of 310 rpm. 29 total injections (one 2 μL injection followed by twenty eight 10 μL injections) were performed at an injection rate of 0.5 μL/sec with a 200 sec space following each injection. The first injection was not included in subsequent analysis.

#### Data were analyzed in R (version 4.1.0). Briefly, raw ITC file were parsed using a custom function. Heats of injections were determined by integrating the raw differential power curves with the polyarea function in pracma (https://CRAN.R-project.org/package=pracma). Curves were fit to a 1:1 Wisman Isotherm1 using the nls non-linear regression function in base R 4.1.0 to determine apparent association and disassociation constants. Usage information and help can be found in Video Tutorial 3.

# 4. References

##### (1) 	Tellinghuisen, J. Isothermal Titration Calorimetry at Very Low c. Anal. Biochem. 2008, 373 (2), 395–397. https://doi.org/10.1016/j.ab.2007.08.039.
##### (2) 	Mattocks, J. A.; Tirsch, J. L.; Cotruvo, J. A. Chapter Two - Determination of Affinities of Lanthanide-Binding Proteins Using Chelator-Buffered Titrations. In Methods in Enzymology; Cotruvo, J. A., Ed.; Rare-Earth Element Biochemistry: Characterization and Applications of Lanthanide-Binding Biomolecules; Academic Press, 2021; Vol. 651, pp 23–61. https://doi.org/10.1016/bs.mie.2021.01.044.
##### (3) 	Scatchard, George. Concentrated Solutions of Strong Electrolytes. Chem. Rev. 1936, 19 (3), 309–327. https://doi.org/10.1021/cr60064a008.
##### (4) 	Martell, A. E.; Smith, R. M. Critical Stability Constants; NewYork ; London : Plenum Press, 1974; Vol. V1.
##### (5) 	Martell, A. E.; Smith, R. M. Critical Stability Constants; New York, Plenum Press, 1974; Vol. V2.
##### (6) 	Martell, A. E.; Smith, R. M. Critical Stability Constants; New York, Plenum Press, 1974; Vol. V3.
##### (7) 	Martell, A. E.; Smith, R. M. Critical Stability Constants.; New York ; London : Plenum, 1982; Vol. V5.
##### (8) 	Smith, R. M. (Robert M.; Martell, A. E. Critical Stability Constants; New York ; London : Plenum, 1989; Vol. V6.
##### (9) 	Martell, A. E.; Smith, R. M. NIST46 Critically Selected Stability Constants of Metal Complexes https://www.nist.gov/srd/nist46 (accessed May 25, 2021).
##### (10) 	Berthon, G. Critical evaluation of the stability constants of metal complexes of amino acids with polar side chains (Technical Report). Pure Appl. Chem. 1995, 67 ##### (7), 1117–1240. https://doi.org/10.1351/pac199567071117.
##### (11) 	Berthon, G. Speciation Studies in Relation to Magnesium Bioavailability. Formation of Mg(I1) Complexes with Glutamate, Aspartate, Glycinate, Lactate, Pyroglutamate, Pyridoxine and Citrate, and Appraisal of Their Potential Significance Towards Magnesium Gastrointestinal Absorption. Inorganica Chim. Acta 1987, 135 (3), 179–181.
##### (12) 	Albert, A. Quantitative Studies of the Avidity of Naturally Occurring Substances for Trace Metals. 2. Amino-Acids Having Three Ionizing Groups. Biochem. J. 1952, 50 (5), 690–697. https://doi.org/10.1042/bj0500690.
##### (13) 	Sovago, I.; Kiss, T.; Gergely, A. Critical Survey of the Stability Constants of Complexes of Aliphatic Amino Acids (Technical Report). Pure Appl. Chem. 1993, 65 (5), 1029–1080. https://doi.org/10.1351/pac199365051029.
##### (14) 	Yamagami, R.; Bingaman, J. L.; Frankel, E. A.; Bevilacqua, P. C. Cellular Conditions of Weakly Chelated Magnesium Ions Strongly Promote RNA Stability and Catalysis. Nat. Commun. 2018, 9 (1), 2149. https://doi.org/10.1038/s41467-018-04415-1.
##### (15) 	Yamagami, R.; Huang, R.; Bevilacqua, P. C. Cellular Concentrations of Nucleotide Diphosphate-Chelated Magnesium Ions Accelerate Catalysis by RNA and DNA Enzymes. Biochemistry 2019, 58 (38), 3971–3979. https://doi.org/10.1021/acs.biochem.9b00578.




