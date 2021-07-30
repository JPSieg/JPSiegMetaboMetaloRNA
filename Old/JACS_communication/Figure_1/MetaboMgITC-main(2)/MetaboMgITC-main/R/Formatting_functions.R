#'Imports a VP-ITC instument ".itc" formatted data file into R
#'
#'Opens a MicroCal VP-ITC instrument ".itc" formatted data file into R. Parses raw text into convenient data
#'frames. Calculates reagent concentrations. Integrates raw differential power curves to determine the heat
#'of each injection using the polyarea function in pracma.
#'
#'@param file.concentrations Vector containing paths to an ITC file then reagent concentration ordered c(Path, syringe X concentration, syringe M concentration, syringe C concentration, cell X concentration, cell M concentration, cell C concentration). Enter a file path to use the concentration embedded in the raw .itc file.
#'@param print.peak.integration.graph Print an aesthetic depiction of the peak integration for each injection. Options = TRUE or FALSE. Default = FALSE. If TRUE will print a png in your working directory.
#'@return A list of two data frames that can be passed to data analysis functions like "MetaboMgITC". The first dataframe contains the raw differential power (microcalories per second) as a function of time (second). The second data frame contains the heat of each injection (kilocalorie per mole injectant) and reagent concentration (molar). Volumes are reported as Liters. M is the concentration of the reagent in the cell. X and dX are the concentration of the ligand in the cell and the moles of ligand added in each injection respectively. dQ and dQ.dX are the heat produced by each injection and the heat produced by each injection divided by the moles added by each injection.
#' @export
read.itc =function(file.concentrations,
                   print.peak.integration.graph = FALSE,
                   integration.folder = NA){
  ####Parses file.concentrations to find path####

  if (length(file.concentrations) == 1){
    path.to.itc.file = file.concentrations
  }else{
    path.to.itc.file = file.concentrations[1]
  }

  ####Open a connection to a ITC file####


  con = file(path.to.itc.file, "r")

  Temp.Found = FALSE
  reagent.lines = c()
  ####Parses ITC data file into convenient R data frames####
  while ( TRUE ) {
    line = readLines(con, n = 1)
    #print(line)
    if ( length(line) == 0 ) {
      break
    }
    #Finds the temperature on the second line of the ITC file
    if (line != "$ITC"){
      if (Temp.Found == FALSE){
        Temp = as.numeric(as.character(strsplit(line, split = " ")[[1]][2]))
        Temp.Found = TRUE
      }
    }
    #Pulls out injection information
    if (strsplit(line, split = ":")[[1]][1] == "$ADCGainCode"){
      inj.V = c()
      inj.dur = c()
      inj.delay =c()
      while(strsplit(line, split = " ")[[1]][1] != "#"){
        line.vector = strsplit(line, split = " ")[[1]]
        if (length(line.vector) > 1){
          if (line.vector[1] != "#"){
            if (line.vector[1] != "$ADCGainCode:"){
              line.vector = line.vector[which(line.vector != "$")]
              line.vector = line.vector[which(line.vector != ",")]
              inj.V = c(inj.V, as.numeric(as.character(line.vector[1])))
              inj.dur = c(inj.dur, as.numeric(as.character(line.vector[2])))
              inj.delay =c(inj.delay, as.numeric(as.character(line.vector[3])))
              #print(line.vector)
            }
          }
        }
        line = readLines(con, n = 1)
      }
    }
    #Pulls out volumes and reagent concentrations

    if (strsplit(line, split = " ")[[1]][1] == "#"){
      #print(line)
      reagent.lines = c(reagent.lines, strsplit(line, split = " ")[[1]][2])
      #print(reagent.lines)
    }

    #Pulls out delta power data
    if (strsplit(line, split = "")[[1]][1] == "@"){
      inj.n = c()
      Time.sec = c()
      dP = c()
      Temperature = c()
      while(TRUE){
        if (strsplit(line, split = "")[[1]][1] == "@"){
          #print(line)
          inj.N = as.numeric(as.character(gsub("@", "", strsplit(line, split = ",")[[1]][1])))
          #print(inj.N)
          line = readLines(con, n = 1)
        }else{
          #print(line)
          line.vector = strsplit(line, split = ",")[[1]]
          inj.n = c(inj.n, inj.N)
          Time.sec = c(Time.sec, as.numeric(as.character(line.vector[1])))
          dP = c(dP, as.numeric(as.character(line.vector[2])))
          Temperature = c(Temperature, as.numeric(as.character(line.vector[3])))
          line = readLines(con, n = 1)
        }
        if ( length(line) == 0 ) {
          break
        }
      }
    }
  }
  ####Close the connection to the ITC data file####

  close(con)

  ####Parse reagent concentration data####

  if(length(file.concentrations) == 1){
    V0 = as.numeric(as.character(reagent.lines[4]))
    M0 = as.numeric(as.character(reagent.lines[3]))/1000
    X0 = as.numeric(as.character(reagent.lines[2]))/1000
  }else{
    V0 = as.numeric(as.character(reagent.lines[4]))
    X0.syringe = as.numeric(as.character(file.concentrations[2]))/1000
    M0.syringe = as.numeric(as.character(file.concentrations[3]))/1000
    C0.syringe = as.numeric(as.character(file.concentrations[4]))/1000
    X0.cell = as.numeric(as.character(file.concentrations[5]))/1000
    M0.cell = as.numeric(as.character(file.concentrations[6]))/1000
    C0.cell = as.numeric(as.character(file.concentrations[7]))/1000
  }

  ####Make a data frame for the raw differential power####

  df.dP = data.frame(inj.n, Time.sec, dP, Temperature)

  ####make a data frame that summarizes each injection####

  N = 1:length(inj.V)
  df.inj = data.frame(N, inj.V, inj.dur, inj.delay)

  ####Determine the concentration of reagents after each injection####

  if(length(file.concentrations) == 1){
    V = c()
    X =c()
    M = c()
    dX = c()

    for (i in df.inj$N){
      V[i] = V0 + sum(df.inj$inj.V[1:i])/1000
      X[i] = (sum(df.inj$inj.V[1:i])/1000)*X0/V[i]
      M[i] = V0*M0/V[i]
    }

    df.inj$V = V
    df.inj$X = X
    df.inj$M = M
  }else{
    V = c()
    X =c()
    M = c()
    C = c()
    dX = c()

    for (i in df.inj$N){
      V[i] = V0 + sum(df.inj$inj.V[1:i])/1000
      X[i] = (sum(df.inj$inj.V[1:i])/1000)*X0.syringe/V[i] + V0*X0.cell/V[i]
      M[i] = (sum(df.inj$inj.V[1:i])/1000)*M0.syringe/V[i] + V0*M0.cell/V[i]
      C[i] = (sum(df.inj$inj.V[1:i])/1000)*C0.syringe/V[i] + V0*C0.cell/V[i]
    }

    df.inj$V = V
    df.inj$X = X
    df.inj$M = M
    df.inj$C = C
  }


  ####Determine the heat produced by each injection####

  dQ = c()
  for (i in 1:length(df.inj$N)){
    df <- subset(df.dP, df.dP$inj.n == df.inj$N[i])
    #pracma::polyarea(c(0,1,1,0), c(0, 0, 1, 1))
    #geometry::polyarea(c(0,1,1,0), c(0, 0, 1, 1))
    #geometry::polyarea(df$Time.sec, df$dP)
    dQ[i] = -pracma::polyarea(df$Time.sec, df$dP)
  }

  if(length(file.concentrations) == 1){
    df.inj$dQ = dQ
    df.inj$dX = (df.inj$inj.V/10^6)*(X0)
    df.inj$dQ.dX = df.inj$dQ/df.inj$dX/(10^9)
  }else{
    df.inj$dQ = dQ
    df.inj$dX = (df.inj$inj.V/10^6)*(X0.syringe)
    df.inj$dQ.dX = df.inj$dQ/df.inj$dX/(10^9)
  }


  ####Print a graph in order to check peak integration####

  if (print.peak.integration.graph){
    if (length(is.na(integration.folder)) != 1){
      ggplot2::ggplot(df.dP, ggplot2::aes(x = Time.sec, y = dP, color = factor(inj.n), fill = factor(inj.n)))+
        ggplot2::geom_line() +
        ggplot2::geom_polygon() +
        ggplot2::scale_fill_manual(name = "Injection", values = viridis::viridis(length(unique(df.dP$inj.n)))) +
        ggplot2::scale_color_manual(name = "Injection", values = viridis::viridis(length(unique(df.dP$inj.n)))) +
        ggplot2::theme_classic() +
        ggplot2::xlab("Time (sec)") +
        ggplot2::ylab("")
      ggplot2::ggsave(paste(integration.folder, "/", path_to_itc_file,"_integration.png", sep = ""))
    }
  }

  ####Make a list containing the data frames####
  output = list(df.dP, df.inj)
  names(output) = c("dP", "inj")
  print(path.to.itc.file)
  print(names(output)[1])
  print(head(output[[1]]))
  print(names(output)[2])
  print(head(output[[2]]))
  output = output
}


