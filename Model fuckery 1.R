library(dplyr)

setwd("/Users/madeline/Desktop/Mideo lab/Parasite maturation")
#directory is set relative to the current working directory, which is set in the commander script

#####Parameters#####

  gamma<-0.005 #adjusts susceptibility based on age
  
  T.R.0<-0.97 #Total proportion of RBCs in circulation that are recipient
  T.D.0<-0.03 #Total proportion of RBCs in circulation that are donor
  
  B.R<-0.4 #Preference factor for recipient RBCs
  B.D<-0.6 #Preference factor for donor RBCs
  
  mat.D<-24 #Number of hours a parasite takes to develop and initiate cell bursting in a donor cell
  mat.R<-24 #Number of hours a parasite takes to develop and initiate cell bursting in a recipient cell
  
  #Concentration of RBCs in healthy host per ul
  R.norm<-8000000
  a<-0.05
  response<-"up"
  
  #Parasites produced from a burst RBC (donor or recipient)
  P.D<-8
  P.R<-8
  
  #Initial concentrations of RBCs
  R.0<-R.norm #if anaemia, can start lower; initial concentration of recipient RBCs
  DI.0<-(R.0*0.005)/(1-0.005) #initial concentration of infected donor cells
  DU.0<-(R.0*0.03)/(1-0.03) #initial concentration of uninfected donor cells

#####Functions#####

  #Function for calculating number of RBCs infected per age group

  
  infect.step<-function(df, B, T, gamma, P){
    
    Prop<-(T*B)/((T*B)+((1-T)*(1-B)))
    PD<-P*Prop
    
    for (i in 1:dim(df)[1]){
      temp.df<-as.data.frame(matrix(nrow=df[i,1], ncol=2))
      temp.df[,1]<-i
      temp.df[,2]<-exp((-i)*gamma)
      if (i == 1) {final.df<-temp.df} else {final.df<-bind_rows(final.df, temp.df)}
    }
    
    inf.temp<-sample(final.df$V1, size=PD, replace = FALSE, prob = final.df$V2)
    inf.table<-as.data.frame(table(inf.temp))
    names(inf.table)<-c("V1", "V2")
    inf.table[,1]<-as.integer(inf.table[,1])
    
    inf.df<-as.data.frame(c(1:dim(df)[1]))
    names(inf.df)<-"V1"
    inf.df[,1]<-as.integer(inf.df[,1])
    inf.ab<-full_join(inf.table, inf.df, by="V1")
    
    inf.ab[is.na(inf.ab)] <- 0
    
    return(inf.ab)
  }
  
  #Function determining the number of reticulocytes added in a time step
  retic.step<-function(res){
    retic.df<-as.data.frame(matrix(nrow=1, ncol=1))
    if (res == "up"){
      retic.df[1,1]<-((0.03*R.norm)/48)*(1+(sum(RI.adj, RU.adj, DI.adj, DU.adj)/R.norm)*a)
    } else if (res == "down"){
      retic.df[1,1]<-((0.03*R.norm)/48)*(1-(sum(RI.adj, RU.adj, DI.adj, DU.adj)/R.norm)*a)
    } else if (res == "none"){
      retic.df[1,1]<-((0.03*R.norm)/48)
    }
  }

#####Dateframes holding concentration of RBCs per age class (hours)#####
  
  #Dataframe holding proportion of uninfected recipient RBCs of age x hours (row number)
  RU.df<-as.data.frame(matrix(nrow=240, ncol=1))
  #Dataframe holding proportion of infected donor RBCs of age x hours (row number)
  DI.df<-as.data.frame(matrix(nrow=240, ncol=mat.D))
  #Dataframe holding proportion of infected donor RBCs of age x hours (row number)
  DU.df<-as.data.frame(matrix(nrow=240, ncol=1))
  #Dataframe holding proportion of infected donor RBCs of age x hours (row number)
  RI.df<-as.data.frame(matrix(nrow=240, ncol=mat.R))
  RI.df[,]<-0 # (initially zero recipient RBCs infected with donor parasite)
  
#####Dateframes holding summary values (one per generation)#####  
  RBC.sum<-as.data.frame(matrix(nrow=72, ncol=4))
  names(RBC.sum)<-c("RU","RI","DI","DU")
  
  P.sum<-as.data.frame(matrix(nrow=72, ncol=3))
  names(P.sum)<-c("Total","R","D")
  
  RI.age<-as.data.frame(matrix(nrow=72, ncol=240))
  DI.age<-as.data.frame(matrix(nrow=72, ncol=240))
  RU.age<-as.data.frame(matrix(nrow=72, ncol=240))
  DU.age<-as.data.frame(matrix(nrow=72, ncol=240))
  
  for (t in 1:72){

    if (t == 1){ #set up starting values
      
      #Total concentration of recipient RBCs
      R.con<-R.0
      
      #Reticulocytes
      RU.df[1:48,1]<-(0.03*R.con)/48 #fills in concentration for first 48 hour age classes
      DI.df[1:48,(1:mat.D)]<-(0.03*DI.0)/(48*mat.D)
      DU.df[1:48,1]<-(0.03*DU.0)/48
      
      #Normocytes
      RU.df[(49:dim(RU.df)[1]),1]<-(0.97*R.con)/(dim(RU.df)[1]-48) #fills in concentration for 49+ hour age classes
      DI.df[(49:dim(DI.df)[1]),(1:mat.D)]<-(0.97*DI.0)/((dim(DI.df)[1]-48)*mat.D)
      DU.df[(49:dim(DU.df)[1]),1]<-(0.97*DU.0)/(dim(DU.df)[1]-48)
      
      #RBC burst of infected donor cells (only donor cells bursting at t = 1)
      Burst.tot.D<-sum(DI.df[,mat.D])
      Burst.tot.R<-0
      Pr<-Burst.tot.D*P.D
      
      #Infect uninfected recipient and donor RBCs
      next.RI<-as.data.frame(infect.step(RU.df, B.R, T.R.0, gamma, Pr))
      next.DI<-as.data.frame(infect.step(DU.df, B.D, T.D.0, gamma, Pr))
      
      #Shifts columns so that parasites age
      DI.adj<-bind_cols(next.DI, DI.df[,1:(mat.D-1)])
      RI.adj<-bind_cols(next.RI, RI.df[,1:(mat.R-1)])
      
      #Removes newly infected RBCs from the uninfected to infected class
      DU.adj<-as.data.frame(as.matrix(DU.df[,1])-as.matrix(next.DI[,2]))
      RU.adj<-as.data.frame(as.matrix(RU.df[,1])-as.matrix(next.RI[,2]))
      
      #Set adjusted dataframes as next-generation dataframes
      RU.df<-RU.adj
      DU.df<-DU.adj
      DI.df<-DI.adj
      RI.df<-RI.adj
      
      #columns "RU","RI","DI","DU"  
      RBC.sum[t,1]<-sum(RU.df)
      RBC.sum[t,2]<-sum(RI.df)
      RBC.sum[t,3]<-sum(DI.df)
      RBC.sum[t,4]<-sum(DU.df)
      
      RI.age[t,]<-rowSums(RI.df)
      DU.age[t,]<-rowSums(DU.df)
      DI.age[t,]<-rowSums(DI.df)
      RU.age[t,]<-rowSums(RU.df)
      
      P.sum[t,1]<-(Burst.tot.D*P.D)+(Burst.tot.R*P.R)
      P.sum[t,2]<-(Burst.tot.R*P.R)
      P.sum[t,3]<-(Burst.tot.D*P.D)
      

    } else {
    
      #Remove oldest row of RBCs
      RU.adj<-as.data.frame(RU.df[(1:(dim(RU.df)[1])-1),])
      names(RU.adj)<-"V1"
      DU.adj<-as.data.frame(DU.df[(1:(dim(RU.df)[1])-1),])
      names(DU.adj)<-"V1"
      RI.adj<-as.data.frame(RI.df[(1:(dim(RU.df)[1])-1),])
      names(RI.adj)<-c(1:mat.R)
      DI.adj<-as.data.frame(DI.df[(1:(dim(RU.df)[1])-1),])
      names(DI.adj)<-c(1:mat.D)
      
      #Shifts rows so that RBCs age
      RU.new<-as.data.frame(retic.step(response))
      names(RU.new)<-"V1"
      DU.new<-as.data.frame(0)
      names(DU.new)<-"V1"
      RU.adj<-bind_rows(RU.new, RU.adj)
      DU.adj<-bind_rows(DU.new, DU.adj)
      
      #Add youngest age class to infected donor and recipient RBCs
      DI.new<-as.data.frame(matrix(nrow=1, ncol=mat.D))
      DI.new[1,]<-0
      names(DI.new)<-c(1:mat.D)
      RI.new<-as.data.frame(matrix(nrow=1, ncol=mat.R))
      RI.new[1,]<-0
      names(RI.new)<-c(1:mat.R)
      DI.adj<-bind_rows(DI.new, DI.adj)
      RI.adj<-bind_rows(RI.new, RI.adj)
      
      #RBC burst of infected cells
      Burst.tot.D<-sum(DI.df[,mat.D])
      Burst.tot.R<-sum(RI.df[,mat.R])
      
      #Number of parasites released
      Pr<-(Burst.tot.D*P.D)+(Burst.tot.R*P.R)
      
      #Total proportion of uninfected RBCs that are either recipient or donor
      T.R<-sum(RU.df)/(sum(RU.df)+sum(DU.df))
      T.D<-1-T.R
      
      #Infect uninfected recipient and donor RBCs
      next.RI<-as.data.frame(infect.step(RU.df, B.R, T.R, gamma, Pr))
      next.DI<-as.data.frame(infect.step(DU.df, B.D, T.D, gamma, Pr))
  
      #Shifts columns so that parasites age
      DI.adj<-bind_cols(next.DI, DI.df[,1:(mat.D-1)])
      RI.adj<-bind_cols(next.RI, RI.df[,1:(mat.R-1)])
    
      #Removes newly infected RBCs from the uninfected to infected class
      DU.adj<-as.data.frame(as.matrix(DU.df[,1])-as.matrix(next.DI[,2]))
      DU.adj[DU.adj<0]<-0
      RU.adj<-as.data.frame(as.matrix(RU.df[,1])-as.matrix(next.RI[,2]))
      RU.adj[RU.adj<0]<-0

      #Set adjusted dataframes as next-generation dataframes
      RU.df<-RU.adj
      DU.df<-DU.adj
      DI.df<-DI.adj
      RI.df<-RI.adj
      
      #columns "RU","RI","DI","DU"  
      RBC.sum[t,1]<-sum(RU.df)
      RBC.sum[t,2]<-sum(RI.df)
      RBC.sum[t,3]<-sum(DI.df)
      RBC.sum[t,4]<-sum(DU.df)
      
      RI.age[t,]<-rowSums(RI.df)
      DU.age[t,]<-rowSums(DU.df)
      DI.age[t,]<-rowSums(DI.df)
      RU.age[t,]<-rowSums(RU.df)
      
      P.sum[t,1]<-(Burst.tot.D*P.D)+(Burst.tot.R*P.R)
      P.sum[t,2]<-(Burst.tot.R*P.R)
      P.sum[t,3]<-(Burst.tot.D*P.D)
      
      
    } #end if/else statement
  
  } #end over t loop
  
  
  
  write.csv(RBC.sum, "RBC.sum.csv", row.names=FALSE)
  write.csv(P.sum, "P.sum.csv", row.names=FALSE)
  write.csv(RI.age, "RI.age.csv", row.names=FALSE)
  write.csv(RU.age, "RU.age.csv", row.names=FALSE)
  write.csv(DI.age, "DI.age.csv", row.names=FALSE)
  write.csv(DU.age, "DU.age.csv", row.names=FALSE)
  

  
