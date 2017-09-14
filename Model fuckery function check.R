infect.step<-function(df, B, T, gamma, P){
  
  df<-DU.df
  B<-0.6
  T<-0.03
  P<-Pr
  Prop<-(T*B)/((T*B)+((1-T)*(1-B)))
  inf.df<-as.data.frame(sapply(X=(1:dim(df)[1]), FUN=function(X) {
    as.numeric(exp(-1*(X)*gamma))
  } ))
  inf.df<-bind_cols(df, inf.df)
  names(inf.df)<-c("V2", "V3")
  inf.df$V2<-(inf.df$V2)/sum(inf.df$V2, na.rm=TRUE)
  inf.df$V3<-(inf.df$V3)/sum(inf.df$V3, na.rm=TRUE)
  inf.corr<-transmute(inf.df, V1=V2*V3)
  inf.corr<-inf.corr/sum(inf.corr, na.rm=TRUE)
  
  inf.ab<-as.data.frame(inf.corr*P*Prop)
  
  return(inf.ab)
  
}

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