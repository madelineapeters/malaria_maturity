library(ggplot2)

p.plot<-ggplot(data=P.sum)+geom_line(aes(x=c(1:72),y=Total))+geom_line(aes(x=c(1:72), y=R))+geom_line(aes(x=c(1:72), y=D))
print(p.plot)

RBC.plot<-ggplot(data=RBC.sum[1:15,])+geom_line(aes(x=c(1:15),y=RU), color="red")+
  geom_line(aes(x=c(1:15), y=RI), color="orange")+
  geom_line(aes(x=c(1:15), y=DU), color="darkblue")+
  geom_line(aes(x=c(1:15), y=DI), color="blue")+
  #geom_line(data=P.sum[1:15,], aes(x=c(1:15), y=Total))+
  theme_bw()
print(RBC.plot)

ratio.df<-as.data.frame(matrix(nrow=72, ncol=2))
ratio.df[,1]<-RBC.sum$RI/RBC.sum$RU
