##Full simulation study, setting 
N=58007
p_M=0.5164
p_F=1-p_M


results=matrix(nrow=1000,ncol=100)

C1=rep(1,N) #forces everyone to have at least 1 kid
C2=rep(1,N) #forces everyone to have at least 2 kids

s_count=1
set.seed(1919)
for (s in seq(from=1,to=1.5,length.out=100) ){
  for (i in 1:1000){
    S1=rbinom(n=N,size=1,p=p_M) #sex random
    S2=rep(NA,N)
    S2[C2==1]=rbinom(n=sum(C2==1),size=1,p=p_M) #sex random
    C3=rep(0,N)
    C3[C2==1]=rbinom(n=sum(C2==1),size=1,p=(0.354*(1+(s-1)*(S1[C2==1]==S2[C2==1]))) ) #s varies from 1-1.5 #0.354 is probability of 3rd if first two different, from the data
    S3=rep(NA,N)
    S3[C3==1]=rbinom(n=sum(C3==1,na.rm=T),size=1,p=p_M) #sex random
    
    #just save association
    #Probability of all 3 same sex conditional on having 3
    results[i,s_count]=sum(S1==1&S2==1&S3==1|S1==0&S2==0&S3==0,na.rm=T)/sum(C3==1)*100
    print(c(s_count,i))
  }
  s_count=s_count+1
}
write.table(results,file="sim_results2.txt",append=FALSE)

#if need to load from file
results=as.matrix(read.table("sim_results2.txt"))


library(ggplot2)
library("latex2exp")

# Example data with pre-calculated mean and CI
summary_data <- data.frame(
  s = seq(from=1,to=1.5,length.out=100),
  mean_val = colMeans(results)/100, #dim 100
  ci_lower = apply(results/100,2,FUN=function(x)quantile(x,0.05) ),
  ci_upper = apply(results/100,2,FUN=function(x)quantile(x,0.95) )
)
single_point_df <- data.frame(x = 1.205, y = 0.275 ) #I used  data to calculate (Excel file)


#check stuff against Judith's dataset (should be equal to my data-based calculations)
#Corollary 3.2: probability of having three kids that instead uses p_M and p_F as in Judith's derivation. It uses law of total probability
cor32=function(s){(p_M^3+p_F^3)/(2*p_M*p_F*p_d/p_s(s)+p_M^2+p_F^2)}

p_d=0.354
p_s=function(s){p_d*s} #s will be the inflation factor greater than 1 #probability of at least third kid if first two are same

#just adding the new line
expectedval <- data.frame(
  s = seq(from=1,to=1.5,length.out=100),
  vals = cor32(s)
)

#install.packages("viridis")
library(viridis)

ggplot(summary_data, aes(x = s, y = mean_val)) +
  theme_light()+
  #theme(panel.grid.major = element_line(color = "lightskyblue1"),panel.grid.minor = element_line(color = "lightskyblue1"))+ # Change gridline color to blue
  ggtitle("Probability of having first 3 same-sex children", subtitle = TeX("conditional on having 3 or more children given $p_S/p_D$"))+ xlab(TeX("$p_S/p_D$"))+ylab("Probability")+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = 1), alpha = 0.6, fill="#440154FF") + # Adds the confidence band, alpha is apacity
  geom_smooth(se = FALSE, span=0.4, linetype="solid", size=1, colour="#22A884FF") + # Connects the means with a loess line, span is % of data points used for local line
  geom_hline(yintercept = (p_M^3+p_F^3), linetype="dashed",size=1.5, colour="#414487FF")+ #theoretical binomial line
  geom_point(data = single_point_df, aes(x = x, y = y), size = 5, shape = 18)+#, color = "grey45")+
  annotate("text", x = 1.4, y = (p_M^3+p_F^3)+0.005, label = "No inflation", size = 4, color = "#414487FF")+
  annotate("text", x = 1.33, y = 0.273, label = "Wang et al data", size = 4, color = "black")+
  annotate("label", x = 1.1, y = 0.3, label = "    Theoretical prob\n    Running average\n     95% density band", color = "black", size = 4)+
  annotate("rect", xmin=c(1.01), xmax=c(1.025), ymin=c(0.2934) , ymax=c(0.2960), alpha=0.6, color="#440154FF", fill="#440154FF")+
  annotate("rect", xmin=c(1.01), xmax=c(1.025), ymin=c(0.2998) , ymax=c(0.3001), alpha=1, color="#22A884FF", fill="#22A884FF")+
  annotate("rect", xmin=c(1.012), xmax=c(1.014), ymin=c(0.3043) , ymax=c(0.3045), alpha=1, color="#FDE725FF", fill="#FDE725FF")+
  annotate("rect", xmin=c(1.019), xmax=c(1.021), ymin=c(0.3043) , ymax=c(0.3045), alpha=1, color="#FDE725FF", fill="#FDE725FF")+
stat_function(fun = cor32, color = "#FDE725FF", linetype = "dotted", size = 1)
#annotate("label", x = 1.2, y = 0.29, label = "running average (loess curve)", color = "black", size = 4)
#geom_vline(xintercept =1.2,colour="red3",linetype="dotdash",size=1.5)
# geom_hline(yintercept =) #should probably be data-informed theoretical line

