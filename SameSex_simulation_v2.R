###Simulation study for the article
###"Sex at birth could well be a biological coin toss.... Beware of conditioning 
###on post-baseline information" 
###by Judith Lok and Mireille Schnitzer

#Version date 2026-07-09

#This article is a response to "Is sex at birth a biological coin toss? Insights 
#from a longitudinal and GWAS analysis" by Wang et al. 2025, Science Advances 11 (29), 
#eadu7402.

#This simulation aims to show that, even if the sex of any offspring is 
#completely random with a constant probability of being male, if a family is
#more likely to have a third child if the first two have the same sex, the 
#probability of having all same-sex children conditional on having at least
#three children, will be greater than what would be predicted by a binomial
#distribution. Specifically, if p_M is the probability of having a male child,
#the probability of having all male or all female children conditional on
#having at least three children will be greater than p_M^3+(1-p_M)^3.

#Here, we explore a setting that is simplified in the following ways: all
#families have at least 2 kids and there is no effect of the sexes of the first
#two kids on whether to have another kid *except* whether the first two kids are 
#of the same sex. (I.e. in the terms of Lok and Schnitzer, we only consider the 
#"second correction factor" as defined in Section 3.)

##load packages
library(ggplot2)
library("latex2exp")
library("viridis")

##Constants and initialization

N=58007 #total sample size from Wang et al. (2025)
p_M=0.5164 #probability of male offspring, value from Wang et al.
p_F=1-p_M #probability of female offspring

#The correction factor s=p_s/p_D is the probability of at least having a third 
#kid if the first two are of same sex

#initialize matrix of results, over 1000 simulations, and 100 values of the 
#correction factor (s)
results=matrix(nrow=1000,ncol=100)

#Ck, k=1,2,3 is an indicator of having k kids
#Sk, k=1,2,3 is an indicator of the sex of the kth kid, NA if kid does not exist

C1=rep(1,N) #forces everyone to have at least 1 kid
C2=rep(1,N) #forces everyone to have at least 2 kids

s_count=1 #index of correction factor values

##Simulation loop
set.seed(1919)
for (s in seq(from=1,to=1.5,length.out=100) ){ #loops over correction factor values (1 (no inflation) to 1.5)
  for (i in 1:1000){ #Monte Carlo simulation loop
    S1=rbinom(n=N,size=1,p=p_M) #sex of 1st kid (1 male, 0 female) -- completely random
    S2=rep(NA,N)
    S2[C2==1]=rbinom(n=sum(C2==1),size=1,p=p_M) #sex of 2nd kid -- completely random
    
    #C3 is whether the couple has a 3rd kid, drawn from a binomial distribution
    #If s>1, the probability will be conditional on whether the first two kids have the same sex
    C3=rep(0,N) #initialize indicator of 3rd kid
    C3[C2==1]=rbinom(n=sum(C2==1),size=1,p=(0.354*(1+(s-1)*(S1[C2==1]==S2[C2==1]))) ) #0.354 is probability of having 3rd kid if first two are of different sex, from Wang et al data
    
    S3=rep(NA,N)
    S3[C3==1]=rbinom(n=sum(C3==1,na.rm=T),size=1,p=p_M) #sex of 3rd kid -- completely random
    
    #Save the probability of all 3 kids having the same sex conditional on having 3 kids
    results[i,s_count]=sum(S1==1&S2==1&S3==1|S1==0&S2==0&S3==0,na.rm=T)/sum(C3==1)
    print(c(s_count,i)) #prints the indices
  }
  s_count=s_count+1 #increment index of correction factor
}
write.table(results,file="sim_results1.txt",append=FALSE) #save results to file


##Summarize results

#if loading a previously generated simulation results data file (in which case, skip the section 
#"Simulation loop"), run the following:
#results=as.matrix(read.table("sim_results1.txt"))

# Summarize the simulation results by taking mean and CI for each value of s
summary_data <- data.frame(
  s = seq(from=1,to=1.5,length.out=100), #iterate over values of s
  mean_val = colMeans(results), #dim 100
  ci_lower = apply(results,2,FUN=function(x)quantile(x,0.025) ), #lower confidence intervals
  ci_upper = apply(results,2,FUN=function(x)quantile(x,0.975) ) #upper confidence intervals
)

#From Lok and Schnitzer, we used the Wang et al data to estimate the correction
#factor: 
#p_s = the probability that a family of 2 with same-sex children will have a 
#third child
#p_d = the probability that a family of 2 with mixed-sex children will have a 
#third child
#p_s=0.426 and p_d=0.354 (see Section 2 of Lok and Schnitzer); 
#s=p_s/p_d=1.203 is the inflation factor (x-axis).
#Also from Lok and Schnitzer (Section 4), the estimated probability of all 
#same-sex conditional on at least three kids is estimated as 0.2743 using the 
#Wang et al data (y-axis).
single_point_df <- data.frame(x = 1.203, y = 0.2751 ) #this is the estimate from the Wang et al data, to contrast with the theoretical values


#Theoretical values from mathematical derivation -- should correspond to the 
#simulation results
#Corollary 3.2: probability of having three kids based on values of p_M and p_F. 
cor32=function(s){(p_M^3+p_F^3)/(2*p_M*p_F*p_d/p_s(s)+p_M^2+p_F^2)}

p_d=0.354 #0.354 is probability of having a 3rd kid if first two are of different sex, from the Wang et al data
p_s=function(s){p_d*s} #probability of at least third kid if first two are of the same sex (s is the correction factor)

#creating the data to later plot the line with the theoretical values
expectedval <- data.frame(
  s = seq(from=1,to=1.5,length.out=100),
  vals = cor32(s)
)

##Plot the simulation results (mean and confidence band), theoretical results, 
##and data point
ggplot(summary_data, aes(x = s, y = mean_val)) +
  theme_light()+theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))+
  #theme(panel.grid.major = element_line(color = "lightskyblue1"),panel.grid.minor = element_line(color = "lightskyblue1"))+ # Change gridline color to blue
  ggtitle("Probability of having first 3 same-sex children", subtitle = TeX("conditional on having 3 or more children given $p_S/p_D$"))+ xlab(TeX("$p_S/p_D$"))+ylab("Probability")+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = 1), alpha = 0.6, fill="#440154FF") + # Adds the confidence band, alpha is opacity
  geom_smooth(se = FALSE, span=0.4, linetype="solid", size=1, colour="#22A884FF") + # Connects the means with a loess line, span is % of data points used for local line
  geom_hline(yintercept = (p_M^3+p_F^3), linetype="dashed",size=1.5, colour="#414487FF")+ #theoretical binomial line
  geom_point(data = single_point_df, aes(x = x, y = y), size = 5, shape = 18)+#, color = "grey45")+ #values derived from Wang et al data
  annotate("text", x = 1.4, y = (p_M^3+p_F^3)+0.005, label = "No inflation", size = 4, color = "#414487FF")+ 
  annotate("text", x = 1.33, y = 0.273, label = "Wang et al. (2025) data", size = 4, color = "black")+
  annotate("label", x = 1.1, y = 0.3, label = "    Theoretical prob\n    Running average\n     95% density band", color = "black", size = 4)+
  annotate("rect", xmin=c(1.01), xmax=c(1.025), ymin=c(0.2934) , ymax=c(0.2960), alpha=0.6, color="#440154FF", fill="#440154FF")+
  annotate("rect", xmin=c(1.01), xmax=c(1.025), ymin=c(0.2998) , ymax=c(0.3001), alpha=1, color="#22A884FF", fill="#22A884FF")+
  annotate("rect", xmin=c(1.012), xmax=c(1.014), ymin=c(0.3043) , ymax=c(0.3045), alpha=1, color="#FDE725FF", fill="#FDE725FF")+
  annotate("rect", xmin=c(1.019), xmax=c(1.021), ymin=c(0.3043) , ymax=c(0.3045), alpha=1, color="#FDE725FF", fill="#FDE725FF")+
stat_function(fun = cor32, color = "#FDE725FF", linetype = "dotted", size = 1)


