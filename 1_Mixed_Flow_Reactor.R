#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Mixed Flow Reactor Model
#Coder: Nate Jones (cnjones7@ua.edu)
#Date: 11/20/2019
#Purpose: Explore biogeochemical processing in mixed flow reactor model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup Workspace------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory (delete this once incoporated into larger workflow)
rm(list=ls())

#download relevant packages
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Create MFR Function -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Mixed flow reactor model (Steady State Solution)
#   Equations taken from Chapra's "Surface Water Quality Modeling" text. 
#   See Lecture 9, Section 9.1.3 on mixed flow reactors.  This formulation 
#   relies on Danckwerts boundary condionts. This assumes, 
#     (i)   steady state conditions (i.e., dc/dt = 0)
#     (ii)  exponentional decay along the length of the reactor, i.e. c = exp(gx)
#     (iii) there is no diffusion at the outout (i.e. dc/dt*L=0)

mfr_fun<-function(k_rxn, k_dif, u_storage, c_in, x){

  #Estimate gamma values
  eta<-k_rxn*k_dif/(u_storage^2)
  g1<-u_storage/(2*k_dif)*(1+sqrt(1+4*eta))
  g2<-u_storage/(2*k_dif)*(1-sqrt(1+4*eta))
  
  #Estimate beta values
  b1<-(u_storage*c_in*g2*exp(g2*x))/
      (((u_storage-(k_dif*g1))*g2*exp(g2*x)) -
       ((u_storage-(k_dif*g2))*g1*exp(g1*x)))
  b2<-(u_storage*c_in*g1*exp(g1*x))/
      (((u_storage-(k_dif*g2))*g1*exp(g1*x)) -
       ((u_storage-(k_dif*g1))*g2*exp(g2*x)))
  
  #Estimate Concentration at output of MFR!
  c_out = b1*exp(g1*x) + b2*exp(g2*x)
  
  #Export % removal
  (c_in-c_out)/c_in*100
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Apply the function---------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create inputs for model run
n<-100
inputs<-tibble(
  k_rxn =     exp(log(10)*seq(log10(1e-4),log10(1),length.out=n)),    #From LINXII study -- assuming 1m2 channel xs for now, units in mg/m2/s
  k_dif =     exp(log(10)*seq(log10(1),log10(1000),length.out=n)),    #Just guessing here...
  u_storage = exp(log(10)*seq(log10(0.00001),log10(1),length.out=n)), #1 m/d to 1 m/s, units are m/s
  c_in      = exp(log(10)*seq(log10(0.001),log10(10),length.out=n))
) 
inputs<-expand.grid(inputs) %>% sample_n(10000)

#Define Reach Length
reach_length <- 100 #meters

#Run model using pmap functionality from purr package
df<-inputs %>%
  #Run model
  mutate(p_removal = 
     pmap_dbl(
       #Model Inputs
       list(k_rxn     = k_rxn, 
            k_dif     = k_dif,
            u_storage = u_storage, 
            c_in      = c_in, 
            x         = reach_length), 
       #Function
       mfr_fun)) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Plot-----------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Plot N Removal Vs model inputs --------------------------------------------
par(mfrow=c(4,1))
par(mar=c(4,4,1,1))
par(ps=12)
par(cex.axis=10/12)
par(cex.lab=14/12)
par(mgp = c(2.2,1,0))
par(col="#4C4C4C80")
par(pch=19)
plot(df$k_rxn, df$p_removal, log="x", ylab="% N Removal", xlab = "Reaction Rate [mg/s]")
plot(df$k_dif, df$p_removal, log="x", ylab="% N Removal", xlab = "Diffusion Coefficent [m2/s]")
plot(df$u_storage, df$p_removal, log="x", ylab="% N Removal",xlab = "Storage Velocity [m/s]")
plot(df$c_in, df$p_removal, log="x", ylab="% N Removal", xlab = "Concentration [mg/L]")

#4.2 Plot N removal vs Damkohler and Peclet numbers-----------------------------
#Estimate removal based on Pe and Da
df %>% 
  as_tibble() %>% 
  mutate(
    Da = k_rxn*reach_length/u_storage,
    Pe  = reach_length*u_storage/k_dif) %>% 
  select(p_removal, Da, Pe) %>% 
  filter(Da<25) %>% 
  filter(Pe<10) %>%
  pivot_longer(-p_removal) %>% 
  ggplot() +
    geom_point(aes(y=p_removal, x = value)) + 
    facet_grid(~name, scales = "free") +
    scale_x_log10() +
    theme_bw() +
      ylab("Solute Removal [%]") +
      xlab(NULL)

#4.3 Plot Removal in Pa-Da Space------------------------------------------------
df %>% 
  as_tibble() %>% 
  mutate(
    Da = k_rxn*reach_length/u_storage,
    Pe  = reach_length*u_storage/k_dif) %>% 
  select(p_removal, Da, Pe) %>% 
  filter(Da<25) %>% 
  filter(Pe<10) %>%
  ggplot(aes(x=Da, y=Pe)) +
    geom_point(aes(color = p_removal)) +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_gradient(low="orange", high="dark red", name="% Removal") +
    theme_bw() 
      
    
  
