#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Initial Toy Model
#Coder: Nate Jones (cnjones7@ua.edu)
#Date: 10/15/2019 
#Purpose: Begin exploring Toy Model to explore drying and its impact on Biogeochemistry
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Potential Workflow
  #-10 reaches w/ user specified hydraulic geometry 
  #-each reach has a no-flow threshold 
  #-when reach volume > no-flow threshold, water is routed using Muskigum-Cunge method 
  #-when reach volume < no-flow threshold, water "infiltrates" at user specified rate 
  #-I'm trying to implement "flux-tracking" to estimate residence time. We'll see if that happens by Thursday.
  #-For biogeochemistry: looking for input! 
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup Workspace------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory (delete this once incoporated into larger workflow)
rm(list=ls())

#download relevant packages
library(tidyverse)

#User defined variables
reach_number <- 10     # 10 reaches
reach_length <- 100    # Each reach is 100 m in length
reach_slope <- 0.001   # Slope of each reach 
reach_bkf_depth <- 0.5 # Bankful depth
reach_wd_ratio <- 2    # Width to depth ratio of bankful channel
reach_noflow <-0.05    # No flow threshold (as a proportion of bankful flow)
reach_roughness<-0.03  # Mannings N

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Estimate Bankfull Flow (Mannings)------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Estimate bankful width (assume rectangular channel)
reach_bkf_width<-reach_wd_ratio*reach_bkf_depth

#Estimate bankfull flow (cms) using mannings equation
Q_bkf<- 
  (1/reach_roughness)*                             #k/n
  (reach_bkf_depth*reach_bkf_width)*               #Bankful Area
  ((reach_bkf_depth*2+reach_bkf_width)^(2/3))*     #Hydraulic Radius ^ (2/3)
  (reach_slope^.5)                                 #slope^0.5

#Estimate Bankful Velocity 
v_bkf<- Q_bkf/(reach_bkf_depth*reach_bkf_width)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Create hydrograph----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create tibble
df<-tibble(
  #Create 28 day period
  time = seq(0,(24*28)),
  #Create flow 
  Q_upstream = c(0,rep(c(rep(Q_bkf, 24), rep(0, (24*6))), 4))
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Create function to estimate outflow [Kinematic Wave Model]-----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mk_routing<-function(Q_in,   #Inflow
                     Q_n,    #Qout from previous time step
                     v_bkf,  #Bankful velocity 
                     w_bkf,  #Bankful widt
                     slope,  #channel slope
                     dt,     #time step
                     dx,     #chanel length
                     threshold = reach_noflow #no flow threshold
                     ){

  #Estimate dynamic Muskigum-Cunge variables
  c <- (5/3)*v_bkf
  k <- dx/c
  q0<- Q_bkf/w_bkf
  x <- 0.5*(1-(q0/slope/c/dx))
  
  #Estimate Muskigum-Cunge Cofficients
  C_0 <- ((-1*k*x)+(0.5*dt))/((k*(1-x))+(0.5*dt))
  C_1 <- ((k*x)+(0.5*dt))/((k*(1-x))+(0.5*dt))
  C_2 <- ((k*(1-x))-(0.5*dt))/((k*(1-x))+(0.5*dt))
  
  #Estimate Q_out
  Q_out<-C_0*Q_in+C_1*Q_in+C_2*Q_n
    
  #Test to see if Q_out<threshold
  Q_out<-if_else(Q_out>threshold, Q_out, 0)
    
  #Export Q_out 
  Q_out
}


