#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Kinematic Wave Model (Surface Water)
#Coder: Nate Jones (cnjones7@ua.edu)
#Date: 10/15/2019 
#Purpose: The goal of this script was to begin exploring surface water routing mechanisms, 
#         As of 11/19/2019 -- we decied to focus more on the storage zone. So this script
#         is really not useful. However, I decided to keep it just incase we came back to it 
#         later. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup Workspace------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory (delete this once incoporated into larger workflow)
rm(list=ls())

#download relevant packages
library(tidyverse)

#User defined variables
#Simulation variables
simulation_length<-2       # days
#Channel characteristics
reach_length <- 10000      # Each reach is 10 km in length
reach_slope <- 0.003       # Slope of each reach 
reach_bkf_depth <- 0.1     # Bankful depth
reach_wd_ratio <- 10       # Width to depth ratio of bankful channel
reach_roughness<-0.07      # Mannings N
#Storage Characteristics
storage_flux_ratio<-0.1    #The ratio of 
storage_volume_ratio<-3    #The ratio of V_storage/V_channel
c_diffusion<-0.5           #Turbulent Diffusion Coefficient (m^2/s)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Estimate hydraulic/hydrologic parameters-----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Hydraulic Geometry---------------------------------------------------------
#Estimate bankful width (assume rectangular channel)
reach_bkf_width<-reach_wd_ratio*reach_bkf_depth

#Estimate bankfull flow (cms) using mannings equation
Q_bkf<- 
  (1/reach_roughness)*                             #k/n
  (reach_bkf_depth*reach_bkf_width)*               #Bankful Area
  (reach_bkf_width^(2/3))*                         #Hydraulic Radius ^ (2/3)
  (reach_slope^.5)                                 #slope^0.5

#Estimate Bankful Velocity 
v_bkf<- Q_bkf/(reach_bkf_depth*reach_bkf_width)

#2.2 Inflow hydrograph ---------------------------------------------------------
#Define time step so that the Courant number (v*dt/dx) is equal to 1
timestep<- floor((reach_length/v_bkf))

#Create tibble
df<-tibble(
      #Create 2 time steps
      time = seq(0,simulation_length*86400, timestep),
      #Create flow 
      Q_in = Q_bkf*sin(pi*time/86400) + Q_bkf/10
    ) %>% 
  mutate(Q_in = if_else(Q_in<=Q_bkf/10, Q_bkf/10, Q_in))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Create function to estimate fluxes-----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function
hydro_fun<-function(
  Q_in,                    #Inflow
  Q_in_1,                  #Qout from previous time step
  dt = timestep,           #time step (seconds)
  dx = reach_length,       #Reach length
  slope =  reach_slope,    #Channel Slope
  width = reach_bkf_width, #Channel Width
  n = reach_roughness      #mannings roughness
){
  
  #Kinematic Wave routing~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Estiate model parameters
  beta <- 3/5  #Common Assumption (See Chow 1988)
  alpha<- (n*(width^(2/3))/(slope^0.5))^beta
  
  #Estimate Reach Outflow
  Q_out<- Q_in + (((Q_in^(1-beta))/(alpha*beta*dx))*(Q_in_1-Q_in))*dt
  
  #Export Qout
  #Q_out<-if_else(Q_out<0, 0, Q_out)
  Q_out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Run routing model===================================-----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run model using pmap functionality from purr package
df<-df %>%
  #Run model
  mutate(Q_out = 
     pmap_dbl(
       #Model Inputs
       list(Q_in =   Q_in, 
            Q_in_1 = lag(Q_in, default = 0)), 
       #Function
       hydro_fun)) %>% 
  #Add time_hr 
  mutate(time_hr = time/3600)

#Plot the results
plot(df$time_hr, df$Q_in, type="l", xlab="Time [hr]", ylab="Flow [cms]", ps=12, cex.lab=14/12, cex.axis=10/12)
points(df$time_hr, df$Q_out, col="red", type="l")

