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
n_reaches <- 10        # 10 reaches
reach_length <- 1000   # Each reach is 1000 m in length
reach_slope <- 0.001   # Slope of each reach 
reach_bkf_depth <- 0.5 #Bankful depth
reach_wd_ratio <- 3    #Width to depth ratio of bankful channel
reach_noflow <-0.05    #No flow threshold (as a proportion of bankful flow)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Hydrqulic geometry---------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

