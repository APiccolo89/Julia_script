"""
LaMEM generation setup 
Main goal: creating a simple julia script that produces .dat file for LaMEM, producing additional plot 
and table related to the initial setup
"""

using Parameters


"""
Abstract type Material Parameters
-> Material
a) Density   -> constant -> thermal -> pressure -> thermal-pressure ->compressible ->Phase diagram 
b) rheology   ->Elastic  -> Viscous                           ->Plastic
                                ->linear                            -> Drucker-Prager
                                    -> Diffusion creep              -> Constant tau yield 
                                    -> Newtonian s.s. 
                                ->non linear
                                    -> Power-law 
                                    -> Dislocation creep
                                    -> Peirl Creep 
c) Thermal -> Basic_Thermal -> Temperature-Dependent Conductivity


"""

abstract type Material end 
# 2nd Level 
abstract type Density <: Material end 
abstract type Rheology <:Material end 
abstract type Thermal <: Material end 
abstract type Diking <: Material end 
# 3rd level
# Density  
abstract type Constant_den <: Density end 
abstract type Thermal_den <: Density end
abstract type Thermal_Pres_den <: Density end
abstract type Phase_Diagram_den <: Density end 
abstract type Compressibility_den <: Density end 
# Rheology 
abstract type Elastic <: Rheology end
abstract type Viscous <: Rheology end
abstract type Plastic <: Rheology end
# Thermal 
abstract type Basic_Thermal <: Thermal end
abstract type T_P_Conductivity <: Thermal end
# 4th level 
abstract type Linear <: Viscous end
abstract type Non_linear <: Viscous end 
abstract type Diffusion_Creep <: Viscous end
abstract type Dislocation_Creep <: Viscous end 
abstract type Peirls_Creep <: Viscous end 












