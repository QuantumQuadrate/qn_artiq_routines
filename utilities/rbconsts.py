from math import pi
from utilities.physconsts import *

## Rb87 constants
mRb = 1.4192261e-25 # [kg]
I = 3/2 # nuclear spin
s = 0.5 # electron spin
nu_hf = 6.83468261090429 # [GHz]
gamma_D2 =2*pi*6.0659e6 ; # [rad/s]
lambda_D2 = 7.8024120968613e-7 # [m]
lambda_D1 = 7.9497885098e-7 # [m]
omega_D2 = 2*pi*c/lambda_D2
omega_D1 = 2*pi*c/lambda_D1
nu_D2 = omega_D2/(2*pi)
gS = 2.00023
gL = 1 # for ground state?
gI = -0.000995
D2_Isat = (5/7)*(hbar*gamma_D2*omega_D2**3
             /(12*pi*c**2)) 
    # saturation intensity for D2 cooling W/m^2

## reduced matrix elements from Steck
D2_MatElem = 3.584e-29 # <J=1/2||er||J'=3/2> [C*m]
D1_MatElem = 2.537e-29 # <J=1/2||er||J'=1/2> [C*m]

# the keys are nested in order n, L, J, F, nu
# excited state frequencies are given as the transition frequency
# to the ground state c.o.m. plus the hyperfine shift for the excited level.
hf_states = {5:
                    {0: 
                        {0.5:
                            {1: -4.271676631815196,
                             2: 2.563005979089114}
                        },
                     1:
                        {0.5:
                            {1: 377106.953053064,
                             2: 377107.769709364},
                         1.5:
                            {0: 384230.1823946245,
                             1: 384230.2546166565,
                             2: 384230.4115571805,
                             3: 384230.6782093585}
                        }
                    }
                }

# hf_states[n][l][j][f]
                
                