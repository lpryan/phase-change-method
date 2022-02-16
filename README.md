# phase-change-method
A method that returns the following based of the input that is given:

      ---- INPUTS ---- 

Required inputs:
  - Compound
  - Temp_Init
  - Temp_Final

Requires ONE of the following:
  - mass
  - moles
  - energy


      ---- TYPES OF INPUTS ----

Compound [custom class inputed by the user that has 9 possible elements]: 
  - molar [molar mass]
  
  - vapor [ of fusion]
  - fusion [ of fusion]
  
  - melt [the melting point of the compound]
  - boil [the boiling point of the compound]
  
  - Cp_g [the heat capacity of the compound as a gas]
  - Cp_1 [the heat capacity of the compound as a liquid]
  - Cp_s [the heat capacity of the compound as a solid]
  - Cp_a [the heat capactiy of the compound when the state of matter is unknown but constant]


Temp_Init [initial temperature of the compound]
Temp_Final [final temperature of the compound]

mass [mass of the compound]
moles [number of moles in the compound]

energy [the amount of energy used in the phase change]






