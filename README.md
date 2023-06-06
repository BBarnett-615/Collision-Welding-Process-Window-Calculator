# Collision-Welding-Calculator
Prediction of collision welding process limits in the weld velocity-impact angle plane using analytical methods for similar and dissimilar welds.


Collision welding minimum velocities are calculated using the greater value of (Ultimate Tensile Strength/Density) between the two weld materials.
Jetting limits are calculated using the procedures described by Walsh (1) and implemented numerically similarly to the methods described by de Rossett (2).
Flat interface - wavy interface transition velocities are not calculated.
Upper welding limits can be calculated using a variety of approaches for modeling both stress propagation and interfacial heating/cooling:
  The Acoustic stress model assumes that propagating stresses have a constant velocity equal to the bulk sound velocity based on work by Wittman (3),
  The elastic model assumes stress propagation dependent on the materials' shear modulus based on work by Efremov et al (4),
  The shock model is adapted from 1-D Rankine-Hugoniot approximations developed by the shock physics community and based largely on the text by Meyers (5).

  Thermal models include a moving heat source model based on Wittman (3) and Carslaw & Jaeger (6), 
  an instantaneous heat source model based on Zakharenko (7) and modified to include dissimilar materials based on work by Karkhin (8), 
  and a shock heating model developed by Signetti and Heine (9).
  
  Shock response data is from Marsh (10), other quasistatic material properties are from associated ASM handbooks (11,12).
  
  Works Cited
  
[1] 	Walsh JM, Shreffler RG, Willig FJ. Limiting conditions for jet formation in high velocity collisions. J Appl Phys. 1953;24:349–359.
[2] 	De Rosset WS. Analysis of explosive bonding parameters. Mater Manuf Process. 2006;21:634–638.
[3] 	Wittman RH. The influence of collision parameters of the strength and microstructure of an explosion welded aluminium alloy. Proc 2nd Int Symp Use an Explos Energy Manuf Met Mater. 1973. p. 153–168.
[4] 	Efremov V V, Zakharenko ID, Division S. Determination of the Upper Limit to Explosive Welding. Fiz Goreniya y Vzryva. 1976;3:226–230.
[5] 	Meyers MA. Dynamic Behavior of Materials. John Wiley & Sons, Inc.; 1994.
[6] 	Carslaw HS, Jaeger JC. Conduction of Heat in Solids. 2nd ed. Oxfordshire: Clarendon Press; 1959.
[7] 	Zakharenko ID, Sobolenko TM. Thermal Effects in the Weld Zone in Explosive Welding. Fiz Goreniya y Vzryva. 1971;7:433–436.
[8] 	Karkhin V. Thermal Processes in Welding. Eng. Mater. Springer; 2019.
[9] 	Signetti S, Heine A. Characterization of the transition regime between high-velocity and hypervelocity impact: thermal effects and energy partitioning in metals. Int J Impact Eng [Internet]. 2021;151:103774. Available from: https://doi.org/10.1016/j.ijimpeng.2020.103774.
[10] 	Marsh SP, editor. LASL Shock Hugoniot Data [Internet]. Los Alamos Ser. Dyn. Mater. Prop. University of California Press; 1980.
[11] 	Klueh RL. Properties and selection: irons, steels, and high performance alloys. ASM International; 2005.
[12] 	Committee A. Properties and Selection: Nonferrous Alloys and Special-Purpose Materials. ASM International; 1990.
