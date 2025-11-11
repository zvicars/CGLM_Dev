#Example Project Description
#Project Name: Flat Slab Contact Angle 

A more involved example illustrating the procedure for adding surfaces to the simulation
box. Simulations of this type can be used to determine the effective contact angle of a 
surface with a given attractive strength (epsilon, as a component of the Lennard-Jones potential).
However, as that is somewhat involved, a simpler output of the average shape of the liquid interface
is shown instead.

To run this example, just execute the "flat_slab.sh" script.

An STL is produced from the "phi" file (an external potential describing the interaction of liquid
water and a surface) called solid.stl. The average shape of the liquid meniscus from the sim can be seen in 
liquid.stl, which is indicative of the macroscopic water contact angle of the surface (if liquid meniscus is convex, surface
is hydrophobic). The best way to visualize it is by importing both .stl files into the same viewport; my 
preferred program for this is Blender. As the liquid and solid surfaces are periodic in x, the meshes are 
hollow in that dimension. Looking along the y,z plane, the solid surface is manifest as a long, rectangular prism 
and the liquid meniscus is the rough surface. The output for the selected value of epsilon (0.2) for the attractive potential
results in a slightly hydrophobic slab.

Admittedly this is unclear without context, but it's a simplified version of the calculation used in Figure S2 of the following publication:
https://www.science.org/doi/10.1126/sciadv.adu8349

With the SI being available at this link:
https://www.science.org/doi/suppl/10.1126/sciadv.adu8349/suppl_file/sciadv.adu8349_sm.pdf

Admittedly, deferring to my paper feels lazy, but it's the cleanest way I can think of to provide some additional clarity as to what I'm showing.
If more information is desired, I would recommend looking at the "Simulations" section in the main text (towards the end), Section S1 in the SI, 
and Section S2 in the SI. On my end, I need to improve the documentation and examples to make them more self-contained for posterity.