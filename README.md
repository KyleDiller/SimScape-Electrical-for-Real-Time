# SimScape-Electrical-Real-Time
Contains realistic power systems models made with SimScape Electrical, ready for parallel execution in real-time.

Since January 2026, the toolbox previously known as SimPowerSystems is not supported anymore by Mathworks.
The global objective of this repository is to provide user with equivalent power systems models in the native SimScape Electrical for Simulink.

SimScape Electrical was not designed by power systems engineers and lacks of some expertise in models, components and methods.
SimPowerSystems in contrast was made with real-time and power systems in mind by Hydro-Quebec.

Nevertheless, SimScape Electrical is a really powerful toolbox, only lacking some twists to make it usable by industries. The SimScape solver engine is a general multi-domain solver and is very similar to a nodal admittance solver like EMTP.

In this repository, we find:

- realistic power systems models, actually inspired by the SimPowerSystems in some case, often ready for real-time deployment in SpeedGoat using concurrent and parallel task processing.
- some optimized model for real-time such as: half-bridge MMC, inverter switching functions, realistic decoupling elements such decouping distributed parameter lines.

I don't have a SpeedGoat simulator so if you want to test the model and report the performance, that would be great. It is also possible to port them to other RTDS with some skills.

Folder description:

CommonFiles: files for commonly used blocks and required by most power system models.

Converter_MultiWinding_9LevelCascadedH-Bridge: a high power converter using shifted multi-winding transformer on the feeder to minimize grid harmonics and 9 level cascaded H-bridge 3-phase inverter. Optimimized for real-time implementation, with decoupled DC stages.

HVDC-MMC-200cellsPerArms: a full HVDC MMC bipolar link with distributed parameter lines with detailed optimized 200 cells per arms MMC model. Ready for concurent/parallel execution on SpeedGoat, using 4 different tasks.

HVDC_bipolar_Thyristor: a full HVDC bipolar link with distributed parameter lines, switched filter/reactive power compensation banks. Ready for concurent/parallel execution on SpeedGoat, using 4 different tasks.

- Cite as Kyle Diller, SimScape Electrical models for real-time, 2026.
