# Modelling _Asparagopsis_ growth in Australian coastal waters

This repository contains all the functions and code that generated the data used in the publication (link to go here).

- The main workflow can be found in `_all_projects_RUN.R`. This runs the scripts in `R-scripts` in the correct order. Note that many of these scripts reference data that is not stored on Github, but is publicly available. For information on what data was used and how to access it, please contact me (stormeyseas).
- The scripts used to parameterise both _Asparagopsis_ species can be found in the `parameterisation` folder.
- The main macroalgae growth model functions are not contained in this repository. They can be found as an R package available here: [macrogrow](https://github.com/stormeyseas/macrogrow.git).
- The minimal code used to generate the figures found in the published paper can be found in `docs-working/all_for_manuscript.qmd`, and the rendered files are located in `docs-rendered/all_for_manuscript.html`.
- All other documents in `docs-working` were used in the data generation and analysis of this project.
