# Marine-Analysis-Verification
Software to do quick analyses/verifications of experiments

How to use the tools:
- Change input parameters
      Study period start and end: yyyy, mm, dd, hh
          NOTE: If you wish to use the tools for a single cycle, the start and end dates should be the same.
      Paths to experiments
      Variables to display: sst, seaice, ssh
      Radiance instruments to use: These are listed in the RUN_TOOLS scripts. Use "all" to compute stats that include all instruments.
      SLURM job parameters
- To run on command line:
      ./RUN_TOOLS...bash

Tools include:
- RUN_TOOLS.innovation.hist.bash
      Generates pdfs of OMB, OMA, obs error on the same figure.
      Run for 1 experiment.
      Gaussian and heavy-tailed Levy distributions are also displayed.
- 
  
