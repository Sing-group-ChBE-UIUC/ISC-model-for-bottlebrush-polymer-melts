1. Target g(r) from the SCMF simulation** 
`*RDF_final_pbc_discrete.c*`
- Input: `“output.txt”, “traj_.xyz”tmax, nbin, ISC_segment(argv[1]), bb_perISC(argv[2])`
- Output: `“RDF0_discrete.txt”`
`*RDF_final_pbc_forblend.c*`
- Input: `”output.txt”, “traj_blend.xyz”, tmax, nbin, Nsc, Nbb, f_branch, ISC_segment(argv[1]), bb_perISC(argv[2])`
- Output: `“RDFAA.txt”, “RDFAB.txt”`
2. Smooth g(r)
`*smooth.py`
- Input: `“RDF__.txt”-ro, gro (argv[1]), spline_factor(s, argv[2])`
- Output: `“RDFsm.txt”, plt_show`
3. g(r) to pmf (Excel or code)
- Output: `“pot0.table”, "potab0.table" - format(index, distance, PMF, MF)`
- Correct extrapolating the short range interaction
4. IBI iteration 
`*IBI_BD.c / Interactions_IBI.h`
- Input: `"Input.txt”, "final0.cfg”, "pot0.table”, "RDF0.txt”`
- Output: `“BD_IBI_.txt”, “thermo_IBI_.txt”, “gri_IBI_.txt”, "pmfi_IBI_.txt”`
- NISC = NISC+NISC-1
`*IBI_BBCP_BD.c / Interaction_IBIBBCP.h`
- Input: `"Input.txt”, "final0.cfg”, "pot0.table”, "potab0.table", "RDFAA.txt”, "RDFAB.txt"`
6. Visualization of convergence 
`gr_vis.py gri_IBI_.txt RDF0.txt`
