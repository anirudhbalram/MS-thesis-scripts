1) It takes a SIMBA snapshot and its caesar catalog, and write the relevant information from all its particles in a pandas dataframe (which will be saved as "Basic_Characteristics*")
2) It creates a file with lists of cloud IDs (which will be saved as "Clouds_per_Core*"). Each line represents the IDs that will be passed to each core.
3) It runs despotic on all the particles (or on a random set of those depending on the given parameters), and outputs a table with their CO, C+ and CII luminosities.
