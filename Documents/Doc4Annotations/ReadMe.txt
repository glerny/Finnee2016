Guillaume Erny
Aveiro, 02/10/2020

This folder contains all the information related to the annotation of untargeted metabolomics clusters by Finnee.

The folder <Documents> contains the databases, either the original files or the .mat files
The folder <Matlab scripts> contains the functions and live sripts used to created the .mat database.

*** Description of files
**  In folder <Documents>
*   hmdb_metabolites.xml
- downloaded or created: 15/06/2020
- size: 3.87 GB
- owner: The Human Metabolome Database (HMDB) (https://hmdb.ca/), 
- copyright: HMDB is offered to the public as a freely available resource. Use and re-distribution of the data, in whole or in part, for commercial purposes requires explicit permission of the authors and explicit acknowledgment of the source material (HMDB) and the original publication (see below). We ask that users who download significant portions of the database cite the HMDB paper in any resulting publications.
Please cite:
Wishart DS, Tzur D, Knox C, et al., HMDB: the Human Metabolome Database. Nucleic Acids Res. 2007 Jan;35(Database issue):D521-6. 17202168 
Wishart DS, Knox C, Guo AC, et al., HMDB: a knowledgebase for the human metabolome. Nucleic Acids Res. 2009 37(Database issue):D603-610. 18953024 
Wishart DS, Jewison T, Guo AC, Wilson M, Knox C, et al., HMDB 3.0 — The Human Metabolome Database in 2013. Nucleic Acids Res. 2013. Jan 1;41(D1):D801-7. 23161693 
Wishart DS, Feunang YD, Marcu A, Guo AC, Liang K, et al., HMDB 4.0 — The Human Metabolome Database for 2018. Nucleic Acids Res. 2018. Jan 4;46(D1):D608-17. 29140435 
- exerg:
  <version>4.0</version>
  <creation_date>2005-11-16 15:48:42 UTC</creation_date>
  <update_date>2019-01-11 19:13:56 UTC</update_date>
  <accession>HMDB0000001</accession>

* hmdb_metabolites.mat
- downloaded or created: 23/06/2020
- size: 2.17 MB
- owner: Guillaume Erny (guillaume@fe.up.pt)
- copyright: None
- Descriptions:
Summary of the key information for the metabolites contain in hmdb_metabolites.xml. If loade in Matlab, the HMDB table will be loaded, with 4 columns and >100,000 lines. 
	{'accesion'}    acess code to hmdb HMDB0000XXX
	{'name'}    	Name of the molecules 
	{'formula'}	Formula
	{'monisotopic'}	Monoisotopic mass as recorded in HMDB

* AdductsRules.mat
- downloaded or created: 02/10/2020
- size: 2.17 MB
- owner: Guillaume Erny (guillaume@fe.up.pt)
- copyright: 
- Descriptions:
AdductsRules is a 5xn table that summarise the rules for generating adducts. The table was created using information from the Fiehn Lab (https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/) and from The University of Washington's Proteomics Resource (UWPR: http://www.proteomicsresource.washington.edu/protocols05/esi_background_ions_repeat_units.php).
The columns are:
	{'Struct'}    		General structure of the adducts, e.g. nM+X where M is the motif, n the number of the repetition of the motif, X the elements lost or gained,
	{'Delta_mass'}  	mass difference due to X (positif or negatif)
  	{'Charge'}  		Charge of the adducts
	{'r'}    		structure array representative of the chemical formulas
	{'IsBaseAdduct'}	1 if the adduct structure could be the most important adducts for a given intrant molecule, zero otherwise.

* IsopicEnveloppes4ESIpos.mat
- downloaded or created: 08/10/2020
- size: 362 MB
- owner: Guillaume Erny (guillaume@fe.up.pt)
- copyright: 
- Descriptions:
IsoSpectra is a 5xn table that contain the isotopic envelop spectra of plausible ions from positive ESI MS
	{'mf'}    	molecular formula of the target ions in the Hill notation
	{'Z'}		mass difference due to X (positif or negatif)
  	{'mi'}  	Monoisotopic mass
	{'BPmz'}    	Mass of the mpss intensed isotopes in the envelop
	{'MS'}		Isotopic enbelop spectre, obtained using ChemCalc.

