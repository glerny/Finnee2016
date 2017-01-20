# Introduction  
Finnee2016 is a Matlab toolbox that was specifically designed to analyse analytical separation techniques hyphenated with high-resolution mass spectrometry (X-HRMS). Finnee2016 is not a software, it is a collection of classes and functions that allow recovering MS scans from mzML files and performed transformation (baseline drift correction, peak picking, centroid algorithms…), As a toolbox Finnee2016 cannot be used without a Matlab working. Finnee2016 has been designed with Matlab2016a, but most of the functions also work with older versions. While a basic knowledge of Matlab is needed, Finnee2016 has been designed as to allow to use it without programming skills. Operations are run by entering the correct line in Matlab command panel.

# Finnee2016, what is the point?
Every commercial instrument comes with software that allows analysing acquired X-HRMS data, so what the point of developing or learning a new software? Well, commercial software is often design for mass spectrometer used as standalone instruments and does not fully take advantage of the separation side. Chemometrics and bioinformatics apply to X-HRMS information are very dynamic fields especially in proteomics and metabolomics. Excellent software such as [xcms] (https://xcmsonline.scripps.edu) and [mzMine2] (http://mzmine.github.io/) have already been developed. The aim of FInnee2016 is not to copy them using another language.  Finnee2016 is specifically designed for high-resolution mass spectrometry with MS^1 experiment but is not specific for a particular field. Moreover, Finnee2016 ambitions to provide a comprehensive analysis from a single dataset. The key approach is to extract relevant and accurate data that perfectly described the results of the analysis. 

# Development
Finnee2016 is not a finished product. It is in constant development and files are likely to be modified. Finnee2016 also aims to be a collaborative work; you can correct, comment or contact me if you would like some special features. You can also contact me or follow this project via [Research Gate] (https://www.researchgate.net/project/Finnee-A-Matlab-toolbox-for-post-processing-of-separation-techniques-hyphenated-high-resolution-mass-spectrometry).
  
# Declarations
## Copyright
Every file in the Finnee2016 repository, including this wiki, are copyright

BSD 3-Clause License

Copyright (c) 2016, Guillaume Erny
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Acknowledgement

**This work was the result of the projects:**      
+ POCI-01-0145-FEDER-006939 (Laboratory for Process Engineering, Environment, Biotechnology and Energy – UID/EQU/00511/2013) funded by the European Regional Development Fund (ERDF), through COMPETE2020 - Programa Operacional Competitividade e Internacionalização (POCI) and by national funds, through FCT - Fundação para a Ciência e a Tecnologia. 
  
+ NORTE‐01‐0145‐FEDER‐000005 – LEPABE-2-ECO-INNOVATION, supported by North Portugal Regional Operational Programme (NORTE 2020), under the Portugal 2020 Partnership Agreement, through the European Regional Development Fund (ERDF).   

+ IF/00528/2013, supported by FEDER funds through the Operational Programme for Human Potential and by National Funds through FCT - Foundation for Science and Technology.    

**This work would not have been possible without the support of:**    
+ The Laboratory for process Engineering, Environment, Biotechnology and Energy ([LEPABE] (https://paginas.fe.up.pt/~lepabe/) ), my host institution.   

+ [The Laboratory of foodomics] (http://www.cial.uam-csic.es/pagperso/foodomics/) from Madrid, especially Professor Alejandro Cifuentes, Dr Carolina Simó and Tanize dos Santos Acunha.   




***
**Up&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:**  
**Next&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:** [Finnee: The User's Manual] (https://github.com/glerny/Finnee2016/wiki/Finnee:-The-User's-Manual)  
**Previous&nbsp;&nbsp;:** 
***
**Related to:**
