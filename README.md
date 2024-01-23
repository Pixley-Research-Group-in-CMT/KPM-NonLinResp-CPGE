# KPM-NonLinResp-CPGE

Physics of the repo:

Model based computation of the circularly polarized photogalvanic effect (CPGE) based in 3D chiral Weyl semimetals using a triple Chebyshev expansion within the kernel polynomial method (KPM) approach in https://arxiv.org/a/pixley_j_1.html. Based on a time reversal broken Weyl semimetal, we break the mirror symmetries making it chiral with Weyl nodes at two different places in energy. As a result, the CPGE

<img width="493" alt="image" src="https://github.com/Pixley-Research-Group-in-CMT/KPM-NonLinResp-CPGE/assets/5032322/62a568e2-9d07-4a20-ab05-1c6b0d6a0c67">

becomes quantized for optical frequencies encompassing the Weyl code above the Fermi energy. The dispersion and CPGE in the clean limit showing the quantization is below.

<img width="604" alt="image" src="https://github.com/Pixley-Research-Group-in-CMT/KPM-NonLinResp-CPGE/assets/5032322/56372383-b42b-43e5-b095-87ca5932f5d3">

This repository evaluates this CPGE within a KPM based appraoch on GPUs. It expands the Kubo expression below (inovlving 3 current operatrors) with 3 Chebyshev based expansions. 

<img width="1120" alt="image" src="https://github.com/Pixley-Research-Group-in-CMT/KPM-NonLinResp-CPGE/assets/5032322/b3b4279e-f508-48d8-9766-5206ca02280a">

This gives rise to the KPM expression

<img width="479" alt="image" src="https://github.com/Pixley-Research-Group-in-CMT/KPM-NonLinResp-CPGE/assets/5032322/1d366e81-a6f3-485b-a172-eaba7f2cef6d">

where $\Lambda$ is defined in the appendix of the draft.


How to use the code base in the repo:
