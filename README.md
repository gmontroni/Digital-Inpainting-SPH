# Digital Inpainting via Block Replication and the SPH method

Digital inpainting techniques are intended to recompose missing parts of an image or remove undesirable objects from it. Considering this study field, one of the main issues faced by inpainting methods is their high computational costs, as well as the visual aspect of the outputs under many specific occasions, leading to poor quality reconstructions. Aiming at exploiting and adapting the classic Smoothed-Particle Hydrodynamics (SPH) method in the context of digital inpainting, this research focuses on the study and development of a new inpainting technique which relies on the traditional SPH formulation. Our goals was to apply the usual SPH framework in the task of digital inpainting by using the patch-based propagation paradigm, i.e., the so-called patch-based image inpainting, where each SPH particle is interpreted as a full patch of pixels in our approach.

One of the first restoration processes was presented by Bertalmio et al, where the original formulation in order to recover an image (a) assumes as input the (b) deteriorated image $I:Ω → ℝ, Ω ⊂ ℝ²$, (c) a guide mask 𝐃 ⊂ Ω, that discriminates the information from the deteriorated region of the image $I$.

<p align="center">
	<img align="center" width="700" height="250" src="https://user-images.githubusercontent.com/96217617/182721437-ab4f87d6-4229-4717-9598-ff0852a9cc49.png">
</p>

Below we present an application of the proposed method (c) compared to the Criminisi (a) and G-SPIR (b) algorithms. 

<p align="center">
	<img align="center" width="700" height="250" src="https://user-images.githubusercontent.com/96217617/182721478-d5de1419-5122-4249-92e7-57e06fbd5412.png">
</p>

More details about these procedures will be presented in the DOI: https://doi.org/10.5540/03.2021.008.01.0424.
