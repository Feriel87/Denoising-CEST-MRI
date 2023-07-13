### Non Local mean Coherence Enhancing Diffusion filter (NLmCED) ###
NLmCED filter is a new method based on combination between two filters such as the Non-Local mean filter and the Anisotropic Diffusion tensor method with an estimator of noise Rician in MR images.
The NLM filter is a technique that utilizes similarity between image patches to reduce noise. It replaces each pixel with a weighted average of similar patches in the image, effectively preserving edges and fine details while reducing noise.

On the other hand, anisotropic diffusion is a diffusion-based filter that selectively smooths an image while preserving strong edges. It achieves this by diffusing the image in a manner that is guided by the local gradient information.

By combining these two filters, NLmCED can benefit from the noise reduction capabilities of the NLM filter while maintaining the edge-preserving characteristics of anisotropic diffusion. This hybrid approach often produces superior results by effectively reducing noise while preserving important image structures and details.

NLmCED filter is an iterative filter with 4 parameters to be set:
-	Wind: A single value to fix the kernel size of the Gaussian window to create a tensor structure and a tensor diffusion matrix. 
-	ρ: A standard deviation of the Gaussian kernel for the creation of the structure tensor.
-	α: A single value to control the diffusion tensor matrix.
-	Iter: Iteration number.

  ### Citation ### 
  [1] F. Romdhane, F. Benzarti, A. Hamid, A new method for three-dimensional magnetic resonance images denoising. International Journal of Computational Vision and Robotics 2018 8:1, 1-17. DOI:10.1504/IJCVR.2018.090012
