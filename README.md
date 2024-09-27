# Fourier-Series-Shape-Analysis-Fibroblasts
Shape analysis of two-dimensional fibroblast shapes using Fourier Series and Principal Component Analysis. This code is based on work from: https://www.palass.org/publications/newsletter/palaeomath-101/palaeomath-part-25-centre-cannot-hold-ii-elliptic-fourier-analysis.

To run this example, download all files and keep them all in the same directory or folder. Simply open the .m file inside the folder that has the same name as the folder (e.g., 'segmentation_fourier_series_fit_DOA.m') in Matlab, press Run, and make sure you change the current directory to the directory containg the code and data files. 

This code reproduces shape analysis of fibroblast cells as done in our pre-print here: https://www.biorxiv.org/content/10.1101/2024.05.15.594418v1.abstract. In general, the folder "1_segmentation_fourier_series_fit_DOA" reproduces figures 1 and 2 in the manuscript, "2_PCA_eigenshape_KDE" reproduces figure 3, "3_logistic_regression_model_k_fold_cross_validation" reproduces figures 4 and 5, and "4_grading_test_dataset" reproduces figure 6. Please note that for size efficiency, only one representative example is included in "1_segmentation_fourier_series_fit_DOA". In practice, this code can be modified to iterate through several imaging datasets. In addition, the Fourier series coefficients for all cells in the test data set are compiled into an Excel sheet in "2_PCA_eigenshape_KDE" for size efficiency as well. Lastly, the test data set is also compiled into an Excel sheet in "4_grading_test_dataset" for size efficiency. 

We have included sufficient comments within every .m file and functions so that our results may be reproduced and that the methodology may be extended to other use cases. 

Please note that the magnitudes of PC axes 1 and 3 were flipped for visual effect in the pre-print. This does not change the meaning of the axes.  

We have also compiled a specific list of applicable references/links that demonstrate the implementation of the aforementioned methods:

Coherence:
- https://doi.org/10.1007/978-3-319-28549-8_3
- https://bigwww.epfl.ch/demo/orientationj/
- 
Logistic Regression Models:
- https://www.mathworks.com/help/stats/mnrfit.html
- 
K-fold cross validation:
- https://www.mathworks.com/discovery/cross-validation.html
