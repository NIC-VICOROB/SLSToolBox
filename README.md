# SLS (SALEM Lesion Segmentation)

SALEM-LS (SLS) is a new algorithm designed to segment WM lesions on MRI FLAIR images. This tool requires either the original T1 weighted image or its tissue segmentation. If the T1 weighted images is found as an input, we firstly co-register T1w image to the FLAIR image. Henceforth all the analysis is performed in this space. We segment the tissues on T1w image, and pre-process the FLAIR image (noise reduction, intensity correction and normalization). Afterwards, the WM lesions are segmented as outliers of the GM distribution. Finally the outliers detected are the potential lesion candidates wich will be filtered by lesion size and contextual information (tissue type and neighbourhood tissue type). Another restriction can be applied consisting of avoiding elongated regions attached to the ventricles. The outlier detection is performed in two steps, firstly we obtain the bigger lesions allowing the more permissive filter. In the second step we apply these filters more conservatively while looking for small lesions. All those parameters can be tunned in order to get a better accuracy, although a set of default values is given.


SLS is currently available to download both as a ready-to-use Matlab function and as a SPM12 library. The useful resources section contains all available software source-codes as zip files.

Please check the Requirements tab before installing it.

Please browse the SALEM webpage for more information (http://eia.udg.edu/salem/slsToolbox).


## Installation as a SPM library

Installing SLStoolbox as a SPM library is also straightforward. Download the toolbox, and extract it into a folder. These are the required steps to run SLS inside the SPM bundle.
 
 *  On SPM, all external libraries are installed into the toolbox folder. This folder is located inside the main SPM directory. Just find where SPM was installed and open the toolbox folder.
 *  Download the toolbox, and extract it into the spm/toolbox folder. Open MATLAB as usual, and type spm in the command-prompt. This will open the main SPM menu. Hence, on the Toolbox list, select the SLStoolbox.
 *  The SPM batch editor will be opened containing the SLStoolbox program.

The program requires two input images:

* A T1-w image to segment the tissues.
* A FLAIR image where the WM lesions will be segemented.

* The program is able to do bias correction and skull-stripping, via SPM internal procedures. We apply them for the T1w image

The resulted lesion segmentation image will be saved in the same directory of the input FLAIR image (FLAIR_PATH/segmLesions.nii). Optionally, intermediate images can also be saved (brain mask, intial candidates, lesion segmentation per iteration and independent tissue masks).
