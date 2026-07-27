import itk
import numpy as np
import opengate_core as g4
from .generic import GenericSource
from ..image import (
    read_image_info,
    get_info_from_image,
    update_image_py_to_cpp,
    compute_image_3D_CDF,
)
from ..utility import ensure_filename_is_str
from ..base import process_cls


class VoxelSource(GenericSource, g4.GateVoxelSource):
    """
    VoxelSource = 3D activity distribution.
    Sampled with cumulative distribution functions.
    """

    # hints for IDE
    image: str

    user_info_defaults = {
        "image": (
            None,
            {
                "doc": "Filename of the image of the 3D activity distribution "
                "(will be automatically normalized to sum=1)",
                "is_input_file": True,
            },
        )
    }

    def __init__(self, *args, **kwargs):
        self.__initcpp__()
        super().__init__(self, *args, **kwargs)
        # the loaded image
        self.itk_image = None

    def __initcpp__(self):
        g4.GateVoxelSource.__init__(self)   

    def set_transform_from_user_info(self):
        # get source image information
        src_info = read_image_info(str(self.image))
        # get the pointer to SPSVoxelPosDistribution
        pg = self.GetSPSVoxelPosDistribution()
        # update cpp image info (no need to allocate)
        update_image_py_to_cpp(self.itk_image, pg.cpp_edep_image, False)
        # set spacing
        pg.cpp_edep_image.set_spacing(src_info.spacing)
        # set origin (half size + translation and half-pixel shift)
        c = (                               
            -src_info.size / 2.0 * src_info.spacing
            + self.position.translation
            + src_info.spacing / 2.0
        )
        pg.cpp_edep_image.set_origin(c)

    def cumulative_distribution_functions(self):
        """
        Compute the Cumulative Distribution Function of the image
        Composed of: CDF_Z = 1D, CDF_Y = 2D, CDF_X = 3D
        """
        cdf_x, cdf_y, cdf_z = compute_image_3D_CDF(self.itk_image)

        # set CDF to the position generator
        pg = self.GetSPSVoxelPosDistribution()
        pg.SetCumulativeDistributionFunction(cdf_z, cdf_y, cdf_x)

    def initialize(self, run_timing_intervals):
        # read source image
        self.itk_image = itk.imread(ensure_filename_is_str(self.image))

        # compute position
        self.set_transform_from_user_info()

        # create Cumulative Distribution Function
        self.cumulative_distribution_functions()

        # FIXME -> check other option in position not used here

        # initialise standard options (particle energy, etc.)
        # we temporarily set the position attribute to reuse
        # the GenericSource verification
        GenericSource.initialize(self, run_timing_intervals)





class VoxelizedPromptGammaTLESource(GenericSource, g4.GateVoxelizedPromptGammaTLESource):
    """
    VoxelizedPromptGammaTLESource = distribution 4D des spectres de prompt-gamma (PG).
    - Echantillonnage spatial : CDF 3D (somme sur les bins d'energie), identique a VoxelSource
    - Echantillonnage en energie : CDF par voxel, construite en C++ 
    """
 
    # indications pour l'IDE
    image: str
    Emin: float
    Emax: float
    nbBins: int
 
    user_info_defaults = {
       # --- image energie ---
       "image": (
           None,
           {
               "doc": "Nom du fichier de l'image 4D (spectres d'énergie PG par voxel) : [bins, Z, Y, X]",
               "is_input_file": True,
               "dynamic": True,
           },
        ),
       "Emin": (
           0,
           {
               "doc": "energie minimale du spectre PG (MeV)",
           },
        ),
       "Emax": (
           10,
           {
               "doc": "energie maximale du spectre PG (MeV)",
           },
        ),
       "nbBins": (
           250,
           {
               "doc": "Nombre de bins en énergie du spectre PG",
           },
        ),
        "n_protons": (
        1000,
        {
            "doc": "Nombre de particules primaires equivalentes en amont "
            "(ex : protons). Le nombre de gammas a emettre (source.n) "
            "est calcule automatiquement a partir de la somme totale de "
            "l'image PG (yield par particule primaire), multiplie par "
            "cette valeur. Ecrase toute valeur de 'n' deja fixee.",
        },
        ),

    }
 
    def __init__(self, *args, **kwargs):
        self.__initcpp__()
        super().__init__(self, *args, **kwargs)
        # l'image 4D d'energie charge (image itk, [bins, Z, Y, X])
        self._current_itk_image = None
        # l'image 3D obtenue en sommant sur les bins d'energie
        # (utilise pour la CDF spatiale + la geomerie, meme role que VoxelSource._current_itk_image)
        self._current_itk_image_3d = None
 
    def __initcpp__(self):
        g4.GateVoxelizedPromptGammaTLESource.__init__(self)
 
    def set_transform_from_user_info(self):
        # identique a VoxelSource mais utilise _current_itk_image_3d (image 3D, pas 4D)
        src_info = get_info_from_image(self._current_itk_image_3d)
        # recupere le pointeur vers SPSVoxelPosDistribution
        pg = self.GetSPSVoxelPosDistribution()
        # met a jour les informations de l'image cote C++ (pas besoin d'allouer)
        update_image_py_to_cpp(self._current_itk_image_3d, pg.cpp_edep_image, False)
        # definit le spacing
        pg.cpp_edep_image.set_spacing(src_info.spacing)
        # definit l'origine (demi-taille + translation + decalage d'un demi-pixel)
        c = (                               # les coordonnees de lorigine de l'image
           -src_info.size / 2.0 * src_info.spacing
           + self.position.translation
           + src_info.spacing / 2.0
        )
        pg.cpp_edep_image.set_origin(c)
 
    def cumulative_distribution_functions(self):
        """
        Calcule la CDF spatiale de l'image de yield 3D
        (obtenue en sommant les spectres d'energie PG 4D sur les bins d'energie).
        Composee de : CDF_Z = 1D, CDF_Y = 2D, CDF_X = 3D
        Identique a VoxelSource mais utilise _current_itk_image_3d
        """
        cdf_x, cdf_y, cdf_z = compute_image_3D_CDF(self._current_itk_image_3d)
        # recupere le pointeur vers SPSVoxelPosDistribution
        pg = self.GetSPSVoxelPosDistribution()
        pg.SetCumulativeDistributionFunction(cdf_z, cdf_y, cdf_x)
 
    def update_pg_images(self, filename_energy):
        # lecture de l'image 4D d'energie PG : [bins, Z, Y, X]
        self._current_itk_image = itk.imread(ensure_filename_is_str(filename_energy))
        arr_4d = itk.array_from_image(self._current_itk_image)


        # calcule automatiquement le nbr de PG a emettre a partir 
        # du yield total contenu dans l'image (somme sur tous les voxels 
        # et tous les bins d'energie), multiplie par le nombre de 
        # particules primaires equivalentes demande
        yield_per_proton = float(arr_4d.sum()) 
        self.n = int(round(self.n_protons * yield_per_proton)) 

 
        # somme sur les bins d'energie -> image de yield 3D (meme role que l'image de VoxelSource)
        arr_3d = np.sum(arr_4d, axis=0)
        self._current_itk_image_3d = itk.image_from_array(arr_3d.astype(np.float32))
        spacing_4d = self._current_itk_image.GetSpacing()
        origin_4d = self._current_itk_image.GetOrigin()
        self._current_itk_image_3d.SetSpacing([spacing_4d[0], spacing_4d[1], spacing_4d[2]])
        self._current_itk_image_3d.SetOrigin([origin_4d[0], origin_4d[1], origin_4d[2]])

 
        # calcul de la geometrie (spacing, origine) - identique a VoxelSource
        self.set_transform_from_user_info()
 
        # calcul de la CDF spatiale - identique a VoxelSource
        self.cumulative_distribution_functions()
 
        # initialisation de la CDF en energie par voxel, cote C++ (option 3 : triple boucle en C++)
        self.InitializeEnergyCDF(arr_4d, self.Emin, self.Emax, self.nbBins)
 
    def initialize(self, run_timing_intervals): 
        # run_timing_intervals : parametre impose par l'interface GeneriSource
        # Non utilise directement ici car pas de decroissance temporelle pour cette sourc
        self.update_pg_images(self.image)
        # note : l'energie n'est pas definie via source.energy ici,
        GenericSource.initialize(self, run_timing_intervals)
 
 
process_cls(VoxelSource)
process_cls(VoxelizedPromptGammaTLESource)
 

