import os
from numpy import source
import nibabel as nib
from config_manager import FLAGS

def dcm_to_nii(dcm_path, nii_path):
    # converts DICOM to Nifti. the input should be only the first dcm file from the subject, e.g., /Image1.dcm)
    convert = 'mri_convert -i ' + dcm_path + ' --in_type dicom --out_type nii -o ' + nii_path + '.nii'  
    os.system(convert)

    
def fov_to_RAS(input_scan):
    # set FOV orientation to RAS+ with nibabel
    scan = nib.load(input_scan)
    scan = nib.as_closest_canonical(scan)
    nifti_ras = input_scan.split('.nii')[0] + '_ras.nii'
    nib.save(scan, nifti_ras)    

    
def skull_strip(input_scan):
    template0 = FLAGS.template0
    template_cerebellum = FLAGS.template_cerebellum
    template_cerebellum_mask = FLAGS.template_cerebellum_mask
    
    extract_brain = 'antsBrainExtraction.sh -d 3 -a ' + input_scan + ' -e ' + template0 + ' -m ' + template_cerebellum +' -f '+ template_cerebellum_mask + ' -c 3x1x2x3 -o ' + input_scan.split('.nii')[0]
    os.system(extract_brain)
    
def apply_mask(input_scan):        
        input_dir = '/'.join(input_scan.split('/')[:-1])
        apply_mask = 'fslmaths ' + input_scan + ' -mul ' + input_dir + '/BrainMask.nii.gz'    
        os.system(apply_mask)    
    
def decompress_dicom(dcm_path):
    # Convert the current DICOM format to a readable format for MRtrix3 using other DICOM tools 
    # First install "dcmtk" tool
    # IMPORTANT: make sure to make a copy of the original DICOM files, 
    # since these commands will overwrite the current one with the decompressed one
    run_ID =  sorted(os.listdir(dcm_path))
    for _, img in enumerate(run_ID):

        in_out_file_name = '/'.join([dcm_path, img])
        out_file_tmp = '/'.join([dcm_path, 
                                (img.split('.')[0])+'_tmp.dcm' 
                                ]) # img[0:img.find('.')]
        command0 = ' '.join(['dcmdjpeg',
                            in_out_file_name, 
                            out_file_tmp
                            ]) 
        command1 = ' '.join(['mv', 
                            out_file_tmp, 
                            in_out_file_name,
                            '--force' 
                            ])
        os.system(command0+' && '+command1)

def dcm_to_mif(dcm_path, des_path):
    # convert list of dcm files to one mif file
    command0 = ' '.join(['mrconvert', 
                        (dcm_path), #+'/' 
                        (des_path+dcm_path.split('/')[-1]+'.mif')
                        ])
    os.system(command0)

def denoising(path_to_dwi_mif):    
    # path_to_dwi_mif is the path to directory where dwi.mif data exist
    # Denosing Phase
    cmd_extract_noise = ' '.join(['dwidenoise',
                                  path_to_dwi_mif+'/dwi.mif',
                                  path_to_dwi_mif+'/dwi_den.mif',
                                  '-noise',
                                  path_to_dwi_mif+'/noise.mif' 
                                ])
    cmd_subtract_noise = ' '.join(['mrcalc',
                                    path_to_dwi_mif+'/dwi.mif',
                                    path_to_dwi_mif+'/dwi_den.mif',
                                    '-subtract',
                                    path_to_dwi_mif+'/residual.mif'                                
                                  ]) # is this necessary? 

    os.system(cmd_extract_noise+' && '+cmd_subtract_noise) 

    # Gibbs ringing removal
def gibbs_ringging_removal(path_to_dwi_mif):
    cmd_extract_gibs = ' '.join(['mrdegibbs', 
                                 path_to_dwi_mif+'/dwi_den.mif', 
                                 path_to_dwi_mif+'/dwi_den_unr.mif', 
                                 '-axes 0,1'
                                 ])
    cmd_subtract_gibs = ' '.join(['dwifslpreproc -rpe_none -pe_dir AP',
                                  path_to_dwi_mif+'/dwi_den_unr.mif',
                                  path_to_dwi_mif+'/dwi_den_unr_pre.mif',
                                  '-nocleanup -eddy_options " --slm=linear "'
                                  ])
    # dwipreproc dwi_raw_den_unr.mif dwi_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -eddy_options " --slm=linear"


    os.system(cmd_extract_gibs+' && '+cmd_subtract_gibs)


    ## detect number of outliers
    #Calculate the number of outlier slices
    #cd dwipreproc-tmp-*
    #totalSlices=`mrinfo dwi.mif | grep Dimensions | awk '{print $6 * $8}'`
    #totalOutliers=`awk '{ for(i=1;i<=NF;i++)sum+=$i } END { print sum }' dwi_post_eddy.eddy_outlier_map`
    #echo "If the following number is greater than 10, you may have to discard this subject because of too much motion or corrupted slices"
    #echo "scale=5; ($totalOutliers / $totalSlices * 100)/1" |bc | tee percentageOutliers.txt

# Bias field correction
def bias_field_correction(path_to_dwi_mif):
    cmd_extract_bias = ' '.join(['dwibiascorrect ants', 
                                 path_to_dwi_mif+'/dwi_den_unr_pre.mif',
                                 path_to_dwi_mif+'/dwi_den_unr_pre_unbis.mif',
                                 '-bias bias.mif'
                                ]) 
    os.system(cmd_extract_bias)

# Up scaling (resize voxel)
def upscaling_mif_file(path_to_dwi_mif):
    cmd_subtract_upscaling =' '.join(['mrgrid',
                                     path_to_dwi_mif+'/dwi_den_unr_unbis.mif',
                                     'regrid -vox 1.3',
                                     path_to_dwi_mif+'/dwi_den_unr_pre_unbia_up.mif' 
                                     ])    
    os.system(cmd_subtract_upscaling)

# Registration

def register_dwi_to_t1c(path_to_dwi_mif, path_to_ni_files, path_to_transform_matrix):
    #!dwiextract dwi_den_unr_pre_unbia_up.mif - -bzero | mrmath - mean mean_b0_AP.nii.gz -axis 3
    #!flirt -dof 6 -cost normmi -in mean_b0_AP.nii.gz -ref t1cBrainExtractionBrain.nii  -omat T_fsl.txt
    #!transformconvert T_fsl.txt mean_b0_AP.nii.gz t1cBrainExtractionBrain.nii flirt_import T_DWItot1.txt && rm T_fsl.txt
    #!mrtransform -linear T_DWItot1.txt dwi_den_unr_pre_unbia_up.mif dwi_den_unr_pre_unbia_up_to_t1.mif --force
    cmd_b0_extract = ' '.join(['dwiextract',
                                     path_to_dwi_mif+'/dwi_den_unr_unbia_up.mif - -bzero | mrmath - mean',
                                     path_to_dwi_mif+'/mean_b0_AP.nii.gz -axis 3' 
                                     ]) 
    os.system(cmd_b0_extract) 
    flirt = ' '.join(['flirt -dof 6 -cost normmi -in',
                                     path_to_dwi_mif+'/mean_b0_AP.nii.gz -ref',                                
                                     path_to_ni_files+'/t1cBrainExtractionBrain.nii -omat',
                                     path_to_transform_matrix+'/T_fsl.txt' 
                                     ]) 
    os.system(flirt) 
    transform_convert = ' '.join(['transformconvert',
                                     path_to_transform_matrix+'/T_fsl.txt',                                
                                     path_to_dwi_mif+'/mean_b0_AP.nii.gz',
                                     path_to_ni_files+'/t1cBrainExtractionBrain.nii flirt_import',
                                     path_to_transform_matrix+'/T_dwi_to_t1c.txt && rm',
                                     path_to_transform_matrix+'/T_fsl.txt'  
                                     ]) 
    os.system(transform_convert) 
    mr_transform = ' '.join(['mrtransform -linear',
                                     path_to_transform_matrix+'/T_dwi_to_t1c.txt',                                
                                     path_to_dwi_mif+'/dwi_den_unr_pre_unbia_up.mif',
                                     path_to_dwi_mif+'/dwi_den_unr_pre_unbia_up_to_t1c.mif',
                                     ]) 
    os.system(mr_transform) 

# mask estimation using RSA stide and upscaled preprocessed registered DWI.
def dwi_mask_estimation(path_to_dwi_mif):
    mask_estimation = ' '.join(['dwi2mask',
                                path_to_dwi_mif+'/dwi_den_unr_pre_unbia_up_to_t1c.mif',
                                path_to_dwi_mif+'/mask_up_to_t1c.nii.gz' 
                                ])
    #mask_nifti = ' '.join(['mrconvert',
    #                            path_to_dwi_mif+'/mask_up_to_t1c.mif',
    #                            path_to_dwi_mif+'/mask_up_to_t1c.nii.gz -strides 1,2,3' 
    #                            ])
    os.system(mask_estimation) 

    
def fod_estimation(path_to_dwi_mif):
    response = ' '.join(['dwi2response dhollander', 
                        path_to_dwi_mif+'/dwi_den_unr_pre_unbia_up_to_t1c.mif',
                        path_to_dwi_mif+'/wm_to_t1c.txt',
                        path_to_dwi_mif+'/gm_to_t1c.txt',
                        path_to_dwi_mif+'/csf_to_t1c.txt' 
                        ])
    fod_estimation = ' '.join(['dwi2fod msmt_csd',
                            path_to_dwi_mif + '/dwi_den_unr_pre_unbia_up_to_t1c.mif', 
                            path_to_dwi_mif + '/wm_to_t1c.txt',
                            path_to_dwi_mif + '/wm_to_t1c.nii.gz',
                            '-mask',
                            path_to_dwi_mif + '/mask_up_to_t1c.nii.gz'
                            ])
    os.system(response+' && '+fod_estimation)

    # !mtnormalise /5_dwi/wm.mif /5_dwi/wm_norm.mif /5_dwi/csf.mif /5_dwi/csf_norm.mif -mask /5_dwi/mask_up.mif;


def export_bvecs_bvals(path_to_dwi_mif):
    # export bvecs and bvals
    # generate the ".nii" files of ".mif" file
    bvecs_bvals = ' '.join(['mrconvert', 
                            path_to_dwi_mif + '/dwi_den_unr_pre_unbia_up_to_t1c.mif',
                            path_to_dwi_mif + '/dwi_den_unr_pre_unbia_up_to_t1c.nii.gz',
                            '-export_grad_fsl',
                            path_to_dwi_mif + '/dwi_den_unr_pre_unbia_up_to_t1c.bvecs',
                            path_to_dwi_mif + '/dwi_den_unr_pre_unbia_up_to_t1c.bvals -stride 1,2,3,4'
                            ])

    os.system(bvecs_bvals)

def generate_peaks(path_to_dwi_mif):
    # generate peaks
    peaks = ' '.join(['sh2peaks',
                      path_to_dwi_mif + '/wm_to_t1c.nii.gz',
                      path_to_dwi_mif + '/peaks_to_t1c.nii.gz' 
                      ]) 
    os.system(peaks)


def tract_seg(path_to_dwi_mif):
    # Tractseg
    # bvecs and bvals should be in the same directory as DWI
    os.system('TractSeg -i ' + path_to_dwi_mif + '/dwi_den_unr_pre_unbia_up_to_t1c.nii.gz --raw_diffusion_input')
    os.system('TractSeg -i ' + path_to_dwi_mif + '/peaks.nii.gz --output_type tract_segmentation')
    os.system('TractSeg -i ' + path_to_dwi_mif + '/peaks.nii.gz --output_type endings_segmentation')
    os.system('TractSeg -i ' + path_to_dwi_mif + '/peaks.nii.gz --output_type TOM')
    os.system('Tracking -i ' + path_to_dwi_mif + '/peaks.nii.gz --bundles CST_right, CST_left')




#def registration_to_t1(path_dwi, path_t1, path_t1_mask):


    # create b0
    # select_dwi_vols ${dwi} ${bvals} nodif.nii.gz -b 1000
#   create_b0 = 'select_dwi_vols ' +  path_dwi + '/dwi_den_unr_unbis_up.nii.gz ' + path_dwi + 'dwi.bvals ' + path_dwi + 'nodif.nii.gz -b 1000'
    
#    convert_4D_3D = 'ExtractSliceFromImage 4 ' +  path_dwi + '/nodif.nii.gz ' + path_dwi + 'nodif_mean.nii.gz 3 0'

    # convert 4D to 3D
    # ExtractSliceFromImage 4 b0.nii b0_volume0.nii.gz 3 0
    #masked_img(path_dwi+'/nodif_mean.nii.gz') # output is nodif_brain
 

#   reg = 'epi_reg --epi=' + path_dwi + 'nodif_brain.nii.gz --t=' +  path_t1 + '--t1brain=' + path_t1_mask

    # epi_reg --epi=nodif_brain --t1=${t1} --t1brain=t1_brain --out=nodif_acpc # use existing white matter segmentation --wmseg=t1brain_wmseg

    # flirt -ref nodif_acpc -in ${dwi} -applyxfm -init nodif_acpc.mat -out dwi_aligned.nii.gz


    





