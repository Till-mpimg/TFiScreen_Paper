#!/usr/bin/env python3

#Read packages
import os
import re
import scipy.stats as stats
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, colorConverter, LinearSegmentedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from skimage import io, measure, util, restoration, exposure
from skimage.feature import peak_local_max
from skimage.morphology import closing, disk
from skimage.segmentation import clear_border, watershed
from skimage.filters import threshold_otsu, gaussian, median, threshold_niblack
from skimage.measure import label, regionprops_table, find_contours
from scipy.ndimage import distance_transform_edt
from aicsimageio import AICSImage
import seaborn as sns



#Function to turn hex codes into rgb
def hex_to_RGB(hex_str):
    """ #FFFFFF -> [255,255,255]"""
    #Pass 16 to the integer function for change of base
    return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

#Function to create custom color maps
def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1))/255
    c2_rgb = np.array(hex_to_RGB(c2))/255
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]

#Function to segment max projections of dapi-stained nuclei
def segment_dapi(dapi_max):
    # Segment DAPI image with Otsu threshold
	dapi_contrast = exposure.equalize_adapthist(dapi_max)
	
	thresh = threshold_otsu(dapi_contrast)
	dapi_otsu = dapi_contrast > thresh

    # Apply watershed to separate cells that touch
    # First find local distances
	distance = distance_transform_edt(dapi_otsu)

    # Find local maxima and return mask
	local_maxima = peak_local_max(distance, min_distance=20)
	peaks_mask = np.zeros_like(distance, dtype=bool)
	peaks_mask[tuple(local_maxima.T)] = True
	markers = label(peaks_mask)

    # Segment with watershed method according to local maxima
	dapi_segmented = watershed(-distance, markers, mask=dapi_otsu)
	
	for region_label in np.unique(dapi_segmented):
		# Create a mask for the current region
		region_mask = (dapi_segmented == region_label)
		# Close holes in the region using closing operation
		closed_region = closing(region_mask, disk(15))
    
		# Update the labeled image with the closed region
		dapi_segmented[region_mask] = 0  # Clear the original region
		dapi_segmented[closed_region] = region_label  # Assign the closed region
    
    # Remove cells bordering the window
	dapi_cleared = clear_border(dapi_segmented)

    # Get table with area and eccentricity of the individual cells
	table = regionprops_table(dapi_cleared, properties=('label', 'area', 'eccentricity'))
	condition = (table['area'] > 1500) & (table['eccentricity'] < 0.8) & (table['area'] < 10000)

    # Zero out labels not meeting the condition
	input_labels = table['label']
	output_labels = input_labels * condition

	dapi_filtered = util.map_array(dapi_segmented, input_labels, output_labels)
	dapi_binary = dapi_filtered > 0

	return dapi_binary, dapi_otsu, dapi_filtered, dapi_segmented

#Function to segment cy5 signal (Xist clouds)
def segment_cy5(cy5_max, dapi_binary):
    # Apply Gaussian Blur to cy5 image
    cy5_gaussian = gaussian(cy5_max, 2)

    # Segment Xist clouds locally with niblack threshold
    cy5_binary = closing(cy5_gaussian > threshold_niblack(cy5_gaussian, 151, k=-3.6), disk(5))

    # Get table with area of the individual clouds
    table = regionprops_table(label(cy5_binary), properties=('label', 'area'))
    condition = (table['area'] > 15) & (table['area'] < 3000)

    # Zero out labels not meeting the condition
    input_labels = table['label']
    output_labels = input_labels * condition

    cy5_filtered = util.map_array(label(cy5_binary), input_labels, output_labels)

    # Remove signal outside of the cells
    cy5_dapi = cy5_filtered * dapi_binary
    cy5_mask = cy5_dapi > 0

    return cy5_dapi, cy5_mask

#Function to loop through cells and return xist cloud properties
def process_cell(row, dapi_filtered, cy5_dapi, cy5_max):
	dapi_area = row['area']
	cell_id = row['label']

    # Select the single cell as a mask
	cell_mask = np.zeros(dapi_filtered.shape)
	cell_mask[:, :][dapi_filtered == cell_id] = 1

    # Check if Xist cloud intersects with the cell
	xist_ol = cell_mask * cy5_dapi
	xist_check = xist_ol[xist_ol != 0]
	out_cloud = pd.DataFrame()
	out_cloud['cell_id'] = [cell_id]
	out_cloud['dapi_area'] = [dapi_area]

	if xist_check.size != 0:
		xist_signal = True
		out_cloud['xist_signal'] = [xist_signal]
        
        
        # If there is a Xist signal, return properties of the background cy5 signal
        # I changed this to only take the biggest area
		bg_ol = cell_mask * np.invert(cy5_dapi)
		bg_label = pd.DataFrame.from_dict(regionprops_table(label(bg_ol), cy5_max, properties=('label', 'area', 'intensity_mean')))
		bg_max_area = bg_label.loc[bg_label['area'].idxmax()]
		out_cloud['bg_mean'] = bg_max_area['intensity_mean']
		out_cloud['bg_area'] = bg_max_area['area']
		out_cloud['bg_sum'] = out_cloud['bg_mean'] * out_cloud['bg_area']


        # Return properties of the cy5 signal
		cy5_label = pd.DataFrame.from_dict(regionprops_table(label(xist_ol), cy5_max, properties=('label', 'area', 'intensity_mean')))
		cy5_label = cy5_label.rename(columns={'label': 'cy5_id', 'area': 'cy5_area', 'intensity_mean': 'cy5_mean'})
		cy5_label['cell_id'] = cell_id
		out_cloud = pd.merge(out_cloud, cy5_label, on='cell_id', how='left')
		out_cloud['cy5_sum'] = out_cloud['cy5_area'] * out_cloud['cy5_mean']
		out_cloud['cy5_norm_sub'] = (out_cloud['cy5_mean'] - out_cloud['bg_mean']) * out_cloud['cy5_area']
		out_cloud['cy5_norm_fc'] = (out_cloud['cy5_mean'] / out_cloud['bg_mean']) * out_cloud['cy5_area']
		out_cloud['n_xist'] = len(cy5_label)

	else:
		out_cloud.loc[0, 'xist_signal'] = False
		out_cloud.loc[0, 'bg_area'] = None
		out_cloud.loc[0, 'bg_mean'] = None
		out_cloud.loc[0, 'bg_sum'] = None
		out_cloud.loc[0, 'cy5_id'] = None
		out_cloud.loc[0, 'cy5_area'] = None
		out_cloud.loc[0, 'cy5_mean'] = None
		out_cloud.loc[0, 'cy5_sum'] = None
		out_cloud.loc[0, 'cy5_norm_sub'] = None
		out_cloud.loc[0, 'cy5_norm_fc'] = None
		out_cloud.loc[0, 'n_xist'] = 0
	
	return out_cloud


#Functions to slice images for plotting
def remove_black_border(a): 
    # Mask of non-zeros
    mask = a!=0 # Use a >tolerance for a tolerance defining black border

    # Mask of non-zero rows and columns
    mask_row = mask.any(1)
    mask_col = mask.any(0)

    # First, last indices among the non-zero rows
    sr0,sr1 = mask_row.argmax(), len(mask_row) - mask_row[::-1].argmax()

    # First, last indices among the non-zero columns
    sc0,sc1 = mask_col.argmax(), len(mask_col) - mask_col[::-1].argmax()

    # Finally slice along the rows & cols with the start and stop indices to get 
    # cropped image. Slicing helps for an efficient operation.
    return a[sr0:sr1, sc0:sc1]

#Function to plot DAPI masks
def plot_segmentation_steps(dapi_max, dapi_otsu, dapi_segmented, dapi_filtered, pixel_length, title_suffix, zoom_in=False):

	
    #Adds scalebar to the plot
    target_length = 10 if zoom_in else 20
    scalebar = ScaleBar(pixel_length, "µm", location='lower right', color='white', frameon=False, fixed_value = target_length, scale_loc='top')
    scalebar._font_properties.set_family('Arial')
    scalebar._font_properties.set_size(6)

    fig, ax = plt.subplots(2, 2)
    ax[0, 0].imshow(dapi_max, cmap=plt.cm.gray)
    ax[0, 0].set_title("Max DAPI")

    ax[0, 1].imshow(dapi_otsu, cmap=plt.cm.gray)
    ax[0, 1].set_title("DAPI Otsu")

    ax[1, 0].imshow(dapi_segmented, cmap=Rainbow)
    ax[1, 0].set_title("DAPI Watershed")

    ax[1, 1].imshow(dapi_filtered, cmap=Rainbow)
    ax[1, 1].set_title("DAPI Final")

    ax[0, 0].add_artist(scalebar)

    if zoom_in:
        for i in range(2):
            for j in range(2):
                ax[i, j].set_xlim(100, 400)
                ax[i, j].set_ylim(100, 400)

    fig.tight_layout()
    filename = os.path.join(wd, 'output_files', f"{sample}_scene_{scene}_dapi_{title_suffix}.pdf")
    fig.savefig(filename)
    plt.close('all')

#Function to plot masks with segmented Cy5 signal
def plot_masks_with_cy5(dapi_max, cy5_max, dapi_binary, cy5_mask, dapi_filtered, pixel_length, title_suffix, zoom_in=False):
    #Adds scalebar to the plot
	target_length = 10 if zoom_in else 20
	scalebar = ScaleBar(pixel_length, "µm", location='lower right', color='white', frameon=False, fixed_value = target_length, scale_loc='top')
	scalebar._font_properties.set_family('Arial')
	scalebar._font_properties.set_size(6)

    # Get unique labels for cells (to draw contours)
	regions = np.unique(dapi_filtered)

	fig, ax = plt.subplots(2, 2)
	ax[0, 0].imshow(dapi_max, cmap=plt.cm.gray)
	ax[0, 0].set_title("Max DAPI")

	ax[0, 1].imshow(cy5_max, cmap=plt.cm.gray, vmin = np.percentile(cy5_max, 50), vmax = np.percentile(cy5_max, 99.9))
	ax[0, 1].set_title("Max Cy5")

	ax[1, 0].imshow(dapi_max, cmap=plt.cm.gray)
	ax[1, 0].imshow(cy5_max, cmap=cmap_rb)
	ax[1, 0].set_title("Max Combi")

	ax[1, 1].imshow(dapi_binary, cmap=plt.cm.gray)
	ax[1, 1].imshow(cy5_mask, cmap=cmap_rb, vmax=1)
	ax[1, 1].set_title("Segmentation")

	ax[0, 0].add_artist(scalebar)

	for region in regions[1:]:
		mask = dapi_filtered == region
		contours = find_contours(mask, 0.5)
		for contour in contours:
			ax[0, 1].plot(contour[:, 1], contour[:, 0], linewidth=0.25, color='white', linestyle='dashed')
			ax[1, 0].plot(contour[:, 1], contour[:, 0], linewidth=0.25, color='white', linestyle='dashed')
			ax[1, 1].plot(contour[:, 1], contour[:, 0], linewidth=0.25, color='red', linestyle='dashed')

	if zoom_in:
		for i in range(2):
			for j in range(2):
				ax[i, j].set_xlim(100, 500)
				ax[i, j].set_ylim(200, 600)

	fig.tight_layout()
	filename = os.path.join(wd, 'output_files', f"{sample}_scene_{scene}_mask_{title_suffix}.pdf")
	fig.savefig(filename)
	plt.close('all')


#Function to process each scene
def process_scene(czi, scene, pixel_length):
	channel_names = czi.channel_names
	dapi_pos = channel_names.index('DAPI')
	cy5_pos = channel_names.index('Cy5')
	
	czi.current_scene
	print(scene)
	czi.set_scene(scene)
	dapi = czi.data[0, dapi_pos, :, :, :] 
	cy5 = czi.data[0, cy5_pos, :, :, :] 

	dapi_max = remove_black_border(np.max(dapi, axis=0))
	cy5_max = remove_black_border(np.max(cy5, axis=0))

	dapi_binary, dapi_cleared, dapi_filtered, dapi_segmented = segment_dapi(dapi_max)

	plot_segmentation_steps(dapi_max, dapi_cleared, dapi_segmented, dapi_filtered, pixel_length, "total")

	cy5_dapi, cy5_mask = segment_cy5(cy5_max, dapi_binary)

	plot_masks_with_cy5(dapi_max, cy5_max, dapi_binary, cy5_mask, dapi_filtered, pixel_length, "total")

	dapi_labels = regionprops_table(dapi_filtered, properties=('label', 'area'))
	dapi_df = pd.DataFrame.from_dict(dapi_labels)
        
	if dapi_df.empty:
		results = None
	else:
		results = pd.concat([process_cell(row, dapi_filtered, cy5_dapi, cy5_max) for _, row in dapi_df.iterrows()], ignore_index=True)
	
	return results

#Creating Colormaps for visualization
Rainbow = ListedColormap(np.random.rand(256,3))
c_white = colorConverter.to_rgba('white',alpha = 0)
c_red= colorConverter.to_rgba('#90134D',alpha = 1)
cmap_rb = LinearSegmentedColormap.from_list('rb_cmap',[c_white,c_red],100)

#Work directory and input directory
wd = sys.argv[1]
input_dir = os.path.join(wd, 'input_files/CRISPRi_cont_FISH/')

# list all files in the folder
file_list = os.listdir(input_dir)
names = [i.split('.cz')[0] for i in file_list]

#Creates empty output dict and todays date for the output file
out_data = pd.DataFrame()

for sample in names:
    #Returns sample information
	file = f"{sample}.czi"
	
    #Read .czi file as AICSImage 
	czi = AICSImage(input_dir + sample + '.czi')

    #Read metadata and extract distance of single pixel
	pixel_length = float(czi.physical_pixel_sizes.Y)

    #Analyze each scene
	for scene in czi.scenes:
		scene_data = process_scene(czi, scene, pixel_length)

		if scene_data is not None:
			scene_data['scene'] = scene
			scene_data['sample'] = sample
			scene_data['file'] = file
			new_column_order = ['sample','file', 'scene'] + scene_data.columns[:-3].tolist()
			scene_data = scene_data[new_column_order]
			out_data = pd.concat([out_data, scene_data], ignore_index=True)


#Converts output to a pandas data frame and prints a txt file
out_data.to_csv(os.path.join(wd, 'output_files', 'CRISPRi_cont_analysis_out.txt'), index = False, sep = '\t')

#Counts Xist percentage and bi-allelic (and mis-segmented) cells
out_df_count = out_data[['sample', 'file', 'scene', 'cell_id', 'n_xist']].drop_duplicates()
out_df_count = out_df_count.groupby(['sample', 'n_xist'])['cell_id'].count().unstack(fill_value=0)

# Determine the number of columns for n_xist and the total amount of cells per sample
num_columns = len(out_df_count.columns)
total_cells = out_df_count.sum(axis=1)
out_perc_df = out_df_count.copy()

# Calculate xist_perc and percent of biallelic xist
print(out_df_count)
out_perc_df['xist_perc'] = out_df_count.loc[:, 1:].sum(axis=1) / total_cells
out_perc_df['mono_perc'] = out_df_count.iloc[:, 1] / total_cells
out_perc_df['bi_perc'] = out_df_count.iloc[:, 2] / total_cells

# Calculate percentage of cells with more than two Xist clouds (mis-segmented)
if num_columns > 3:
    out_perc_df['multi_perc'] = out_df_count.iloc[:, 3:].sum(axis=1) / total_cells


out_perc_df['total_cells'] = total_cells
out_perc_df.reset_index(inplace = True)

#Average Xist level
out_true = out_data[out_data['xist_signal'] == True]

out_intense = out_data.groupby(['sample']).apply(lambda group: group[['cy5_norm_sub', 'cy5_norm_fc']].apply(
	lambda col: stats.gmean(col, axis=0, nan_policy='omit')
	)).reset_index()
out_intense.columns = ['sample', 'cy5_norm_sub_geomean', 'cy5_norm_fc_geomean']

#Combine together and print plot
out_sum = out_perc_df.merge(out_intense)
out_sum.to_csv(os.path.join(wd, 'output_files', 'CRISPRi_cont_analysis_sum.txt'), index = False, sep = '\t')

#Plot Xist percentage and Xist level per cloud (with both background normalizations)
sns.set(style='whitegrid')
sns.barplot(x='sample', y='xist_perc', hue = 'sample', data = out_sum, alpha = 0.5, edgecolor = 'black', legend = False)
sns.barplot(x='sample', y='bi_perc', hue = 'sample', data = out_sum, alpha = 0.5, edgecolor = 'blue', legend = False)
sns.barplot(x='sample', y='multi_perc', hue = 'sample', data = out_sum, alpha = 0.5, edgecolor = 'red', legend = False)
sns.stripplot(x='sample', y='xist_perc', hue = 'sample', data = out_sum, dodge = True)
plt.ylim(0, 1)
plt.savefig(os.path.join(wd, 'output_files', '_CRISPRi_cont_xist_perc.pdf'))
plt.close('all')

sns.pointplot(data = out_true,  x = 'sample', hue = 'sample', y = 'cy5_norm_sub', dodge = 0.4, estimator = stats.gmean, errorbar = None,
                linestyles = "none", markers = "_", palette='dark:black', legend = False,  markersize=20, markeredgewidth=3, zorder = 10)
sns.swarmplot(data = out_true,  x = 'sample', hue = 'sample', y = 'cy5_norm_sub', dodge = 1, zorder = 1, size = 3)
plt.yscale('log')
plt.savefig(os.path.join(wd, 'output_files', '_CRISPRi_cont_cy5_norm_sub_swarm.pdf'))
plt.close('all')


sns.pointplot(data = out_true,  x = 'sample', hue = 'sample', y = 'cy5_norm_fc', dodge = 0.4, estimator = stats.gmean, errorbar = None,
                linestyles = "none", markers = "_", palette='dark:black', legend = False,  markersize=20, markeredgewidth=3, zorder = 10)
sns.swarmplot(data = out_true,  x = 'sample', hue = 'sample', y = 'cy5_norm_fc', dodge = 1, zorder = 1, size = 3)
sns.swarmplot(data = out_true,  x = 'sample', hue = 'sample', y = 'cy5_norm_fc', dodge = 1, zorder = 1, size = 3)
plt.yscale('log')
plt.savefig(os.path.join(wd, 'output_files', '_CRISPRi_cont__cy5_norm_fc_swarm.pdf'))
plt.close('all')
