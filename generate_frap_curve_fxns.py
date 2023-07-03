import numpy as np
import matplotlib.pyplot as plt
import skimage
from skimage import io, filters, draw
from skimage.measure import label, regionprops, regionprops_table
from matplotlib.patches import Circle
import os
import sys



def generate_frap_curves(read_dir, save_dir, time_array, frap_frame, radius, saveoutput):
    ###########################################
    
    #this function takes in a set of droplet FRAP videos and generates normalized frap curves from them
    #inputs:
        #read_dir    - where your individual droplet videos are saved (nothing else should be in this folder!)
        #save_dir    - directory where to save stuff - burn spots will be saved in a subfolder called './burn_spots_python/'
        #time_array  - timekeeping array of frap videos - units are up to you
        #frap_frame  - the frame at which the frap event first is present (in 0 based counting)
        #radius      - radius of the frap spot in pixels (input for get_burn_center)
        #saveoutput  - whether to save the burn spot output from get_burn_center

     #outputs:
        #image output of get_burn_center, which is saved in location specified above
            #shows raw image, burn spot pixels, burn spot centroid, and FRAPped droplet pixels with FRAP spot
        #normalized_intensities - n-droplets by t-timepoints of normalized intensities
        
    ###########################################
    
    isExist = os.path.exists(save_dir)
    if not isExist:
        os.makedirs(save_dir)
    else:
        quer = 0
        while quer < 1:
            user_input = input('Save directory already exists, files may be overwritten. Proceed? [y/n]')
            if user_input.lower() == 'y':
                print('Proceeding')
                quer = 1
            elif user_input.lower() == 'n':
                print('Aborting on user command')
                sys.exit()
            else:
                print('Please answer with [y] or [n]')
                quer = 0
         
    files = os.listdir(read_dir)
    normalized_intensities = np.zeros((len(files), len(time_array)))

    count = 0
    for file in files:
        filename = file.split('.')
        filename = filename[0]
        I = io.imread(read_dir + file)
        save_name = filename +'_burnspots.png'
        [yw,xw] = get_burn_center(I, frap_frame, radius, saveoutput, save_dir + './burn_spots_python', save_name)
        avg_intensity_inside = np.zeros(len(I-1))
        avg_intensity_outside = np.zeros(len(I-1))

        for t in np.arange(0,len(I-1)):
            #First, get droplet boundaries at time t

            # Computing Otsu's thresholding value
            threshold = filters.threshold_otsu(I[t])

            # Computing binarized values using the obtained threshold
            Ithresh = (I[frap_frame-1] > 1.25*threshold)*1

            #Now we find boundaries - largest boundary should correspond to our droplet of interest
            label_img = label(Ithresh)
            regions = regionprops(label_img)

            #note: if your droplet of interest for whatever reason is NOT the largest region, crop image to make this true
            areas = [i.area for i in regions]
            roi = areas.index(max(areas))

            #make mask with just our biggest drop included
            mask = np.zeros(np.shape(I[frap_frame]))
            mask[regions[roi].coords[:,0], regions[roi].coords[:,1]] = 1
            masked_droplet = mask*I[t];

            #get values within the burn spot radius of image at time t
            meanSpot = mean_circle(I[t], [xw,yw], radius)

            #now want to get pixels inside the droplet boundary but outside spot_radius
            h, w = np.shape(I[t])
            Y, X = np.ogrid[:h, :w]
            dR = np.sqrt((X - xw)**2 + (Y-yw)**2)
            maskspot = dR > radius
            final_mask = np.logical_and(masked_droplet, maskspot)
            cframe = I[t]
            meanOutside = np.mean(cframe[final_mask]);
            avg_intensity_inside[t] = meanSpot
            avg_intensity_outside[t] = meanOutside


        #renormalize intensities inside burn spot relative to outside the burn spot
        #this is based on normalization from Taylor, N.O., et al. (2019) Biophys. J. DOI: 10.1016/j.bpj.2019.08.030
        renorm = (avg_intensity_inside - avg_intensity_inside[frap_frame])/(avg_intensity_outside - avg_intensity_inside[frap_frame])
        normalized_intensities[count][:] = renorm 

        print('Completed: ' + file)
        count += 1
    
    return(normalized_intensities)


def get_burn_center(I, frap_frame, frap_radius, saveoutput, save_dir, save_name):
    ###########################################
    
    #this function fetches the burn center of a condensate frap timecourse
    #inputs:
        #Istack      - image stack, loaded with skimage.io.imread
        #frap_frame  - the frame at which the frap event first is present (in 0 based counting)
        #frap_radius - radius of the frap spot in pixels
        #saveoutput  - whether to save the burn spot output - this is highly recommended
        #save_dir    - directory where to save 
        #save_name   - name of save file

     #outputs:
        #image       - which is saved in location specified above
        #              shows raw image, burn spot pixels, burn spot centroid,
        #              and FRAPped droplet pixels with FRAP spot
        #center      - array with the burn spot weighted centroid location: [yw,xw]
        
    ###########################################


    #pixels which were particularly affected by the frap event
    frap_difference = I[frap_frame-1]-I[frap_frame]

    #exclude pixels outside droplet --> need to identify droplet boundary

    # Computing Otsu's thresholding value
    threshold = filters.threshold_otsu(I[frap_frame-1])

    # Computing binarized values using the obtained threshold
    Ithresh = (I[frap_frame-1] > 1.25*threshold)*1

    #Now we find boundaries - largest boundary should correspond to our droplet of interest
    label_img = label(Ithresh)
    regions = regionprops(label_img)

    #note: if your droplet of interest for whatever reason is NOT the largest region,
    #      crop image to make this true
    areas = [i.area for i in regions]
    roi = areas.index(max(areas))

    #make mask with just our biggest drop included
    mask = np.zeros(np.shape(I[frap_frame]))
    mask[regions[roi].coords[:,0], regions[roi].coords[:,1]] = 1

    masked_diff = mask*frap_difference

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=4, sharex=True,
                                        figsize=(12, 4))

    ax0.set_aspect('equal')
    ax0.imshow(I[frap_frame])
    ax0.set_title('Raw Image')

    ax1.set_aspect('equal')
    ax1.imshow(masked_diff)
    ax1.set_title('Droplet Pixels-Frap Frame')
    
    #Now from the pixels most affected by the FRAP, want to identify the actual burn spot
    #So another round of segmentation

    # Computing Yen's thresholding value
    thresh_in = filters.threshold_yen(masked_diff)
    # Computing binarized values using the obtained threshold
    Ithresh_in = (masked_diff > 1.2*thresh_in)*1

    #Now we find boundaries - largest boundary should correspond to our actual burn spot
    label_spot = label(Ithresh_in)
    potential_spots = regionprops(label_spot, masked_diff)
    spot_areas = [i.area for i in potential_spots]
    spot_no = spot_areas.index(max(spot_areas))


    #make mask with just the burn spot
    mask2 = np.zeros(np.shape(I[frap_frame]))
    mask2[potential_spots[spot_no].coords[:,0], potential_spots[spot_no].coords[:,1]] = 1
    burn_spot = masked_diff*mask2

    #cool, now want the center of the burn spot - weighted center of the burned pixels
    burn_center_raw = potential_spots[spot_no].centroid
    burn_center_weighted = potential_spots[spot_no].centroid_weighted
    yr,xr = burn_center_raw
    yw,xw = burn_center_weighted

    #plots the pixels of burn spot with the identified spot center
    ax2.set_aspect('equal')
    ax2.imshow(burn_spot)
    ax2.set_title('Burn spot mask and center')
    circ_weight = Circle((xw,yw),frap_radius, facecolor='none', edgecolor='blue')
    circ_excl = Circle((xw,yw),1.5*frap_radius, facecolor='none', edgecolor='red')
    ax2.scatter(xw,yw, s=5,facecolor='red')
    ax2.add_patch(circ_weight)
    ax2.add_patch(circ_excl)

    #plots the pixels of the droplet overlayed on the raw image file,
    #     and includes the burn spot pixels
    ax3.set_aspect('equal')
    ax3.imshow(I[frap_frame], cmap='gray')
    ax3.imshow(mask, alpha=0.5)
    ax3.set_title('Droplet Pixels with Burn Spot')
    circ_weight = Circle((xw,yw),frap_radius, facecolor='none', edgecolor='blue')
    circ_excl = Circle((xw,yw),1.5*frap_radius, facecolor='none', edgecolor='red')
    ax3.scatter(xw,yw, s=5, facecolor='red')
    ax3.add_patch(circ_weight)
    ax3.add_patch(circ_excl)

    if saveoutput == 1:
        plt.savefig(save_dir+save_name)

    return [yw,xw]


def mean_circle(img, center, radius):
    #pretty simple, just gets the average of all pixels within some radius of a center
    #inputs:
        #img    - 2d image that you want this performed on
        #center - coordinates (x,y) of center you want calculated on
        #radius - radius of circle for the operation
    #outputs:
        #mean   - mean intensity of the circular region
    
    h, w = np.shape(img)
    Y, X = np.ogrid[:h, :w]
    dR = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dR <= radius
    cmean = np.mean(img[mask])
    return(cmean)