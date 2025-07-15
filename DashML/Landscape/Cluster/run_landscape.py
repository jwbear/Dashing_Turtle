import os, sys, re
import pandas as pd
from DashML.Landscape.Cluster import Cluster_means
from DashML.Landscape.Cluster import Cluster_mode
from DashML.Landscape.Cluster import Cluster_native
from DashML.Landscape.Cluster import Centroid_MFE
from DashML.Landscape.Cluster import Centroid_Putative
from DashML.Landscape.Cluster import Native_MFE
from DashML.Landscape.Cluster import Cluster_Num as clust_num
from DashML.Landscape.Cluster import Centroid_Fold as centroid_prob
from DashML.Landscape.Cluster import Centroid_ConservedRegions as centroid_regions
from importlib.resources import files
from platformdirs import user_documents_path
from pathlib import Path
import DashML.Database_fx.Select_DB as dbsel

pd.options.mode.chained_assignment = None
print(user_documents_path())
output = user_documents_path() / "DTLandscape_Output"
output.mkdir(parents=True, exist_ok=True)


def save_path(subdir=""):
    save_path = output / "Figures" / subdir
    print(save_path)
    save_path.mkdir(parents=True, exist_ok=True)
    return str(save_path.resolve())+'/'


plot = True
native_lid = None
lids = None

### generate clusters, calculate centroids, remove short reads from clusters
#### from experimentally modified (predicted reactivities) data for kmeans/hamming
### opt plot dendrogram, heatmap, corrmatrix
def generate_clusters(lids, clust_num):
    #removes reads shorter than average
    df = dbsel.select_read_depth_full(lids)
    Cluster_means.save_path = save_path("Generated_Clusters/") # save figures
    img_dend, img_heatmap, img_corr = Cluster_means.den_array(df, clust_num=clust_num, plot=plot)
    Cluster_mode.save_path = save_path(subdir="Generated_Clusters/")  #save figures
    img_dendh, img_heatmaph, img_corrh = Cluster_mode.den_array(df, clust_num=clust_num, plot=plot)
    return img_dendh, img_heatmaph, img_corr


### generate distance matrices between native secondary structure and cluster centroids
### requires control structures, skip if no controls
### opt plot dendrogram, heatmap, corrmatrix
def get_native_structure_distances(lids, native_lid, plot=plot):
    Cluster_native.save_path= save_path(subdir="Native/")
    img_path = Cluster_native.den_array(lids, native_lid, plot=plot)
    return img_path

### generate structures of cluster centroids
### get mfe of cluster centroids (putative structures)
# TODO check that when creating centroids we have make sense, it seems okay but the cutoff in kmeans should be about .4
# to emulate hamming because modifications are rarely dominant
def generate_putative_structures(lids):
    Centroid_MFE.get_putative_structure(lids)

### generate varna images of predicted secondary structures from ViennaRNA
def draw_putative_structures(lids):
    Centroid_Putative.save_path= save_path(subdir="Putative_Structures/")
    Centroid_Putative.get_vplot(lids)


#calculate optimal number of clusters for a set of predictions
# save estimates, average over all sequences is used.
def get_optimal_clustnums(lids):
    # removes reads shorter than average
    df = dbsel.select_read_depth_full(lids)
    #print(df.head())
    clust_num.save_path = save_path(subdir="Cluster_Numbers/")  # save figures
    img_sil_k, img_wss_k = clust_num.den_array(df)
    return img_sil_k, img_wss_k

############# Generate Landscape of Reactivities
def run_landscape(temp=37,type1="dmso", type2="acim", complex=0, contig= "HCV", optimize_clusters=False):
    ### Extract Library IDs for controls (optionals if controls exist)
    native_lid = dbsel.select_lids(temp, type1, type1, complex, contig)
    native_lid = native_lid['ID'].unique()
    native_lid = re.sub("\[|\(|\]|\)", "", str(native_lid))

    ### get mfes for native secondary structures (controls)
    ### optional add lids for specific structures, default is all
    if native_lid is not None:
        if len(native_lid) > 0:
            Native_MFE.get_mfes(native_lid)

    ### Extract Library IDs
    lids = dbsel.select_lids(temp, type1, type2, complex, contig)
    lids = lids['ID'].unique()
    lids = re.sub("\[|\(|\]|\)", "", str(lids))

    ### Analysis Only: uses optimal numver of clusters###
    if optimize_clusters==True:
        img_sil_k, img_wss_k = get_optimal_clustnums(lids)
        return

    ### generate clusters
    generate_clusters(lids, clust_num=20)

    #### calculate controls/native if available
    if native_lid is not None:
        if len(native_lid) > 0:
            img_corr_means = get_native_structure_distances(lids, native_lid, plot=plot)

    ##### generate putative secondary structures
    generate_putative_structures(lids)

    ##### draw putative structures of generated secondary structures
    draw_putative_structures(lids)

    ### get interactions probabilities between centroids
    ### calculate possible interactions between alternate foldings
    centroid_prob.get_probabilities(lids)

    ##### calulate conserved regions
    centroid_regions.save_path = save_path(subdir="Conserved_Regions/")
    centroid_regions.get_conserved_regions(lids, region_size=3, gap=1, centroid_threshold=None)

    return


############# Generate Landscape of Reactivities
def run_landscape(lid, unmod_lid=None, optimize_clusters=False):
    try:
        ## image variables
        images = native_lid = img_heatmaph = img_dendh = img_corr_means = img_wss_k = None

        ### Extract Library IDs for controls (optionals if controls exist)
        if unmod_lid is not None:
            native_lid = dbsel.select_librarybyid(unmod_lid)
            #print(native_lid)
            if native_lid['secondary'].isnull().any():
                native_lid=None
            else:
                native_lid = native_lid['ID'].unique()
                native_lid = re.sub("\[|\(|\]|\)", "", str(native_lid))

        ### get mfes for native secondary structures (controls)
        ### optional add lids for specific structures, default is all
        if native_lid is not None:
            if len(native_lid) > 0:
                Native_MFE.get_mfes(native_lid)

        ### Extract Library IDs
        lids = dbsel.select_librarybyid(lid)
        lids = lids['ID'].unique()
        lids = re.sub("\[|\(|\]|\)", "", str(lids))

        ### Analysis Only: get optimal numver of clusters###
        if optimize_clusters==True:
            img_sil_k, img_wss_k = get_optimal_clustnums(lids)

        ### generate clusters
        img_dendh, img_heatmaph, img_corrh = generate_clusters(lids, clust_num=20)

        #### calculate controls/native if available
        if native_lid is not None:
            if len(native_lid) > 0:
                img_corr_means = get_native_structure_distances(lids, native_lid, plot=plot)

        ##### generate putative secondary structures
        generate_putative_structures(lids)

        ##### draw putative structures of generated secondary structures
        draw_putative_structures(lids)

        ### get interactions probabilities between centroids
        ### calculate possible interactions between alternate foldings
        centroid_prob.get_probabilities(lids)

        ##### calulate conserved regions
        centroid_regions.save_path = save_path("Conserved_Regions/")
        centroid_regions.get_conserved_regions(lids, region_size=3, gap=1, centroid_threshold=None)

        ### return sample images
        images = [img_heatmaph, img_dendh, img_corrh,img_wss_k]
        for i, img in enumerate(images):
            if img is None:
                if i == 3:
                    images[i] = "default8.png"
                else:
                    images[i] = "default7.png"

        #print(images)
        return images
    except Exception as err:
        raise Exception(str(err))
        return None
    finally:
        print("Landscape Analysis Complete")
        # TODO cleanup temp files


### optitmize clusters is only for getting cluster number
#run_landscape(unmod_lid=52, lid=37, optimize_clusters=True)

# if __name__ == "__main__":
#     print("Starting main......")
#     i = int(sys.argv[1])  # get the value of the $SLURM_ARRAY_TASK_ID
#     sequences = ['RNAse_P',"cen_3'utr", "cen_3'utr_complex", 'cen_FL', 'cen_FL_complex',
#                  "ik2_3'utr_complex", 'ik2_FL_complex', 'T_thermophila', 'ik2_FL', 'HCV',
#                  "ik2_3'utr"]
#     print(sequences[i])
#     run_landscape(sequences[i], temp, type1, type2, complex)



#sys.exit(0)
