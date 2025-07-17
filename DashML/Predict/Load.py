import sys
import Predict_Fold as pbpp
#import Metric as mx
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel

#### Functions to preprocess/load data before predicting
#### review



## DESCRIPTION
### get predicted unmodified secondary structures for controls for comparison against known structure and
### other modified predictions, how far does vienna differ on unmodified sequence
### should match the non-psuedoknotted structure
### get predicted unmodified secondary structures for controls
### read index == -2 indicates it's a prediction on unmodified data
### saves to structure_secondary_interactions as base probabilities are recorded
def get_predicted_structures():
    # no reactivity for unmodified and averaged if unmodified reactivity not calculated
    # uses unmodified lids

    #get unmod lids
    df = dbsel.select_unmod_ssi()
    lids = df['LID'].unique()

    for lid in lids:
        print(lid)
        pbpp.get_probabilities(lid, lid, reactivity=None, read=-2)

get_predicted_structures()
