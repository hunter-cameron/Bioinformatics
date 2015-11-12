
import argparse
import pandas
from mypyli import picrust




def run_full_correlation(pred_f, obs_f, metadata_f, pa_f, out_f):
    ko_categories = picrust.get_ko_by_function(metadata_f, level=2)
    plant_associated = picrust.get_plant_associated_kos(pa_f)

    ko_cat_names = sorted(ko_categories.keys())
    plant_lineages = sorted(plant_associated.keys())

    headers = ["Overall;Overall"] +  ko_cat_names + plant_lineages

    category_dict = ko_categories.copy()    
    category_dict.update(plant_associated)

    with open(out_f, 'w') as OUT:
        OUT.write("\t".join(["genome"] + ["Overall;Overall"] +  ko_cat_names + ["pa_" + p for p in plant_lineages] + ["NSTI"]) + "\n")

        pred_ttm = picrust.TraitTableManager(pred_f)
        obs_ttm = picrust.TraitTableManager(obs_f)


        obs_dict = {obs.name: obs for obs in obs_ttm}

        for pred in pred_ttm:
 
            try:
                obs = obs_dict[pred.name]
            except KeyError:
                continue

            try:
                nsti = pred.traits["metadata_NSTI"]
            except KeyError:
                nsti = None

            # calculate the correlation by group
            correl_dict = get_correlation_dict(obs, pred, category_dict)

            # calculate the overall correlation
            cor = obs.correlation(pred, traits=pred_ttm.traits)
            correl_dict["Overall;Overall"] = cor

            OUT.write("\t".join([obs.name] + [str(correl_dict[name]) for name in headers] + [str(nsti)]) + "\n")

def get_correlation_dict(obs, pred, category_dict):
    correl_dict = {}

    for header in category_dict:
        try:
            cor = obs.correlation(pred, traits=category_dict[header])
            
            # this is checking for NaN
            if cor > 0:
                correl_dict[header] = cor
            else:
                correl_dict[header] = "NA"

        except ValueError as e:
            print(e)
            correl_dict[header] = "NA"

    return correl_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculates correlations between PICRUSt tables")
    parser.add_argument("-pred", help="predicted table", required=True)
    parser.add_argument("-obs", help="observed table (this will be a subset of the predicted table)", required=True)
    parser.add_argument("-metadata", help="trait metadata file to parse pathways from")
    parser.add_argument("-pa", help="table of plant associated genes")
    parser.add_argument("-out", help="output path for the correlations")

    args = parser.parse_args()

    run_full_correlation(args.pred, args.obs, args.metadata, args.pa, args.out)    
