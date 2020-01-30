import csv
import json

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        with open(snakemake.input['json'], 'r') as f:
            library_meta = json.load(f)[snakemake.wildcards['library']]
        
        columns = list(library_meta.keys())
        values = [library_meta[key] for key in columns]
        
        with open(snakemake.output['csv'], 'w') as csv_io:
            out_csv = csv.writer(csv_io, delimiter=',')
            out_csv = csv.writer(csv_io, delimiter=',')
            out_csv.writerow(["Barcode"] + columns)
            with open(snakemake.input['counts'], 'r') as f:
                for line in f:
                    row = [line.split('\t')[0]] + values
                    out_csv.writerow(row)