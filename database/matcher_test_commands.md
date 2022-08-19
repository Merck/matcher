# Matcher test dataset

## Dependencies

- rdkit
- psycopg2
- mmpdb
- scipy
- numpy

## Creating fragments
Fragment the input SMILES into a fragment set, no DB connection required.  

```
python ../mmpdb/mmpdb.py fragment matcher_test_structures.smi -o matcher_test_structures.fragments --cut-smarts 'matcher_alpha'
```

## Indexing fragments
Index the fragments to define the MMPs, writing all data to DB (except for property data).  
The final argument is DB connection info, having format `<server_address>$<port>$<db_name>$postgres`  
If necessary for connection to DB, set environment variables **MATCHER_USER** and **MATCHER_PASSWORD**
```
python ../mmpdb/mmpdb.py index matcher_test_structures.fragments -o 'localhost$5432$matcher_db_dev$postgres'
```

## Property data
Load property data for compounds to the DB, and calculate statistics
```
python ../mmpdb/mmpdb.py loadprops -p matcher_test_data.csv 'localhost$5432$matcher_db_dev$postgres'
```
