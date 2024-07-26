"""Load data for rare diseases for moco-makers meetup"""
import pandas as pd
import sqlite3
import os
import logging

def load_data_to_db():
    """Loads data from CSV into an sqlite data base"""
    # open an in memory database to do things quickly
    connection = sqlite3.connect(':memory:')
    # read the data and manipulate it and write it into the memory database
    # the response curve parameters has to be converted to multiple rows splitting the targets to be normalized
    logging.info("creating response curve table")
    response = pd.read_csv("data" + os.sep + "secondary-screen-dose-response-curve-parameters.csv")
    response.rename(columns={'disease.area': 'disease_area'}, inplace=True)
    response['targets'] = response['target']
    response['target'] = response['targets'].str.split(', ')
    response = response.explode('target')
    # the mutation matrix needs to be converted from a pivot table to a relational table
    response.to_sql(name="response", con=connection)
    # create indexes on the table for future use
    cursor = connection.cursor()
    # define columns to index in the database for future faster processing
    index_column_names = \
        ['depmap_id', 'target', 'moa', 'disease_area', 'ccle_name', 'name', 'indication', 'phase', 'smiles']
    for column in index_column_names:
        logging.info(f"indexing response curve table {column=}")
        sql = f"CREATE INDEX index_response_{column} ON response ({column});"
        cursor.execute(sql)
    # free memory
    del response
    logging.info("creating mutation matrix table as relational table")
    mutation_matrix_pivot = pd.read_csv("data" + os.sep + "OmicsSomaticMutationsMatrixDamaging.csv")
    mutation_matrix_pivot.rename(columns={'Unnamed: 0': 'depmap_id'}, inplace=True)
    mutation_matrix = mutation_matrix_pivot.melt(id_vars='depmap_id')
    mutation_matrix['gene'] = mutation_matrix['variable'].str.split(' ').str[0]
    mutation_matrix.to_sql(name="mutation_matrix", con=connection)
    # create indexes on the table for future use
    for column in ['depmap_id', 'gene', 'value']:
        logging.info(f"indexing mutation matrix table {column=}")
        sql = f"CREATE INDEX index_mutation_matrix_{column} ON mutation_matrix ({column});"
        cursor.execute(sql)
    # create tables as a select statement with mutated numbers with mutation_level level of mutation of 0,1,2
    # free memory
    del mutation_matrix_pivot, mutation_matrix
    for mutation_level in [2, 0]:
        logging.info(f"create combined mutation table for mutation level {mutation_level}")

        # TODO: Refactor this:
        # Instead of an inner join, this might a many-to-many join of some kind
        # https://stackoverflow.com/questions/17774373/sql-join-many-to-many
        sql = f"CREATE TABLE mutated{mutation_level} AS " \
              f"SELECT * FROM response INNER JOIN mutation_matrix " \
              f"ON response.depmap_id = mutation_matrix.depmap_id " \
              f"WHERE mutation_matrix.value = {mutation_level}"
        cursor.execute(sql)
        for column in index_column_names:
            logging.info(f"index combined mutation table for mutation level {mutation_level} {column=}")
            sql = f"CREATE INDEX index_mutated{mutation_level}_{column} ON mutated{mutation_level} ({column});"
            cursor.execute(sql)
        logging.info(f"create aggregated mutation table for mutation level {mutation_level}")
        sql = f"CREATE TABLE aggregated{mutation_level} AS " \
              f"select depmap_id, ccle_name, name, moa, disease_area, indication, count(target), " \
              f"avg(upper_limit), avg(lower_limit), avg(slope), avg(r2), avg(auc), avg(ec50), avg(ic50), " \
              f"min(upper_limit), min(lower_limit), min(slope), min(r2), min(auc), min(ec50), min(ic50), " \
              f"max(upper_limit), max(lower_limit), max(slope), max(r2), max(auc), max(ec50), max(ic50) " \
              f"from mutated{mutation_level} " \
              f"group by depmap_id, ccle_name, name, moa, disease_area, indication"
        cursor.execute(sql)
        index_agg_column_names = \
            ['depmap_id', 'moa', 'disease_area', 'ccle_name', 'name', 'indication']
        for column in index_agg_column_names:
            logging.info(f"index aggregated mutation table for mutation level {mutation_level} {column=}")
            sql = f"CREATE INDEX aggregated{mutation_level}_{column} ON aggregated{mutation_level} ({column});"
            cursor.execute(sql)

    # copy the memory database to the output database
    out_db = sqlite3.connect('rare_disease.db')
    connection.backup(out_db)
    connection.close()
    out_db.close()


if __name__ == '__main__':
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO)
    load_data_to_db()
