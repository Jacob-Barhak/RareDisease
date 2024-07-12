"""Load data for rare diseases for moco-makers meetup"""
import pandas as pd
import sqlite3
import os

def load_data_to_db():
    """Loads data from CSV into an sqlite data base"""
    # open an in memory database to do things quickly
    connection = sqlite3.connect(':memory:')
    # read the data and manipulate it and write it into the memory database
    # the response curve parameters has to be converted to multiple rows splitting the targets to be normalized
    response_curve_parameters = pd.read_csv("data" + os.sep + "secondary-screen-dose-response-curve-parameters.csv")
    response_curve_parameters['targets'] = response_curve_parameters['target']
    response_curve_parameters['target'] = response_curve_parameters['targets'].str.split(', ')
    response_curve_parameters = response_curve_parameters.explode('target')
    # the mutation matrix needs to be converted from a pivot table to a relational table
    response_curve_parameters.to_sql(name="response_curve_parameters", con=connection)
    # create indexes on the table for future use
    cursor = connection.cursor()
    for column in ['ccle_name', 'depmap_id', 'target', 'name']:
        sql = f"CREATE INDEX index_response_curve_parameters_{column} ON response_curve_parameters ({column});"
        cursor.execute(sql)
    del response_curve_parameters
    mutation_matrix_pivot = pd.read_csv("data" + os.sep + "OmicsSomaticMutationsMatrixDamaging.csv")
    mutation_matrix_pivot.rename(columns={'Unnamed: 0': 'depmap_id'}, inplace=True)
    mutation_matrix = mutation_matrix_pivot.melt(id_vars='depmap_id')
    mutation_matrix['target'] = mutation_matrix['variable'].str.split(' ').str[0]
    mutation_matrix.to_sql(name="mutation_matrix", con=connection)
    # create indexes on the table for future use
    for column in ['depmap_id', 'target', 'value']:
        sql = f"CREATE INDEX index_mutation_matrix_{column} ON mutation_matrix ({column});"
        cursor.execute(sql)
    # create a table as a select statement with only mutated numbers
    sql = "CREATE TABLE mutated AS " \
          "SELECT * FROM response_curve_parameters INNER JOIN mutation_matrix " \
          "ON response_curve_parameters.depmap_id = mutation_matrix.depmap_id " \
          "AND response_curve_parameters.target = mutation_matrix.target " \
          "WHERE mutation_matrix.value = 2"
    cursor.execute(sql)
    for column in ['depmap_id', 'target']:
        sql = f"CREATE INDEX index_mutated_{column} ON mutated ({column});"
        cursor.execute(sql)
    # copy the memory database to the output database
    out_db = sqlite3.connect('rare_disease.db')
    connection.backup(out_db)
    connection.close()
    out_db.close()


if __name__ == '__main__':
    load_data_to_db()
