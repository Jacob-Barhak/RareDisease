"""Load data for rare diseases for moco-makers meetup"""
import pandas as pd
import sqlite3
import os
import logging


def load_data_to_db():
    """Loads data from CSV into an sqlite data base"""
    # open an in memory database to do things quickly
    connection = sqlite3.connect('rare_disease.db')
    cursor = connection.cursor()
    # load ligand data raw
    logging.info("creating ligand table")
    ligand = pd.read_csv("data" + os.sep + "family_target_ligand_all_gtopdb.csv")
    ligand.to_sql(name="family_target_ligand_raw", con=connection)
    del ligand
    # created a minimal version of this table to have only date we need without duplicates
    sql = \
        "CREATE TABLE family_target_ligand_min AS " \
        "SELECT ligand_name, GROUP_CONCAT(DISTINCT family_type) as family_types, " \
        "GROUP_CONCAT(DISTINCT family_name) as family_names " \
        "FROM (SELECT DISTINCT family_type, family_name, ligand_name FROM family_target_ligand_raw) " \
        "GROUP BY ligand_name"
    cursor.execute(sql)
    logging.info(sql)

    logging.info("loading approved compounds table")
    with open("data" + os.sep + "approved_compounds.txt") as file:
        text = file.read()
    diseases = text.split('Targeted therapy approved for ')[1:]
    records = []
    for disease in diseases:
        disease_data = [line for line in disease.split('\n') if line != '']
        disease_name = disease_data[0]
        compounds = disease_data[1:]
        for compound in compounds:
            [generic, brand] = (compound.rstrip(')') + ' (').split(' (')[:2]
            records.append([disease_name, generic, brand])

    approved_compounds = pd.DataFrame(records, columns=['disease', 'generic', 'brand'])
    approved_compounds.to_sql(name="approved_compounds_raw", con=connection)

    sql = \
        "CREATE TABLE approved_compounds AS " \
        "SELECT DISTINCT generic, brand " \
        "FROM approved_compounds_raw"
    cursor.execute(sql)
    logging.info(sql)

    # read the data and manipulate it and write it into the memory database
    # the response curve parameters has to be converted to multiple rows splitting the targets to be normalized
    logging.info("creating response curve table")
    response = pd.read_csv("data" + os.sep + "secondary-screen-dose-response-curve-parameters.csv")
    response.rename(columns={'disease.area': 'disease_area', 'target': 'targets'}, inplace=True)
    response['tissue'] = response['ccle_name'].str.split(pat='_', n=1).str[1]
    # the mutation matrix needs to be converted from a pivot table to a relational table
    response.to_sql(name="response_raw", con=connection)

    logging.info("joining raw response curve, ligand and approved compounds tables into a new "
                 "response table without failed records")
    sql = \
        "CREATE TABLE response AS " \
        "SELECT response_raw.*, " \
        "family_target_ligand_min.*, " \
        "ifnull(response_raw.name = approved_compounds.generic " \
        " OR response_raw.name = approved_compounds.brand, 0) as approved " \
        "FROM response_raw " \
        "LEFT JOIN family_target_ligand_min " \
        "ON response_raw.name = ligand_name " \
        "LEFT JOIN approved_compounds " \
        "ON response_raw.name = approved_compounds.generic " \
        "OR response_raw.name = approved_compounds.brand " \
        "WHERE response_raw.passed_str_profiling = 1"
    cursor.execute(sql)
    logging.info(sql)

    for column_to_drop in ['smiles', 'upper_limit', 'passed_str_profiling', 'row_name', 'phase']:
        logging.info(f"removed the redundant column {column_to_drop} from the response table")
        sql = \
            f"ALTER TABLE response DROP {column_to_drop}"
        cursor.execute(sql)
        logging.info(sql)

    # create indexes on the table for future use
    # define columns to index in the database for future faster processing
    index_column_names = \
        ['depmap_id', 'targets', 'moa',
         'tissue', 'screen_id', 'name',
         "family_types", "family_names", "approved"]
    for column in index_column_names:
        logging.info(f"indexing response curve table {column=}")
        sql = f"CREATE INDEX index_response_{column} ON response ({column});"
        cursor.execute(sql)
        logging.info(sql)

    # free memory
    del response
    logging.info("creating mutation matrix table as relational table")
    mutation_matrix_pivot = pd.read_csv("data" + os.sep + "OmicsSomaticMutationsMatrixDamaging.csv")
    mutation_matrix_pivot.rename(columns={'Unnamed: 0': 'depmap_id'}, inplace=True)
    mutation_matrix = mutation_matrix_pivot.melt(id_vars='depmap_id')
    mutation_matrix['gene'] = mutation_matrix['variable'].str.split(' ').str[0]
    mutation_matrix.to_sql(name="mutation_matrix_full", con=connection)
    # create indexes on the table for future use
    for column in ['depmap_id', 'gene', 'value']:
        logging.info(f"indexing mutation matrix table {column=}")
        sql = f"CREATE INDEX index_mutation_matrix_full_{column} ON mutation_matrix_full ({column});"
        cursor.execute(sql)
        logging.info(sql)

    # create a smaller mutation table with relevant records
    # free memory
    del mutation_matrix_pivot, mutation_matrix

    for mutation_level in [2]:
        logging.info(f"extracting targets that have mutation level {mutation_level}")
        sql = f"CREATE TABLE select_genes{mutation_level} AS " \
              f"SELECT DISTINCT gene from mutation_matrix_full WHERE value = {mutation_level} "
        cursor.execute(sql)
        logging.info(sql)

        for column in ['gene']:
            logging.info(f"index select_genes for mutation level {mutation_level} {column=}")
            sql = \
                f"CREATE INDEX index_select_genes{mutation_level}_{column} ON select_genes{mutation_level} ({column});"
            cursor.execute(sql)
            logging.info(sql)

        # create a smaller mutation table ignoring columns where the mutation level was not mutation_level
        logging.info(f"creating a mutation matrix that have mutation level {mutation_level}")
        sql = f"CREATE TABLE select_mutation_matrix{mutation_level} AS " \
              f"SELECT * FROM mutation_matrix_full inner join select_genes{mutation_level} " \
              f"WHERE mutation_matrix_full.gene = select_genes{mutation_level}.gene"
        cursor.execute(sql)
        logging.info(sql)

        for column in ['depmap_id', 'gene', 'value']:
            logging.info(f"indexing select_mutation_matrix {mutation_level} table {column=}")
            sql = f"CREATE INDEX index_mutation_matrix_{column} ON select_mutation_matrix{mutation_level} ({column});"
            cursor.execute(sql)
            logging.info(sql)

        logging.info(f"create combined mutation table for mutation level {mutation_level}")
        sql = f"CREATE TABLE mutated{mutation_level} AS " \
              f"SELECT * FROM response INNER JOIN select_mutation_matrix{mutation_level} " \
              f"ON response.depmap_id = select_mutation_matrix{mutation_level}.depmap_id "
        cursor.execute(sql)
        logging.info(sql)

        for column in index_column_names + ['gene', 'value']:
            logging.info(f"index combined mutation table for mutation level {mutation_level} {column=}")
            sql = f"CREATE INDEX index_mutated{mutation_level}_{column} ON mutated{mutation_level} ({column});"
            cursor.execute(sql)
            logging.info(sql)

        logging.info(f"create aggregate_keys for mutated{mutation_level}")

        for agg_columns in [
            ['screen_id', 'gene', 'tissue', 'moa', 'name', 'family_types', 'family_names', 'approved'],
        ]:
            agg_columns_comma_str = ', '.join(agg_columns)
            agg_columns_underscore_str = '_'.join(agg_columns)
            logging.info(
                f"create the aggregated mutation table for {agg_columns} for mutation level {mutation_level}")

            only_calculations = []
            target_calculations = []
            target_columns = []
            for agg_func in ['avg', 'min', 'max']:
                for agg_column in ['lower_limit', 'slope', 'r2', 'auc', 'ec50', 'ic50']:
                    target_calculations += [f'{agg_func}({agg_column}) as {agg_func}_{agg_column}']
                    only_calculations += [f'{agg_func}({agg_column})']
                    target_columns += [f'{agg_func}_{agg_column}']

            comma_separated_target_calculations = ', '.join(target_calculations)
            sql = f"CREATE TABLE aggregated_{agg_columns_underscore_str}_{mutation_level} AS " \
                  f"SELECT value, {agg_columns_comma_str}, count(value) as number_of_records, " \
                  f"{comma_separated_target_calculations} " \
                  f"FROM mutated{mutation_level} " \
                  f"GROUP BY value, {agg_columns_comma_str}"
            cursor.execute(sql)
            logging.info(sql)

            index_agg_column_names = \
                ['value'] + agg_columns
            for column in index_agg_column_names:
                logging.info(f"index aggregated mutation table "
                             f"aggregated_{agg_columns_underscore_str}_{mutation_level} "
                             f"for mutation level {mutation_level} {column}")
                sql = f"CREATE INDEX aggregated_{agg_columns_underscore_str}_{mutation_level}_{column} " \
                      f"ON aggregated_{agg_columns_underscore_str}_{mutation_level} ({column});"
                cursor.execute(sql)
                logging.info(sql)

            new_tables = {}
            for filter_mutation_level in [0, 1, 2]:
                new_table_name = \
                    f"select_aggregated_{agg_columns_underscore_str}_{mutation_level}_{filter_mutation_level}"
                new_tables[filter_mutation_level] = new_table_name
                sql = f"CREATE TABLE {new_table_name} AS " \
                      f"SELECT *" \
                      f"FROM aggregated_{agg_columns_underscore_str}_{mutation_level} " \
                      f"WHERE value = {filter_mutation_level}"
                cursor.execute(sql)
                logging.info(sql)

                index_agg_column_names = agg_columns
                for column in index_agg_column_names:
                    logging.info(f"index select_aggregated mutation table "
                                 f"{new_table_name} "
                                 f"for mutation level {mutation_level} {filter_mutation_level} {column}")
                    sql = f"CREATE INDEX {new_table_name}_{column} " \
                          f"ON {new_table_name} ({column});"
                    cursor.execute(sql)
                    logging.info(sql)

                logging.info(f"index select_aggregated mutation table "
                             f"{new_table_name} "
                             f"for mutation level {mutation_level} {filter_mutation_level} {agg_columns_comma_str}")
                sql = f"CREATE INDEX {new_table_name}_keys " \
                      f"ON {new_table_name} ({agg_columns_comma_str});"
                cursor.execute(sql)
                logging.info(sql)

            logging.info(
                f"create the diff aggregated mutation table for {agg_columns} for mutation level {mutation_level}")

            table_columns = {}
            table_columns_str = {}
            for filter_mutation_level in [0, 1, 2]:
                table_columns[filter_mutation_level] = \
                    [f"{new_tables[filter_mutation_level]}.{column}"
                     for column in target_columns + ["number_of_records"]]
                table_columns_str[filter_mutation_level] = \
                    ", ".join(table_columns[filter_mutation_level])

            diff_columns = [f"{column_table2} - {column_table0} as diff_{column_name_only}"
                            for (column_table2, column_table0, column_name_only) in
                            zip(table_columns[2][:-1], table_columns[0][:-1], target_columns)]
            comma_separated_diff_columns = ", ".join(diff_columns)

            join_clause = [f"{new_tables[2]}.{column} = {new_tables[0]}.{column}"
                           for column in agg_columns]
            join_clause_str = " AND ".join(join_clause)

            agg_table_and_columns = [f"{new_tables[2]}.{column}" for column in agg_columns]
            agg_table_and_columns_comma_str = ", ".join(agg_table_and_columns)

            new_table_name = f"diff_aggregated_{agg_columns_underscore_str}_{mutation_level}"
            sql = f"CREATE TABLE {new_table_name} AS " \
                  f"SELECT {agg_table_and_columns_comma_str}, " \
                  f"{table_columns_str[2]}, {table_columns_str[0]} , " \
                  f"{comma_separated_diff_columns} " \
                  f"FROM {new_tables[2]} JOIN {new_tables[0]} ON " \
                  f"{join_clause_str}"
            cursor.execute(sql)
            logging.info(sql)

            index_agg_column_names = agg_columns
            for column in index_agg_column_names:
                logging.info(f"index diff_aggregated mutation table "
                             f"{new_table_name} "
                             f"for mutation level {mutation_level} {column}")
                sql = f"CREATE INDEX {new_table_name}_{column} " \
                      f"ON {new_table_name} ({column});"
                cursor.execute(sql)
                logging.info(sql)

    connection.close()


if __name__ == '__main__':
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        filename="log.txt",
        level=logging.INFO)
    logger = logging.getLogger()
    stream = logging.StreamHandler()
    logger.addHandler(stream)

    load_data_to_db()
