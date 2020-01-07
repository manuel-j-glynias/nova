import sys
import os
import pymongo
import dns # required for connecting with SRV
from pymongo.errors import ConnectionFailure


def get_mongo_client():
    client = pymongo.MongoClient(
        "mongodb+srv://omni_user:700Ellicott@omniseq-unjsy.gcp.mongodb.net/test?retryWrites=true&w=majority"
    )

    try:
        print("MongoDB version is %s" %
                client.server_info()['version'])
    except pymongo.errors.OperationFailure as error:
        print(error)
    return client


def get_database(client,db_name):
    mydb = client[db_name]
    return mydb


def maybe_create_collection(db,collection_name):
    collections = db.list_collection_names(include_system_collections=False)
    if collection_name not in collections:
        collection = db.create_collection(collection_name)
    else:
        collection = db[collection_name]
    return collection

def get_list_of_files(path):
    files = []
    for entry in os.scandir(path):
        if entry.is_file():
            files.append(entry.path)
    return files

