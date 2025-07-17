import os
import mariadb
import sys
from sqlalchemy import create_engine

engine = None
def get_pd_remote_server():
    try:
        user = "jwbear"
        password = "fU8VYQSv9!GGH9vM"
        host = "cedar-mysql-vm.int.cedar.computecanada.ca"
        port = 3306
        database = "jwbear_DASH"

        global engine
        if engine is None:
            engine = create_engine(
                url="mysql+pymysql://{0}:{1}@{2}:{3}/{4}".format(
                    user, password, host, port, database
                ), pool_size=100, max_overflow=0, pool_timeout=120
            )
            return engine
        else:
            return engine

    except Exception as e:
        print(f"Error connecting to DB Platform: {e}")
        sys.exit(0)

def get_pd_remote():
    try:
        user = "jwbear"
        password = "fU8VYQSv9!GGH9vM"
        host = "127.0.0.1"
        port = 28901
        database = "jwbear_DASH"

        global engine
        if engine is None:
            engine = create_engine(
                url="mysql+pymysql://{0}:{1}@{2}:{3}/{4}".format(
                    user, password, host, port, database
                ), pool_size=100, max_overflow=0, pool_timeout=120
            )
            return engine
        else:
            return engine

    except Exception as e:
        print(f"Error connecting to DB Platform: {e}")
        sys.exit(0)

def get_local():
    try:
        user = "root"
        password = ""
        host = "127.0.0.1"
        port = 3306
        database = "jwbear_DASH"

        engine = create_engine(
            url="mysql+pymysql://{0}:{1}@{2}:{3}/{4}".format(
                user, password, host, port, database
            ) , pool_size=20, max_overflow=0, pool_timeout=60
        )

        return engine

    except Exception as e:
        print(f"Error connecting to DB Platform: {e}")
        sys.exit(1)
