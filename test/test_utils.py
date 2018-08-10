

def get_test_dir(key):

    import datetime
    import os

    time_stamp_string = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d%H%M%S")
    working_dir = os.path.join('/tmp', '{0}_{1}'.format(key, time_stamp_string))
    os.makedirs(working_dir)

    return working_dir

