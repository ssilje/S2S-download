import os
from paramiko import AutoAddPolicy, SSHClient
from netCDF4 import Dataset

# Requires ssh key

def make_dir(path):
    """
    Creates directory if it does not exist.

    args:
        path: str
    """
    if not os.path.exists(path):
        os.makedirs(path)


def pull(uni_path,filename,username,prnt=True):
    """
    Download /nird/projects/NS9853K/DATA to local.

    args:
        path:       str
        filename:   list of str
        username:   str

    returns:
        local_path
    """

    remote_path = '/nird/projects/NS9853K/DATA'+uni_path
    local_path  = './data/nird_copy'+uni_path

    # connect to nird
    ssh_client = SSHClient()
    ssh_client.load_system_host_keys()

    ssh_client.connect(
        'login-tos.nird.sigma2.no',
        username=username
    )

    make_dir(local_path)

    sftp_client = ssh_client.open_sftp()

    # copy file from nird to local
    for fn in filename:
        if not os.path.exists(local_path+fn):
            if prnt:
                print('downloading '+fn)
            sftp_client.get(remote_path+fn,local_path+fn)

    sftp_client.close()
    ssh_client.close()

    return local_path

# TODO: Consider making read() iterable
def read(uni_path,filename,username):
    """
    open and read /nird/projects/NS9853K/DATA at local.

    args:
        path:       str
        filename:   str
        username:   str

    returns:
        data

    - Only works with netCDF files
    """

    remote_path = '/nird/projects/NS9853K/DATA'+uni_path

    # connect to nird
    ssh_client = SSHClient()
    ssh_client.load_system_host_keys()

    ssh_client.connect(
        'login-tos.nird.sigma2.no',
        username=username
    )

    sftp_client = ssh_client.open_sftp()

    # read file from nird to local
    # for fn in filename:
    file = sftp_client.open(remote_path+filename)
    file.prefetch()
    raw_data = file.read()

    data = Dataset('temp.nc', memory=raw_data)

    sftp_client.close()
    ssh_client.close()

    return data
