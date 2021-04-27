"""
Perform tasks against a remote host.
"""

from .remote_config import (
    host,
    local_file_directory,
    remote_path,
    username,
)

from .client import RemoteClient
# from .files import fetch_local_files

def main():
    """Initialize remote host client and execute actions."""
    remote = RemoteClient(host, username, remote_path)
