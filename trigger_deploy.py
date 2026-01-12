"""
trigger_deploy.py
----------------
Script to run tests and, if successful, tag and push a new release version to git.

Modernized for Python 3.10+ and pytest >=6.0.
"""

import subprocess
import tomllib

import pytest

with open("pyproject.toml", "rb") as f:
    pyproject = tomllib.load(f)
version = pyproject["project"]["version"]

failures = pytest.main()
if failures == 0:
    subprocess.check_call(["git", "tag", version, "-m", f"version {version}"])
    subprocess.check_call(["git", "push", "--tags"])

# git push --tags origin pypi
# git tag 0.5.2 -m "version 0.5.2"
