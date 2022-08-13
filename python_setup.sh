#!/usr/bin/env bash
pip install virtualenv
virtualenv --python=/usr/local/bin/python3.8 ac-metro-schematization-framework/.venv38_pg22
source ac-metro-schematization-framework/.venv38_pg22/bin/activate
cd ac-metro-schematization-framework
pip install -r requirements_pg22.txt
pip install -e .