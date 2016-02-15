import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks', help="""Input peaks file (Required).""")
    parser.add_argument('-b', '--bed', help="""Folder of bed files.""")
    parser.add_argument('-c', '--config', help="""Folder of config.py.""")
    args = parser.parse_args()
    return args
