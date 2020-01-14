import argparse
import json


def main():

    parser = argparse.ArgumentParser(
        description='Removes jupytext metadata from Jupyter notebooks that establish pairing.')

    parser.add_argument('input',
        help='the Jupyter notebook to have pairing removed')

    parser.add_argument('--output',
        help='the output Jupyter notebook')

    args = parser.parse_args()

    if args.output == None:
        args.output = args.input

    with open(args.input, "r") as f:
        data = json.load(f)
        data["metadata"].pop("jupytext", None)

    with open(args.output, "w") as f:
        json.dump(data, f)


if __name__ == "__main__":
    main()

