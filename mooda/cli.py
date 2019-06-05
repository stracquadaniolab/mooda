import sys

import argh

from .commands import echo

def main():
    parser = argh.ArghParser()
    parser.add_commands([echo])
    parser.dispatch()

if __name__ == "__main__":
    sys.exit(main())
